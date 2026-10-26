use anyhow::Result;
use bio::bio_types::strand::Strand;
use genebears::{AnnotateOptions, AnnotatedVariant, GeneBearError, GeneBears, Genome, Variant};
use itertools::Itertools;
use log::warn;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use crate::graph::node::Node;
use crate::graph::transcript::Transcript;

#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct Annotation {
    pub(crate) revel_score: Option<f64>,
    pub(crate) acmg_score: Option<f64>,
    pub(crate) spliceai_score: Option<f64>,
    pub(crate) alphamissense_score: Option<f64>,
    #[serde(default)]
    pub(crate) gnomad_frequencies: HashMap<String, Option<f64>>,
}

impl Annotation {
    /// Aggregates the annotations of the variants of a haplotype. Missense scores only count
    /// variants in the original frame and before the stop codon at `stop`.
    pub(crate) fn from_haplotype(
        haplotype: &[Node],
        transcript: &Transcript,
        annotations: &HashMap<String, AnnotatedVariant>,
        stop: Option<usize>,
    ) -> Self {
        let mut variants: Vec<_> = haplotype
            .iter()
            .filter(|n| n.node_type.is_variant())
            .collect();
        if transcript.strand == Strand::Reverse {
            variants.reverse();
        }
        let mut all = Vec::new();
        let mut translated = Vec::new();
        let mut shift = 0;
        for node in variants {
            // Indels before the variant shift its codon in the altered protein.
            let before_stop = match (stop, transcript.position_in_transcript(node.pos as usize)) {
                (Some(stop), Ok(position)) => (position as i64 + shift) / 3 < stop as i64,
                _ => true,
            };
            let in_frame = shift % 3 == 0 && before_stop;
            shift += node.frameshift();
            if let Some(annotation) = annotations.get(&node.hgvsg_token()) {
                all.push((node.hgvsg_token(), annotation));
                if in_frame {
                    translated.push(annotation);
                }
            }
        }

        Self {
            revel_score: probabilistic_or(translated.iter().map(|a| a.revel_score)),
            acmg_score: max_score(all.iter().map(|(_, a)| a.acmg_score)),
            spliceai_score: probabilistic_or(all.iter().map(|(_, a)| a.spliceai_max_score)),
            alphamissense_score: probabilistic_or(translated.iter().map(|a| a.alphamissense_score)),
            gnomad_frequencies: all
                .iter()
                .map(|(token, annotation)| (token.clone(), gnomad_af(annotation)))
                .collect(),
        }
    }
}

/// GeneBe annotations by HGVS token, skipping variants GeneBe warns about.
pub(crate) fn fetch_annotations(
    variants: &[&Node],
    target: &str,
    genome_build: Genome,
    genebe_client: &Arc<Mutex<GeneBears>>,
) -> Result<HashMap<String, AnnotatedVariant>> {
    if variants.is_empty() {
        return Ok(HashMap::new());
    }
    let requests = variants
        .iter()
        .map(|node| {
            Variant::new(
                target.to_string(),
                node.pos as u64 + 1, // GeneBe expects 1-based positions
                node.reference_allele.to_string(),
                node.alternative_allele.to_string(),
            )
        })
        .collect_vec();

    let results = {
        let client = genebe_client.lock().unwrap();
        let rt = tokio::runtime::Runtime::new()?;
        rt.block_on(client.annotate_variants_chunked(
            &requests,
            genome_build,
            AnnotateOptions::default(),
        ))
    };

    let results = match results {
        Ok(r) => r,
        Err(GeneBearError::ApiServerError { message, .. }) => {
            warn!(
                "GeneBe failed to annotate variants for {}: {}. Scores will be absent.",
                target, message
            );
            vec![]
        }
        Err(e) => return Err(anyhow::Error::from(e)),
    };

    Ok(variants
        .iter()
        .zip(results)
        .filter(|(node, result)| {
            if let Some(message) = &result.warning {
                warn!(
                    "Ignoring GeneBe annotation for {} on {}: {}",
                    node.hgvsg_token(),
                    target,
                    message
                );
            }
            result.warning.is_none()
        })
        .map(|(node, result)| (node.hgvsg_token(), result))
        .collect())
}

fn gnomad_af(variant: &AnnotatedVariant) -> Option<f64> {
    max_score([variant.gnomad_exomes_af, variant.gnomad_genomes_af].into_iter())
}

fn max_score(iter: impl Iterator<Item = Option<f64>>) -> Option<f64> {
    iter.flatten().reduce(f64::max)
}

fn probabilistic_or(iter: impl Iterator<Item = Option<f64>>) -> Option<f64> {
    iter.flatten()
        .map(|score| 1.0 - score)
        .reduce(|acc, complement| acc * complement)
        .map(|product| 1.0 - product)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::node::NodeType;
    use crate::graph::paths::Cds;
    use genebears::ClientConfig;
    use std::collections::HashMap;

    #[test]
    fn probabilistic_or_of_a_single_score_is_unchanged() {
        let result = probabilistic_or([Some(0.6)].into_iter()).unwrap();
        assert!((result - 0.6).abs() < 1e-9);
    }

    #[test]
    fn probabilistic_or_combines_as_complement_of_product() {
        let result = probabilistic_or([Some(0.5), Some(0.5)].into_iter()).unwrap();
        assert!((result - 0.75).abs() < 1e-9);
    }

    #[test]
    fn probabilistic_or_ignores_absent_scores() {
        let result = probabilistic_or([Some(0.5), None, Some(0.5)].into_iter()).unwrap();
        assert!((result - 0.75).abs() < 1e-9);
    }

    #[test]
    fn probabilistic_or_without_scores_is_none() {
        assert_eq!(probabilistic_or([None, None].into_iter()), None);
    }

    fn variant(pos: i64, reference_allele: &str, alternative_allele: &str) -> Node {
        Node::new(
            NodeType::Variant,
            pos,
            reference_allele.to_string(),
            alternative_allele.to_string(),
        )
    }

    fn annotated(revel_score: f64) -> AnnotatedVariant {
        AnnotatedVariant {
            revel_score: Some(revel_score),
            ..Default::default()
        }
    }

    #[test]
    fn fetch_annotations_returns_scores() {
        let client = Arc::new(Mutex::new(GeneBears::new(ClientConfig::default()).unwrap()));
        let node = variant(11_796_320, "G", "T");

        let annotations = fetch_annotations(&[&node], "1", Genome::Hg38, &client).unwrap();

        let annotation = &annotations["11796321G>T"];
        assert!(annotation.alphamissense_score.is_some());
        assert!(annotation.revel_score.is_some());
    }

    #[test]
    fn fetch_annotations_ignores_variants_with_mismatching_reference() {
        let client = Arc::new(Mutex::new(GeneBears::new(ClientConfig::default()).unwrap()));
        let node = variant(11_796_320, "A", "T");

        let annotations = fetch_annotations(&[&node], "1", Genome::Hg38, &client).unwrap();

        assert!(annotations.is_empty());
    }

    #[test]
    fn from_haplotype_on_reverse_strand() {
        let transcript =
            Transcript::new("test".to_string(), "1".to_string(), Strand::Reverse, vec![]);
        // In genomic order the deletion comes first, but the transcript starts at the SNV.
        let haplotype = vec![variant(10, "AC", "A"), variant(20, "G", "T")];
        let annotations = HashMap::from([
            ("11AC>A".to_string(), annotated(0.5)),
            ("21G>T".to_string(), annotated(0.8)),
        ]);

        let annotation = Annotation::from_haplotype(&haplotype, &transcript, &annotations, None);

        assert!((annotation.revel_score.unwrap() - 0.9).abs() < 1e-9);
        assert_eq!(annotation.gnomad_frequencies.len(), 2);
    }

    #[test]
    fn missense_scores_after_restored_frame() {
        let transcript =
            Transcript::new("test".to_string(), "1".to_string(), Strand::Forward, vec![]);
        let haplotype = vec![
            variant(10, "AC", "A"),
            variant(15, "G", "T"),
            variant(20, "A", "AC"),
            variant(30, "G", "T"),
        ];
        let annotations = HashMap::from([
            ("11AC>A".to_string(), annotated(0.5)),
            ("16G>T".to_string(), annotated(0.3)),
            ("21A>AC".to_string(), annotated(0.4)),
            ("31G>T".to_string(), annotated(0.8)),
        ]);

        let annotation = Annotation::from_haplotype(&haplotype, &transcript, &annotations, None);

        // Only the deletion and the SNV after the restoring insertion are in frame.
        assert!((annotation.revel_score.unwrap() - 0.9).abs() < 1e-9);
        assert_eq!(annotation.gnomad_frequencies.len(), 4);
    }

    #[test]
    fn missense_scores_skip_variants_after_stop() {
        let transcript = Transcript::new(
            "test".to_string(),
            "1".to_string(),
            Strand::Forward,
            vec![Cds::new(0, 99, 0)],
        );
        // Missense variants in codon 1 and 13.
        let haplotype = vec![variant(4, "G", "T"), variant(40, "G", "T")];
        let annotations = HashMap::from([
            ("5G>T".to_string(), annotated(0.5)),
            ("41G>T".to_string(), annotated(0.8)),
        ]);

        let stopped = Annotation::from_haplotype(&haplotype, &transcript, &annotations, Some(3));

        assert_eq!(stopped.revel_score, Some(0.5));
        assert_eq!(stopped.gnomad_frequencies.len(), 2);
    }
}
