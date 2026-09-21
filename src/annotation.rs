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
    pub(crate) fn from_haplotype(
        haplotype: &[Node],
        transcript: &Transcript,
        genome_build: Genome,
        genebe_client: &Arc<Mutex<GeneBears>>,
    ) -> Result<Self> {
        let mut variants: Vec<_> = haplotype
            .iter()
            .filter(|n| n.node_type.is_variant())
            .collect();
        if transcript.strand == Strand::Reverse {
            variants.reverse();
        }
        let in_phase: Vec<&Node> = variants
            .into_iter()
            .scan(0i64, |running_fs, node| {
                let fs = node.frameshift();
                let in_phase = *running_fs % 3 == 0;
                *running_fs += fs;
                Some((node, in_phase))
            })
            .filter_map(|(node, in_phase)| in_phase.then_some(node))
            .collect();
        let variants = in_phase
            .iter()
            .map(|node| {
                Variant::new(
                    transcript.target.to_string(),
                    node.pos as u64 + 1, // GeneBe expects 1-based positions
                    node.reference_allele.to_string(),
                    node.alternative_allele.to_string(),
                )
            })
            .collect_vec();

        let results = {
            let client = genebe_client.lock().unwrap();
            let rt = tokio::runtime::Runtime::new()?;
            rt.block_on(client.annotate_variants(
                &variants,
                genome_build,
                AnnotateOptions::default(),
            ))
        };

        let results = match results {
            Ok(r) => r,
            Err(GeneBearError::ApiServerError { message, .. }) => {
                warn!(
                    "GeneBe failed to annotate variants for {}: {}. Scores will be absent.",
                    transcript.target, message
                );
                vec![]
            }
            Err(e) => return Err(anyhow::Error::from(e)),
        };

        let frequency_by_variant: HashMap<(u64, String, String), Option<f64>> = results
            .iter()
            .filter_map(|r| {
                Some((
                    (r.pos?, r.ref_allele.clone()?, r.alt.clone()?),
                    gnomad_af(r),
                ))
            })
            .collect();
        let gnomad_frequencies = in_phase
            .iter()
            .map(|node| {
                let frequency = frequency_by_variant
                    .get(&(
                        node.pos as u64 + 1,
                        node.reference_allele.clone(),
                        node.alternative_allele.clone(),
                    ))
                    .copied()
                    .flatten();
                (node.hgvsg_token(), frequency)
            })
            .collect();

        Ok(Self {
            revel_score: probabilistic_or(results.iter().map(|r| r.revel_score)),
            acmg_score: max_score(results.iter().map(|r| r.acmg_score)),
            spliceai_score: probabilistic_or(results.iter().map(|r| r.spliceai_max_score)),
            alphamissense_score: probabilistic_or(results.iter().map(|r| r.alphamissense_score)),
            gnomad_frequencies,
        })
    }
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
    use crate::graph::EventProbs;
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

    #[test]
    fn from_haplotype_annotates_a_single_variant() {
        let client = Arc::new(Mutex::new(GeneBears::new(ClientConfig::default()).unwrap()));
        let transcript =
            Transcript::new("test".to_string(), "1".to_string(), Strand::Forward, vec![]);
        let haplotype = vec![Node {
            node_type: NodeType::Variant,
            reference_allele: "G".to_string(),
            alternative_allele: "T".to_string(),
            vaf: HashMap::new(),
            probs: EventProbs(HashMap::new()),
            pos: 11_796_320,
            index: 0,
        }];

        let annotation =
            Annotation::from_haplotype(&haplotype, &transcript, Genome::Hg38, &client).unwrap();

        assert!(annotation.alphamissense_score.is_some());
        assert!(annotation.revel_score.is_some());
    }
}
