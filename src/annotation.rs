use anyhow::{bail, Context, Result};
use bio::bio_types::strand::Strand;
use genebears::{AnnotateOptions, GeneBearError, GeneBears, Genome, Store, Variant};
use itertools::Itertools;
use log::{info, warn};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;
use std::thread;
use std::time::Duration;

use crate::graph::duck::variant_nodes;
use crate::graph::node::Node;
use crate::graph::transcript::Transcript;

/// GeneBe Hub databases with the hg38 scores, the column read from each and the field it fills.
pub(crate) const HUB_DATABASES: [(&str, &str, Field); 5] = [
    ("@genebe/revel:0.0.1", "score", |a| &mut a.revel_score),
    ("@genebe/alpha_missense:0.2.2", "am_pathogenicity", |a| {
        &mut a.alphamissense_score
    }),
    ("@genebe/spliceai:0.0.1", "max", |a| {
        &mut a.spliceai_max_score
    }),
    ("@genebe/gnomad_exomes4:0.0.4-4.1.0", "AF", |a| {
        &mut a.gnomad_exomes_af
    }),
    ("@genebe/gnomad_genomes4:0.0.4-4.1.0", "AF", |a| {
        &mut a.gnomad_genomes_af
    }),
];

type Field = fn(&mut VariantAnnotation) -> &mut Option<f64>;

/// Scores of a single variant, named like the fields of the GeneBe API.
#[derive(Debug, Clone, Default, Deserialize)]
pub(crate) struct VariantAnnotation {
    pub(crate) revel_score: Option<f64>,
    pub(crate) alphamissense_score: Option<f64>,
    pub(crate) spliceai_max_score: Option<f64>,
    pub(crate) gnomad_exomes_af: Option<f64>,
    pub(crate) gnomad_genomes_af: Option<f64>,
}

/// A GeneBe API response, with a warning if the variant could not be annotated.
#[derive(Deserialize)]
struct GeneBeResponse {
    warning: Option<String>,
    #[serde(flatten)]
    scores: VariantAnnotation,
}

#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct Annotation {
    pub(crate) revel_score: Option<f64>,
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
        annotations: &HashMap<String, VariantAnnotation>,
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
            spliceai_score: probabilistic_or(all.iter().map(|(_, a)| a.spliceai_max_score)),
            alphamissense_score: probabilistic_or(translated.iter().map(|a| a.alphamissense_score)),
            gnomad_frequencies: all
                .iter()
                .map(|(token, annotation)| (token.clone(), gnomad_af(annotation)))
                .collect(),
        }
    }
}

/// Where the scores of the variants come from.
pub(crate) enum Scores {
    /// The local GeneBe Hub databases, for hg38.
    Hub(Store),
    Api(Box<GeneBears>, Genome),
}

impl Scores {
    /// The GeneBe Hub databases in the store for hg38 if they are installed, otherwise the API.
    pub(crate) fn new(genome: Genome, store: Store, client: GeneBears) -> Result<Self> {
        let installed = HUB_DATABASES
            .iter()
            .filter(|(id, _, _)| id.parse().is_ok_and(|id| store.get(&id).is_ok()))
            .count();
        match genome {
            Genome::Hg38 if installed == HUB_DATABASES.len() => Ok(Scores::Hub(store)),
            Genome::Hg38 if installed > 0 => bail!(
                "GeneBe Hub databases in {} are incomplete, rerun `haplodon cache`",
                store.root().display()
            ),
            _ => Ok(Scores::Api(Box::new(client), genome)),
        }
    }
}

/// Annotates the coding variants of the transcripts once per target, by HGVS token.
pub(crate) fn annotate_coding_variants(
    graph: &Path,
    transcripts: &[Transcript],
    scores: &Scores,
) -> Result<HashMap<String, HashMap<String, VariantAnnotation>>> {
    transcripts
        .iter()
        .into_group_map_by(|transcript| transcript.target.as_str())
        .into_iter()
        .map(|(target, transcripts)| {
            let nodes = coding_variants(variant_nodes(graph, target)?, &transcripts);
            let variants = nodes
                .iter()
                .map(|node| {
                    Variant::new(
                        target,
                        node.pos as u64 + 1, // GeneBe expects 1-based positions
                        node.reference_allele.as_str(),
                        node.alternative_allele.as_str(),
                    )
                })
                .collect_vec();
            info!("Annotating {} coding variants on {target}", variants.len());
            let annotations = match scores {
                Scores::Hub(store) => lookup(&variants, store)?,
                Scores::Api(client, genome) => fetch(&variants, *genome, client)
                    .with_context(|| format!("GeneBe failed to annotate variants on {target}"))?,
            };
            let annotations = nodes
                .iter()
                .zip(annotations)
                .filter_map(|(node, annotation)| Some((node.hgvsg_token(), annotation?)))
                .collect();
            Ok((target.to_string(), annotations))
        })
        .collect()
}

/// The variants, sorted by position, that lie in a coding sequence of the transcripts.
fn coding_variants(variants: Vec<Node>, transcripts: &[&Transcript]) -> Vec<Node> {
    let mut coding = vec![false; variants.len()];
    for cds in transcripts.iter().flat_map(|t| &t.coding_sequences) {
        let start = variants.partition_point(|node| node.pos < cds.start as i64);
        let end = variants.partition_point(|node| node.pos <= cds.end as i64);
        coding[start..end].fill(true);
    }
    variants
        .into_iter()
        .zip(coding)
        .filter_map(|(node, coding)| coding.then_some(node))
        .collect()
}

/// Scores from the local GeneBe Hub databases, the highest of several rows, e.g. per gene.
fn lookup(variants: &[Variant], store: &Store) -> Result<Vec<Option<VariantAnnotation>>> {
    let mut annotations = vec![VariantAnnotation::default(); variants.len()];
    for (id, column, field) in HUB_DATABASES {
        let rows: Vec<Vec<HashMap<String, Option<f64>>>> =
            store
                .get(&id.parse()?)?
                .lookup(variants, Genome::Hg38, &[column])?;
        for (annotation, rows) in annotations.iter_mut().zip(rows) {
            *field(annotation) = rows.iter().filter_map(|row| row[column]).reduce(f64::max);
        }
    }
    Ok(annotations.into_iter().map(Some).collect())
}

/// Scores from the GeneBe API, without the ones GeneBe warns about, e.g. for a wrong reference
/// allele. Server and network errors are retried per batch.
fn fetch(
    variants: &[Variant],
    genome: Genome,
    client: &GeneBears,
) -> Result<Vec<Option<VariantAnnotation>>> {
    let rt = tokio::runtime::Runtime::new()?;
    let mut responses = Vec::with_capacity(variants.len());
    // As many variants as GeneBe takes per request.
    for batch in variants.chunks(1000) {
        // Seconds to wait before asking GeneBe again after a server or network error.
        let mut delays = [10, 30, 60].into_iter();
        loop {
            let options = AnnotateOptions::default();
            match (
                rt.block_on(client.annotate::<GeneBeResponse>(batch, genome, options)),
                delays.next(),
            ) {
                (Ok(batch), _) => break responses.extend(batch),
                (
                    Err(GeneBearError::ApiServerError { .. } | GeneBearError::Http(_)),
                    Some(delay),
                ) => {
                    warn!("GeneBe failed, retrying in {delay} s");
                    thread::sleep(Duration::from_secs(delay));
                }
                (Err(e), _) => return Err(e.into()),
            }
        }
    }
    Ok(variants
        .iter()
        .zip(responses)
        .map(|(variant, response)| match response.warning {
            Some(warning) => {
                warn!("Ignoring GeneBe annotation for {variant:?}: {warning}");
                None
            }
            None => Some(response.scores),
        })
        .collect())
}

fn gnomad_af(variant: &VariantAnnotation) -> Option<f64> {
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

    fn annotated(revel_score: f64) -> VariantAnnotation {
        VariantAnnotation {
            revel_score: Some(revel_score),
            ..Default::default()
        }
    }

    #[test]
    fn fetch_skips_variants_with_wrong_reference() {
        let client = GeneBears::new(ClientConfig::default()).unwrap();
        let variants = [
            Variant::new("1", 11_796_321, "G", "T"),
            Variant::new("1", 11_796_321, "A", "T"),
        ];

        let annotations = fetch(&variants, Genome::Hg38, &client).unwrap();

        let annotation = annotations[0].as_ref().unwrap();
        assert!(annotation.revel_score.is_some() && annotation.alphamissense_score.is_some());
        assert!(annotations[1].is_none());
    }

    #[test]
    fn lookup_takes_the_highest_row_per_database() {
        let dir = tempfile::tempdir().unwrap();
        let db = duckdb::Connection::open_in_memory().unwrap();
        for (id, column, _) in HUB_DATABASES {
            let id: genebears::DatabaseId = id.parse().unwrap();
            let version = id.version.unwrap();
            let path = dir.path().join(&id.owner).join(&id.name).join(&version);
            let rows = "VALUES ('1', 99, 1, 'T', 0.25), ('1', 99, 1, 'T', 0.5)";
            std::fs::create_dir_all(&path).unwrap();
            db.execute_batch(&format!(
                "COPY (SELECT * FROM ({rows}) t(_seq, _pos, _del, _ins, \"{column}\")) \
                 TO '{}' (FORMAT PARQUET, PARTITION_BY (_seq))",
                path.join("parquet").display()
            ))
            .unwrap();
            let description = format!("columns = [{{ name = \"{column}\" }}]");
            std::fs::write(path.join("description.toml"), description).unwrap();
        }
        let variants = [
            Variant::new("chr1", 100, "A", "T"),
            Variant::new("1", 100, "A", "G"),
        ];

        let annotations = lookup(&variants, &Store::new(dir.path())).unwrap();

        let found = annotations[0].as_ref().unwrap();
        assert_eq!(found.revel_score, Some(0.5));
        assert_eq!(found.gnomad_genomes_af, Some(0.5));
        assert_eq!(annotations[1].as_ref().unwrap().revel_score, None);
    }

    #[test]
    fn coding_variants_of_overlapping_transcripts() {
        let transcript = |coding| Transcript::new("t".into(), "1".into(), Strand::Forward, coding);
        let first = transcript(vec![Cds::new(10, 20, 0), Cds::new(40, 50, 0)]);
        let second = transcript(vec![Cds::new(15, 30, 0)]);
        let variants = [5, 10, 18, 30, 35, 50, 51].map(|pos| variant(pos, "A", "T"));

        let coding = coding_variants(variants.to_vec(), &[&first, &second]);

        assert_eq!(
            coding.iter().map(|node| node.pos).collect_vec(),
            [10, 18, 30, 50]
        );
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
