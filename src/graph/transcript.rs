use crate::annotation::Annotation;
use crate::cli::Interval;
use crate::graph::duck::{feature_graph, variants_on_graph};
use crate::graph::node::{Node, NodeType};
use crate::graph::paths::{Cds, HaplotypePath};
use crate::graph::peptide::Peptide;
use crate::graph::score::EffectScore;
use crate::graph::score::HaplotypeMetric;
use crate::graph::score::HaplotypeScore;
use crate::graph::{EventProbs, VariantGraph};
use crate::translation::amino_acids::Protein;
use crate::translation::distance::DistanceMetric;
use anyhow::{bail, Result};
use bio::alphabets::dna;
use bio::bio_types::strand::Strand;
use bio::io::gff::{self};
use bio::stats::LogProb;
use genebears::GeneBears;
use genebears::Genome;
use itertools::Itertools;
use log::info;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub(crate) struct Transcript {
    pub(crate) feature: String,
    pub(crate) target: String,
    pub(crate) strand: Strand,
    pub(crate) coding_sequences: Vec<Cds>,
}

type Haplotypes = Vec<(Vec<Node>, Vec<HashMap<String, u32>>)>;

impl Transcript {
    pub(crate) fn new(
        feature: String,
        target: String,
        strand: Strand,
        coding_sequences: Vec<Cds>,
    ) -> Transcript {
        Transcript {
            feature,
            target,
            strand,
            coding_sequences,
        }
    }

    pub(crate) fn name(&self) -> String {
        format!("{}:{}", self.target, self.feature)
    }

    pub(crate) fn start(&self) -> Result<u64> {
        self.coding_sequences
            .iter()
            .map(|cds| cds.start)
            .min()
            .ok_or_else(|| anyhow::anyhow!("No CDS found for transcript {}", self.name()))
    }

    pub(crate) fn end(&self) -> Result<u64> {
        self.coding_sequences
            .iter()
            .map(|cds| cds.end)
            .max()
            .ok_or_else(|| anyhow::anyhow!("No CDS found for transcript {}", self.name()))
    }

    /// Returns an iterator over the coding sequences of the transcript
    /// The iterator is sorted by start position in ascending order
    /// If the strand is reverse, the iterator is reversed
    pub(crate) fn cds(&self) -> Box<dyn Iterator<Item = &Cds> + '_> {
        match self.strand {
            Strand::Reverse => Box::new(
                self.coding_sequences
                    .iter()
                    .sorted_by_key(|cds| cds.start)
                    .rev(),
            ),
            _ => Box::new(self.coding_sequences.iter().sorted_by_key(|cds| cds.start)),
        }
    }

    /// Returns true if the reference coding sequence consists solely of A, C, G
    /// and T bases and can therefore be translated.
    pub(crate) fn is_translatable(&self, reference: &HashMap<String, Vec<u8>>) -> bool {
        self.untranslatable_reference_base(reference).is_none()
    }

    /// Returns the first coding reference base that is not A, C, G or T, if any.
    /// Such bases (`N` or IUPAC ambiguity codes) cannot be translated.
    pub(crate) fn untranslatable_reference_base(
        &self,
        reference: &HashMap<String, Vec<u8>>,
    ) -> Option<u8> {
        let target = reference.get(&self.target)?;
        self.cds()
            .flat_map(|cds| &target[cds.start as usize..=cds.end as usize])
            .copied()
            .find(|&base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
    }

    /// Returns the 0-based position of the variant from the transcript start,
    /// where "start" means the 5' end of the coding sequence — i.e. the highest
    /// genomic coordinate for reverse-strand transcripts.
    pub(crate) fn position_in_transcript(&self, pos: usize) -> Result<usize> {
        let mut offset = 0;
        for cds in self.cds() {
            if ((cds.start as usize)..=(cds.end as usize)).contains(&pos) {
                let cds_offset = if self.strand == Strand::Reverse {
                    cds.end as usize - pos
                } else {
                    pos - cds.start as usize
                };
                return Ok(offset + cds_offset);
            }
            offset += cds.length();
        }
        bail!(
            "Position {} not found in CDS of transcript {}",
            pos,
            self.name()
        );
    }

    fn paths(&self, graph: &VariantGraph) -> Result<Vec<HaplotypePath>> {
        match self.strand {
            Strand::Forward => Ok(graph.paths()),
            Strand::Reverse => Ok(graph.reverse_paths()),
            Strand::Unknown => Err(anyhow::anyhow!(
                "Strand is unknown for transcript {}",
                self.name()
            )),
        }
    }

    fn haplotypes(
        &self,
        graph: &PathBuf,
        max_haplotypes: usize,
    ) -> Result<(Haplotypes, HashSet<String>)> {
        let graph = match feature_graph(
            graph.to_owned(),
            self.target.clone(),
            self.start()?,
            self.end()?,
        ) {
            Ok(g) => g,
            Err(_) => return Ok((vec![], HashSet::new())),
        };

        let samples: HashSet<String> = graph
            .graph
            .node_indices()
            .flat_map(|i| graph.graph[i].vaf.keys().cloned())
            .collect();

        let cds_intervals: Vec<(u64, u64)> = self.cds().map(|c| (c.start, c.end)).collect();
        let in_cds = |pos: i64| {
            cds_intervals
                .iter()
                .any(|(s, e)| pos >= *s as i64 && pos <= *e as i64)
        };

        let variant_count = graph
            .graph
            .node_indices()
            .filter(|&i| graph.graph[i].node_type.is_variant())
            .filter(|&i| in_cds(graph.graph[i].pos))
            .count();

        info!(
            "Calculating paths for transcript {} with {} variant nodes",
            self.name(),
            variant_count
        );

        let raw_paths = match self.strand {
            Strand::Forward => graph.top_k_paths(max_haplotypes),
            Strand::Reverse => graph.reverse_top_k_paths(max_haplotypes),
            Strand::Unknown => bail!("Strand is unknown for transcript {}", self.name()),
        };

        let haplotypes = raw_paths
            .iter()
            .map(|p| {
                let nodes: Vec<Node> =
                    p.0.iter()
                        .map(|v| graph.graph.node_weight(*v).unwrap().to_owned())
                        .collect();
                let edge_reads = graph.edge_reads(&p.0);

                let mut filtered_nodes: Vec<Node> = Vec::new();
                let mut filtered_edges: Vec<HashMap<String, u32>> = Vec::new();
                let mut last_cds_idx: Option<usize> = None;

                for (i, node) in nodes.iter().enumerate() {
                    if !in_cds(node.pos) {
                        continue;
                    }
                    if let Some(prev_idx) = last_cds_idx {
                        let span = &edge_reads[prev_idx..i];
                        let merged = if span.is_empty() {
                            HashMap::new()
                        } else {
                            let samples: std::collections::HashSet<&String> =
                                span.iter().flat_map(|e| e.keys()).collect();
                            samples
                                .into_iter()
                                .map(|s| {
                                    let min = span
                                        .iter()
                                        .map(|e| *e.get(s).unwrap_or(&0))
                                        .min()
                                        .unwrap_or(0);
                                    (s.clone(), min)
                                })
                                .collect()
                        };
                        filtered_edges.push(merged);
                    }
                    filtered_nodes.push(node.clone());
                    last_cds_idx = Some(i);
                }

                (filtered_nodes, filtered_edges)
            })
            .unique_by(|(nodes, _)| {
                nodes
                    .iter()
                    .map(|n| {
                        (
                            n.pos,
                            n.node_type.clone(),
                            n.reference_allele.clone(),
                            n.alternative_allele.clone(),
                        )
                    })
                    .collect_vec()
            })
            .collect();

        Ok((haplotypes, samples))
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn scores(
        &self,
        graph: &PathBuf,
        reference: &HashMap<String, Vec<u8>>,
        haplotype_metric: HaplotypeMetric,
        max_haplotypes: usize,
        distance_metric: DistanceMetric,
        genome_build: Genome,
        genebe_client: &Arc<Mutex<GeneBears>>,
    ) -> Result<Vec<HaplotypeScore>> {
        let (haplotypes, samples) = self.haplotypes(graph, max_haplotypes)?;
        let mut scores = Vec::with_capacity(haplotypes.len());
        let original_protein = Protein::from_transcript(reference, self)?;
        for (haplotype, supporting_reads) in haplotypes {
            let realign = haplotype.iter().map(|node| node.frameshift()).sum::<i64>() != 0;
            let effect_score = EffectScore::from_haplotype(
                reference,
                self,
                &haplotype,
                original_protein.clone(),
                distance_metric,
                realign,
            )?;
            let frequency = haplotype_metric.calculate(&haplotype, &samples);
            let annotation =
                Annotation::from_haplotype(&haplotype, self, genome_build, genebe_client)?;
            scores.push((effect_score, frequency, supporting_reads, annotation));
        }
        Ok(scores)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn peptides(
        &self,
        graph: &PathBuf,
        reference: &HashMap<String, Vec<u8>>,
        interval: Interval,
        sample: &str,
        events: &[String],
        min_event_prob: LogProb,
        background_events: &[String],
        max_background_event_prob: LogProb,
    ) -> Result<Vec<Peptide>> {
        let rnas = self.rna(graph, reference, sample)?;
        let mut peptides = Vec::new();
        for i in interval {
            for rna in &rnas {
                let p = rna.peptides(i)?;
                for peptide in p {
                    peptides.push(peptide);
                }
            }
        }
        peptides = peptides
            .iter()
            .filter(|peptide| {
                peptide.prob(events).unwrap() >= min_event_prob
                    && peptide.prob(background_events).unwrap() <= max_background_event_prob
            })
            .cloned()
            .collect_vec();
        Ok(peptides)
    }

    // Generate RNA sequences for the transcript, create all possible paths by combining the CDSs paths and then apply the variants to the RNA sequences
    pub(crate) fn rna(
        &self,
        graph: &PathBuf,
        reference: &HashMap<String, Vec<u8>>,
        sample: &str,
    ) -> Result<Vec<RnaPath>> {
        let reference = reference.get(&self.target).unwrap();
        let mut rna = Vec::new();
        for cds in self.cds() {
            let mut cds_rna = Vec::new();
            if let Ok(graph) = feature_graph(
                graph.to_owned(),
                self.target.to_string(),
                cds.start,
                cds.end,
            ) {
                info!(
                    "Subgraph for CDS ({}-{}) of transcript {} has {} nodes",
                    cds.start,
                    cds.end,
                    self.name(),
                    graph.graph.node_count()
                );
                let mut sequence = reference[cds.start as usize..=cds.end as usize].to_vec();
                if self.strand == Strand::Reverse {
                    sequence = dna::revcomp(&sequence);
                }
                sequence = sequence[cds.phase as usize..].to_vec();
                for path in self.paths(&graph)? {
                    let mut variants = HashMap::new();
                    let mut frameshift = 0;
                    let mut path_sequence = sequence.clone();
                    for node_index in path.0.iter() {
                        let node = graph.graph.node_weight(*node_index).unwrap();
                        if node.node_type == NodeType::Variant {
                            let position_in_cds = match self.strand {
                                Strand::Forward => {
                                    (node.pos - cds.start as i64 + frameshift + cds.phase as i64)
                                        as usize
                                }
                                Strand::Reverse => {
                                    (cds.end as i64 - node.pos + frameshift - cds.phase as i64)
                                        as usize
                                }
                                Strand::Unknown => {
                                    return Err(anyhow::anyhow!("Strand is unknown"))
                                }
                            };
                            variants.insert(
                                position_in_cds,
                                (
                                    graph.graph.node_weight(*node_index).unwrap().probs.clone(),
                                    node.frameshift(),
                                    *node.vaf.get(sample).unwrap(),
                                ),
                            );
                            match self.strand {
                                Strand::Forward => {
                                    path_sequence.splice(
                                        position_in_cds..position_in_cds + 1,
                                        node.alternative_allele.bytes(),
                                    );
                                }
                                Strand::Reverse => {
                                    path_sequence.splice(
                                        position_in_cds..position_in_cds + 1,
                                        String::from_utf8_lossy(
                                            dna::revcomp(node.alternative_allele.as_bytes())
                                                .as_slice(),
                                        )
                                        .bytes(),
                                    );
                                }
                                Strand::Unknown => unreachable!(),
                            }
                            frameshift += node.frameshift();
                        }
                    }
                    let rna_path = RnaPath::new(path_sequence, variants, self.clone());
                    cds_rna.push(rna_path);
                }
            }
            if rna.is_empty() {
                rna = cds_rna;
            } else if cds_rna.is_empty() {
                continue;
            } else {
                let mut new_rna = Vec::new();
                for rna_path in &rna {
                    for cds_rna_path in &cds_rna {
                        new_rna.push(rna_path.merge(cds_rna_path.clone())?);
                    }
                }
                rna = new_rna;
            }
        }
        Ok(rna)
    }
}

/// Represents an RNA sequence for a transcript and holds all variants on it
#[derive(Clone, Debug)]
pub(crate) struct RnaPath {
    pub(crate) sequence: Vec<u8>,
    pub(crate) variants: HashMap<usize, (EventProbs, i64, f32)>,
    pub(crate) transcript: Transcript,
}

impl RnaPath {
    pub(crate) fn new(
        sequence: Vec<u8>,
        variants: HashMap<usize, (EventProbs, i64, f32)>,
        transcript: Transcript,
    ) -> RnaPath {
        RnaPath {
            sequence,
            variants,
            transcript,
        }
    }

    pub(crate) fn merge(&self, other: RnaPath) -> Result<RnaPath> {
        assert_eq!(
            self.transcript, other.transcript,
            "Transcripts must be the same to merge both paths"
        );
        let mut sequence = self.sequence.clone();
        let mut variants = self.variants.clone();
        let offset = self.sequence.len();
        sequence.extend(other.sequence);
        for (position, probs) in other.variants {
            variants.insert(position + offset, probs);
        }
        Ok(RnaPath {
            sequence,
            variants,
            transcript: self.transcript.clone(),
        })
    }

    pub(crate) fn peptides(&self, length: u32) -> Result<Vec<Peptide>> {
        let mut peptides = Vec::new();
        for (i, p) in self
            .sequence
            .windows(length as usize * 3)
            .step_by(3)
            .enumerate()
        {
            let interval = (i * 3, i * 3 + length as usize * 3);
            // We consider a variant to influence the peptide if it is in the interval or if it is before the interval and causes a frameshift not equal to 0
            let probs_variants_in_interval = self
                .variants
                .iter()
                .filter(|(pos, (_, fs, _))| {
                    (*pos >= &interval.0 && *pos < &interval.1) || (*pos < &interval.0 && fs != &0)
                })
                .map(|(_, (probs, _, af))| (probs, *af))
                .collect_vec();
            // We skip peptides that have no variants in the interval or a frameshift upstream as they are not of any interest
            if probs_variants_in_interval.is_empty() {
                continue;
            }
            let mut probs = HashMap::new();
            for (variant, _) in &probs_variants_in_interval {
                for (event, prob) in variant.0.iter() {
                    // Add event or multiply probability if it already exists. We are summing because we are in log space
                    match probs.get_mut(event) {
                        Some(p) => *p += *prob,
                        None => {
                            probs.insert(event.clone(), *prob);
                        }
                    }
                }
            }
            let afs = probs_variants_in_interval
                .iter()
                .map(|(_, af)| *af)
                .collect();
            let peptide =
                Peptide::from_rna(p.to_vec(), EventProbs(probs), afs, self.transcript.clone())?;
            peptides.push(peptide);
        }
        Ok(peptides)
    }
}

pub(crate) fn transcripts(gff_file: &PathBuf, graph: &PathBuf) -> Result<Vec<Transcript>> {
    let variants = variants_on_graph(graph)?;
    let mut feature_reader = gff::Reader::from_file(gff_file, gff::GffType::GFF3)?;
    let mut transcripts = HashMap::new();
    let mut variant_coverage = HashMap::new();
    for record in feature_reader
        .records()
        .filter_map(Result::ok)
        .filter(|record| record.feature_type() == "CDS")
    {
        let ensp = record.attributes().get("ID").ok_or_else(|| {
            anyhow::anyhow!("No ID found for CDS in sequence {}", record.seqname())
        })?;
        let target = record.seqname().to_string();
        let start = *record.start() - 1;
        let end = *record.end() - 1;
        let phase = record.phase().clone().try_into().unwrap();
        let strand = record.strand().ok_or_else(|| {
            anyhow::anyhow!("No strand found for CDS in sequence {}", record.seqname())
        })?;
        let cds = Cds::new(start, end, phase);
        if !variant_coverage.get(ensp).copied().unwrap_or(false) {
            let has_variant = cds.contains_variant(&target, &variants);
            if has_variant {
                variant_coverage.insert(ensp.to_string(), true);
            }
        }
        let transcript = transcripts.entry(ensp.to_string()).or_insert_with(|| {
            Transcript::new(ensp.to_string(), target.clone(), strand, Vec::new())
        });
        transcript.coding_sequences.push(cds);
    }
    Ok(transcripts
        .into_values()
        .filter(|t| matches!(variant_coverage.get(&t.feature), Some(true)))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::duck::write_graphs;
    use crate::graph::node::Node;
    use crate::graph::Edge;
    use crate::translation::amino_acids::AminoAcid;
    use bio::stats::Prob;
    use petgraph::{Directed, Graph};

    #[test]
    fn name_formats_transcript_correctly() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "chr1".to_string(),
            Strand::Forward,
            vec![Cds::new(1, 10, 0)],
        );
        assert_eq!(transcript.name(), "chr1:ENSP00000493376");
    }

    fn setup_graph() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let alt_node_vaf_1 = HashMap::from([("s1".to_string(), 0.5)]);
        let alt_node_vaf_2 = HashMap::from([("s1".to_string(), 0.3)]);
        let event_probs = EventProbs(HashMap::from([
            ("PROB_GERMLINE".to_string(), LogProb::from(Prob(0.3))),
            ("PROB_SOMATIC".to_string(), LogProb::from(Prob(0.8))),
        ]));
        let alt_node_1 = graph.add_node(Node {
            node_type: NodeType::Variant,
            reference_allele: "A".to_string(),
            alternative_allele: "".to_string(),
            vaf: alt_node_vaf_1,
            probs: event_probs.clone(),
            pos: 1,
            index: 0,
        });
        let alt_node_2 = graph.add_node(Node {
            node_type: NodeType::Variant,
            reference_allele: "G".to_string(),
            alternative_allele: "G".to_string(),
            vaf: alt_node_vaf_2.clone(),
            probs: event_probs.clone(),
            pos: 4,
            index: 1,
        });
        let alt_node_3 = graph.add_node(Node {
            node_type: NodeType::Variant,
            reference_allele: "C".to_string(),
            alternative_allele: "A".to_string(),
            vaf: alt_node_vaf_2.clone(),
            probs: event_probs.clone(),
            pos: 14,
            index: 3,
        });
        let ref_node_1 = graph.add_node(Node {
            node_type: NodeType::Reference,
            reference_allele: "".to_string(),
            alternative_allele: "".to_string(),
            vaf: alt_node_vaf_2,
            probs: event_probs.clone(),
            pos: 12,
            index: 2,
        });
        graph.add_edge(
            alt_node_1,
            alt_node_2,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        graph.add_edge(
            ref_node_1,
            alt_node_3,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        VariantGraph {
            graph,
            start: 0,
            end: 15,
            target: "chr1".to_string(),
        }
    }

    #[test]
    fn rna_generates_correct_paths_for_forward_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(0, 10, 0), Cds::new(12, 15, 0)],
        );
        let graph = setup_graph();
        let tmp = tempfile::tempdir().unwrap();
        let graph_path = tmp.path().join("graph.duckdb");
        write_graphs(HashMap::from([("test".to_string(), graph)]), &graph_path).unwrap();
        let reference = HashMap::from([(
            "test".to_string(),
            vec![
                b'A', b'T', b'G', b'C', b'A', b'T', b'G', b'C', b'A', b'T', b'T', b'T', b'T', b'T',
                b'T', b'T', b'T', b'T',
            ],
        )]);
        let rna_paths = transcript.rna(&graph_path, &reference, "s1").unwrap();
        assert_eq!(rna_paths.len(), 1);
        assert_eq!(
            rna_paths[0].sequence,
            vec![
                b'A', b'G', b'C', b'G', b'T', b'G', b'C', b'A', b'T', b'T', b'T', b'T', b'A', b'T'
            ]
        );
    }

    #[test]
    fn peptides_generate_correctly_for_valid_input() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(0, 10, 0)],
        );
        let graph = setup_graph();
        let tmp = tempfile::tempdir().unwrap();
        let graph_path = tmp.path().join("graph.duckdb");
        write_graphs(HashMap::from([("test".to_string(), graph)]), &graph_path).unwrap();
        let reference = HashMap::from([(
            "test".to_string(),
            vec![
                b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T',
            ],
        )]);
        let interval = Interval { start: 2, end: 2 };
        let events = vec!["SOMATIC".to_string()];
        let background_events = vec!["GERMLINE".to_string()];
        let peptides = transcript
            .peptides(
                &graph_path,
                &reference,
                interval,
                "s1",
                &events,
                LogProb::from(Prob(0.4)),
                &background_events,
                LogProb::from(Prob(0.1)),
            )
            .unwrap();
        let peptide_sequences = peptides
            .iter()
            .map(|p| p.sequence.clone())
            .unique()
            .collect::<Vec<Vec<AminoAcid>>>();
        let expected = vec![
            vec![AminoAcid::Phenylalanine, AminoAcid::Valine],
            vec![AminoAcid::Valine, AminoAcid::Phenylalanine],
        ];
        assert_eq!(peptide_sequences, expected);
    }

    #[test]
    fn transcripts_parses_gff_file_correctly() {
        let gff_content = "\
##gff-version 3
chr1\tsource\tCDS\t1\t100\t.\t+\t0\tID=ENSP00000493376
chr1\tsource\tCDS\t200\t300\t.\t+\t0\tID=ENSP00000493376
chr1\tsource\tCDS\t400\t500\t.\t-\t0\tID=ENSP00000493377
        ";
        let tmp = tempfile::tempdir().unwrap();
        let gff_path = tmp.path().join("test.gff");
        std::fs::write(&gff_path, gff_content).unwrap();
        let graph_path = tmp.path().join("graph.duckdb");
        write_graphs(
            HashMap::from([("chr1".to_string(), setup_graph())]),
            &graph_path,
        )
        .unwrap();

        let transcripts = transcripts(&gff_path, &graph_path).unwrap();
        assert_eq!(transcripts.len(), 1);

        let transcript1 = transcripts
            .iter()
            .find(|t| t.feature == "ENSP00000493376")
            .unwrap();
        assert_eq!(transcript1.coding_sequences.len(), 2);
        assert_eq!(transcript1.strand, Strand::Forward);
    }

    #[test]
    fn start_returns_min_cds_start() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(10, 20, 0), Cds::new(5, 15, 0)],
        );
        assert_eq!(transcript.start().unwrap(), 5);
    }

    #[test]
    fn end_returns_max_cds_end() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(10, 20, 0), Cds::new(5, 25, 0)],
        );
        assert_eq!(transcript.end().unwrap(), 25);
    }

    #[test]
    fn cds_returns_sorted_cds_for_forward_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(10, 20, 0), Cds::new(5, 15, 0)],
        );
        let cds: Vec<&Cds> = transcript.cds().collect();
        assert_eq!(cds.len(), 2);
        assert_eq!(cds[0].start, 5);
        assert_eq!(cds[1].start, 10);
    }

    #[test]
    fn cds_returns_sorted_cds_for_reverse_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Reverse,
            vec![Cds::new(10, 20, 0), Cds::new(5, 15, 0)],
        );
        let cds: Vec<&Cds> = transcript.cds().collect();
        assert_eq!(cds.len(), 2);
        assert_eq!(cds[0].start, 10);
        assert_eq!(cds[1].start, 5);
    }

    #[test]
    fn cds_returns_empty_iterator_when_no_cds() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![],
        );
        let cds: Vec<&Cds> = transcript.cds().collect();
        assert!(cds.is_empty());
    }

    #[test]
    fn is_translatable_detects_non_acgt_in_coding_sequence() {
        let mut reference = HashMap::new();
        reference.insert("test".to_string(), b"ACGTNACGT".to_vec());
        let transcript = Transcript::new(
            "ENSP1".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(0, 8, 0)],
        );
        assert!(!transcript.is_translatable(&reference));
        assert_eq!(
            transcript.untranslatable_reference_base(&reference),
            Some(b'N')
        );
    }

    #[test]
    fn is_translatable_ignores_bases_outside_coding_sequence() {
        let mut reference = HashMap::new();
        reference.insert("test".to_string(), b"NACGTACGTN".to_vec());
        let transcript = Transcript::new(
            "ENSP1".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(1, 8, 0)],
        );
        assert!(transcript.is_translatable(&reference));
    }

    #[test]
    fn paths_returns_paths_for_forward_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Forward,
            vec![Cds::new(1, 10, 0)],
        );
        let graph = setup_graph();
        let paths = transcript.paths(&graph).unwrap();
        assert_eq!(paths, graph.paths());
    }

    #[test]
    fn paths_returns_reverse_paths_for_reverse_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "test".to_string(),
            Strand::Reverse,
            vec![Cds::new(1, 10, 0)],
        );
        let graph = setup_graph();
        let paths = transcript.paths(&graph).unwrap();
        assert_eq!(paths, graph.reverse_paths());
    }

    #[test]
    fn paths_returns_error_for_unknown_strand() {
        let transcript = Transcript::new(
            "ENSP00000493376".to_string(),
            "chr1".to_string(),
            Strand::Unknown,
            vec![Cds::new(1, 10, 0)],
        );
        let graph = setup_graph();
        let result = transcript.paths(&graph);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Strand is unknown for transcript chr1:ENSP00000493376"
        );
    }

    #[test]
    fn test_position_in_transcript() {
        let transcript = Transcript {
            feature: "Test".to_string(),
            target: "chr1".to_string(),
            strand: Strand::Forward,
            coding_sequences: vec![
                Cds {
                    start: 100,
                    end: 104,
                    phase: 0,
                },
                Cds {
                    start: 200,
                    end: 202,
                    phase: 0,
                },
            ],
        };

        assert_eq!(transcript.position_in_transcript(100).unwrap(), 0);
        assert_eq!(transcript.position_in_transcript(104).unwrap(), 4);
        assert_eq!(transcript.position_in_transcript(200).unwrap(), 5);
        assert_eq!(transcript.position_in_transcript(202).unwrap(), 7);
        assert!(transcript.position_in_transcript(150).is_err());
    }

    #[test]
    fn test_transcript_reader() {
        let gff_file = PathBuf::from("tests/resources/ENSP00000355304.gff3");
        let graph = PathBuf::from("tests/resources/graph.duckdb");
        let transcripts = transcripts(&gff_file, &graph).unwrap();
        assert_eq!(transcripts.len(), 1)
    }
}
