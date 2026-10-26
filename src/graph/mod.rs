pub(crate) mod combinations;
pub(crate) mod duck;
pub(crate) mod hgvs;
pub(crate) mod node;
pub(crate) mod paths;
pub(crate) mod peptide;
pub(crate) mod score;
pub(crate) mod transcript;

use crate::cli::ObservationFile;
use crate::graph::node::{Node, NodeType};
use crate::graph::paths::HaplotypePath;
use crate::utils::bcf::extract_event_names;
use crate::utils::NUMERICAL_EPSILON;
use anyhow::Result;
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::bayesian::BayesFactor;
use bio::stats::{LogProb, PHREDProb};
use itertools::Itertools;
use log::{info, warn};
use petgraph::graph::NodeIndex;
use petgraph::{Directed, Graph};
use rust_htslib::bcf::{Read, Reader, Record};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::ops::RangeInclusive;
use std::path::PathBuf;
use varlociraptor::calling::variants::preprocessing::read_observations;
use varlociraptor::utils::collect_variants::collect_variants;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct VariantGraph {
    pub(crate) graph: Graph<Node, Edge, Directed>,
    start: i64,
    end: i64,
    target: String,
}

impl VariantGraph {
    pub(crate) fn build(
        calls_file: &PathBuf,
        observation_files: &[ObservationFile],
        target: &str,
        min_prob_present: LogProb,
        min_vaf: f32,
    ) -> Result<VariantGraph> {
        let mut calls_reader = Reader::from_path(calls_file)?;
        let header = calls_reader.header().clone();
        let mut nodes_by_index = HashMap::new();
        let mut observation_readers: HashMap<_, _> = observation_files
            .iter()
            .map(|o| (o.sample.to_string(), Reader::from_path(&o.path).unwrap()))
            .collect();
        let mut observations_records = observation_readers
            .iter_mut()
            .map(|(sample, reader)| (sample.clone(), reader.records()))
            .collect::<HashMap<_, _>>();

        let mut samples = calls_reader
            .header()
            .samples()
            .iter()
            .map(|s| String::from_utf8(s.to_vec()).unwrap())
            .collect_vec();

        let observation_samples = observation_files
            .iter()
            .map(|o| o.sample.clone())
            .collect_vec();

        for sample in &observation_samples {
            if !samples.contains(sample) {
                warn!("Sample {sample} in observations file is not present in calls file");
                samples.retain(|s| s != sample);
            }
        }

        let event_names = extract_event_names(calls_file);

        let tags = event_names
            .iter()
            .map(|event| format!("PROB_{event}"))
            .collect();

        let mut supporting_reads = HashMap::new();

        let mut variant_graph = Graph::<Node, Edge, Directed>::new();
        let mut index = 0;
        let mut last_position = -1;
        let mut start = 0;

        info!("Adding nodes for target {target}.");

        for calls_record in calls_reader.records() {
            let mut calls_record = calls_record?;
            let position = calls_record.pos();
            if last_position == -1 {
                start = position; // Set start position
            }
            if last_position != position {
                index += 1;
            }
            let mut observations_records = observations_records
                .iter_mut()
                .map(|(sample, records)| {
                    let record = records
                        .next()
                        .unwrap_or_else(|| {
                            panic!(
                                "Missing observation record for calls record at position {position}"
                            )
                        })
                        .unwrap();
                    (sample, record)
                })
                .collect::<HashMap<_, _>>();

            if header.rid2name(calls_record.rid().unwrap()).unwrap() != target.as_bytes() {
                continue;
            }

            let _variants = collect_variants(&mut calls_record, false, None, None, None)?;
            let observations = observations_records
                .iter_mut()
                .map(|(sample, record)| {
                    let observations = read_observations(record).unwrap();
                    (sample, observations.pileup.read_observations().clone())
                })
                .collect::<HashMap<_, _>>();

            let _fragment_ids: HashSet<_> = observations
                .iter()
                .flat_map(|(s, v)| v.iter().map(move |o| (s, o.fragment_id)))
                .collect();

            let alleles = calls_record.alleles();
            let alt_allele = String::from_utf8(alleles[1].to_vec())?;
            if alt_allele == "*" {
                continue;
            }
            if !alt_allele
                .bytes()
                .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            {
                warn!(
                    "Skipping variant at {target}:{} with non-ACGT alternative allele {alt_allele}",
                    position + 1
                );
                continue;
            }

            let event_probs = EventProbs::from_record(&calls_record, &tags);
            if !event_probs.is_valid() {
                return Err(anyhow::anyhow!(
                    "Invalid event probabilities in record at position {}",
                    position
                ));
            } else if event_probs.all_nan() || event_probs.prob_present()? < min_prob_present {
                continue;
            }

            let var_node = Node::from_records(
                &calls_record,
                &observations,
                &event_probs,
                NodeType::Variant,
                &samples,
                index,
            );
            if var_node.max_vaf() < min_vaf {
                continue;
            }
            let var_node_index = variant_graph.add_node(var_node);
            nodes_by_index
                .entry(index)
                .or_insert(Vec::new())
                .push(var_node_index);

            let mut ref_node_index = None;
            if last_position != position {
                let ref_node = Node::from_records(
                    &calls_record,
                    &observations,
                    &event_probs,
                    NodeType::Reference,
                    &samples,
                    index,
                );
                ref_node_index = Some(variant_graph.add_node(ref_node));
                nodes_by_index
                    .entry(index)
                    .or_insert(Vec::new())
                    .push(ref_node_index.unwrap());
            }

            for (sample, observations) in observations {
                for observation in observations {
                    let supports_allele = |numerator, denominator| {
                        matches!(
                            BayesFactor::new(numerator, denominator).evidence_kass_raftery(),
                            KassRaftery::Positive | KassRaftery::Strong | KassRaftery::VeryStrong
                        )
                    };
                    if supports_allele(observation.prob_alt(), observation.prob_ref()) {
                        // Read supports variant
                        supporting_reads
                            .entry((sample.to_string(), observation.fragment_id))
                            .or_insert(Vec::new())
                            .push(var_node_index);
                    } else if supports_allele(observation.prob_ref(), observation.prob_alt()) {
                        // Read supports reference
                        if let Some(index) = ref_node_index {
                            supporting_reads
                                .entry((sample.to_string(), observation.fragment_id))
                                .or_insert(Vec::new())
                                .push(index);
                        }
                    }
                }
            }
            last_position = position;
        }

        let mut variant_graph = VariantGraph {
            graph: variant_graph,
            start,
            end: last_position,
            target: target.to_string(),
        };

        // For each (sample, fragment_id), if this fragment supports
        // any variant node at position P, remove its ref node entries at position P.
        // This fixes false ref support from normalized multi-allelic sites where a
        // fragment supporting alt allele B is incorrectly recorded as ref support
        // when processing alt allele A at the same position.
        for nodes in supporting_reads.values_mut() {
            let variant_positions: HashSet<i64> = nodes
                .iter()
                .filter(|&&n| variant_graph.graph[n].node_type == NodeType::Variant)
                .map(|&n| variant_graph.graph[n].pos)
                .collect();
            nodes.retain(|&n| {
                variant_graph.graph[n].node_type != NodeType::Reference
                    || !variant_positions.contains(&variant_graph.graph[n].pos)
            });
        }

        let mut possible_node_pairs: HashSet<(NodeIndex, NodeIndex)> = HashSet::new();
        let mut sorted_indices: Vec<_> = nodes_by_index.keys().collect();
        sorted_indices.sort_unstable();

        for (idx_a, idx_b) in sorted_indices.iter().tuple_windows() {
            for &node_a in &nodes_by_index[idx_a] {
                for &node_b in &nodes_by_index[idx_b] {
                    possible_node_pairs.insert((node_a, node_b));
                }
            }
        }

        info!("Finished adding nodes for target {target}.");
        variant_graph.connect_consecutive_positions(&possible_node_pairs);

        info!("Adding read support for target {target}.");
        variant_graph.add_read_support(&supporting_reads, &possible_node_pairs)?;

        info!("Removing edges without evidence for target {target}.");
        variant_graph.prune_edges_without_evidence(&nodes_by_index)?;

        Ok(variant_graph)
    }

    pub(crate) fn prune_edges_without_evidence(
        &mut self,
        nodes_by_index: &HashMap<u32, Vec<NodeIndex>>,
    ) -> Result<()> {
        let sorted_positions: Vec<_> = nodes_by_index.keys().sorted().collect();
        for (&idx_a, &idx_b) in sorted_positions.iter().tuple_windows() {
            // Remove all edges without read support if and only if there is at least one edge with read support between any nodes at idx_a and idx_b
            if self.has_read_support(&nodes_by_index[idx_a], &nodes_by_index[idx_b]) {
                for &node_a in &nodes_by_index[idx_a] {
                    for &node_b in &nodes_by_index[idx_b] {
                        if let Some(edge) = self.graph.find_edge(node_a, node_b) {
                            let edge_weight = self.graph.edge_weight(edge).unwrap();
                            if edge_weight.supporting_reads.values().sum::<u32>() == 0 {
                                self.graph.remove_edge(edge);
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    pub(crate) fn connect_consecutive_positions(
        &mut self,
        possible_node_pairs: &HashSet<(NodeIndex, NodeIndex)>,
    ) {
        for &(l, r) in possible_node_pairs {
            if self.graph.find_edge(l, r).is_none() {
                self.graph.add_edge(
                    l,
                    r,
                    Edge {
                        supporting_reads: HashMap::new(),
                    },
                );
            }
        }
    }

    pub(crate) fn add_read_support(
        &mut self,
        supporting_reads: &HashMap<(String, Option<u64>), Vec<NodeIndex>>,
        possible_node_pairs: &HashSet<(NodeIndex, NodeIndex)>,
    ) -> Result<()> {
        for ((sample, _), nodes) in supporting_reads {
            for [a, b] in nodes.iter().copied().array_combinations::<2>() {
                let (a, b) = if possible_node_pairs.contains(&(a, b)) {
                    (a, b)
                } else if possible_node_pairs.contains(&(b, a)) {
                    (b, a)
                } else {
                    continue;
                };
                let edge = self.graph.find_edge(a, b);
                if let Some(edge) = edge {
                    self.graph
                        .edge_weight_mut(edge)
                        .unwrap()
                        .supporting_reads
                        .entry(sample.to_string())
                        .and_modify(|v| *v += 1)
                        .or_insert(1);
                } else {
                    self.graph.add_edge(
                        a,
                        b,
                        Edge {
                            supporting_reads: HashMap::from([(sample.to_string(), 1)]),
                        },
                    );
                }
            }
        }
        Ok(())
    }

    /// Finds all paths starting from all nodes with the lowest graph index.
    ///
    /// This method performs a depth-first search (DFS) starting from from every node whose graph index is the minimum index in the graph.
    /// It collects all possible paths, including those that reach leaf nodes
    /// (nodes with no outgoing edges). The paths are represented as vectors of `NodeIndex`.
    ///
    /// Returns:
    ///     A vector of vectors, where each inner vector represents a path through the graph
    ///     starting from one of the lowest-index start nodes. Each path is a sequence of `NodeIndex`
    ///     elements, representing the nodes in the order they are visited
    pub(crate) fn paths(&self) -> Vec<HaplotypePath> {
        let mut all_paths = Vec::new();
        let min_index = match self.graph.node_indices().map(|i| self.graph[i].index).min() {
            Some(min_index) => min_index,
            None => return Vec::new(),
        };
        let start_nodes = self
            .graph
            .node_indices()
            .filter(|&i| self.graph[i].index == min_index)
            .collect::<Vec<_>>();
        for start_node in start_nodes {
            let mut stack = vec![(start_node, vec![start_node])];

            while let Some((node, mut path)) = stack.pop() {
                // Use the forward-only filter to only traverse neighbors with higher indices.
                let mut found_forward = false;
                for neighbor in self
                    .graph
                    .neighbors(node)
                    .filter(|&n| self.graph[n].index > self.graph[node].index)
                {
                    if path.contains(&neighbor) {
                        continue;
                    }
                    found_forward = true;
                    path.push(neighbor);
                    stack.push((neighbor, path.clone()));
                    path.pop();
                }
                // If no valid forward neighbor is found, we've reached a "leaf" in terms of forward traversal.
                if !found_forward {
                    all_paths.push(HaplotypePath(path));
                }
            }
        }

        all_paths.retain(|path| {
            let nodes = path
                .0
                .iter()
                .map(|n| self.graph.node_weight(*n).unwrap())
                .collect_vec();
            nodes.iter().all(|n| n.pos != -1)
        });

        all_paths.into_iter().unique().collect_vec()
    }

    pub(crate) fn reverse_paths(&self) -> Vec<HaplotypePath> {
        self.paths()
            .iter()
            .map(|path| HaplotypePath(path.0.iter().rev().cloned().collect()))
            .collect()
    }

    /// Whether a read links any of the `from` nodes to any of the `to` nodes.
    fn has_read_support(&self, from: &[NodeIndex], to: &[NodeIndex]) -> bool {
        from.iter().cartesian_product(to).any(|(&a, &b)| {
            self.graph
                .find_edge(a, b)
                .is_some_and(|edge| self.graph[edge].supporting_reads.values().sum::<u32>() > 0)
        })
    }

    /// Runs of positions linked by reads, as ranges of position indices.
    pub(crate) fn phase_blocks(&self) -> Vec<RangeInclusive<u32>> {
        let mut nodes_by_index: BTreeMap<u32, Vec<NodeIndex>> = BTreeMap::new();
        for node in self.graph.node_indices() {
            nodes_by_index
                .entry(self.graph[node].index)
                .or_default()
                .push(node);
        }
        let mut blocks: Vec<RangeInclusive<u32>> = Vec::new();
        for (&index, nodes) in &nodes_by_index {
            match blocks.last_mut() {
                Some(block) if self.has_read_support(&nodes_by_index[block.end()], nodes) => {
                    *block = *block.start()..=index
                }
                _ => blocks.push(index..=index),
            }
        }
        blocks
    }

    /// The `k` best supported paths through the phase block.
    pub(crate) fn top_k_paths(&self, block: &RangeInclusive<u32>, k: usize) -> Vec<HaplotypePath> {
        let (first, last) = (*block.start(), *block.end());
        let start_nodes: Vec<NodeIndex> = self
            .graph
            .node_indices()
            .filter(|&i| self.graph[i].index == first)
            .collect();

        // (path, no_zero_edges, min_nonzero_reads)
        type ScoredPath = (Vec<NodeIndex>, bool, u32);

        let score = |path: &[NodeIndex]| -> (bool, u32) {
            let mut min_nz = u32::MAX;
            let mut has_zero = false;
            for w in path.windows(2) {
                if let Some(e) = self
                    .graph
                    .find_edge(w[0], w[1])
                    .or_else(|| self.graph.find_edge(w[1], w[0]))
                {
                    let total: u32 = self
                        .graph
                        .edge_weight(e)
                        .unwrap()
                        .supporting_reads
                        .values()
                        .sum();
                    if total == 0 {
                        has_zero = true;
                    } else {
                        min_nz = min_nz.min(total);
                    }
                }
            }
            let min_nz = if min_nz == u32::MAX { 0 } else { min_nz };
            (!has_zero, min_nz)
        };

        let mut beam: Vec<ScoredPath> = start_nodes
            .iter()
            .map(|&n| (vec![n], true, u32::MAX))
            .collect();

        let mut complete: Vec<ScoredPath> = Vec::new();

        while !beam.is_empty() {
            let mut next_beam: Vec<ScoredPath> = Vec::new();

            for (path, _, _) in beam {
                let current = *path.last().unwrap();
                let neighbors: Vec<NodeIndex> = self
                    .graph
                    .neighbors(current)
                    .filter(|&n| self.graph[n].index > self.graph[current].index)
                    .filter(|&n| self.graph[n].index <= last)
                    .filter(|&n| !path.contains(&n))
                    .collect();

                if neighbors.is_empty() {
                    if self.graph[current].pos != -1 {
                        let (no_zero, min_nz) = score(&path);
                        complete.push((path, no_zero, min_nz));
                    }
                } else {
                    for neighbor in neighbors {
                        let mut new_path = path.clone();
                        new_path.push(neighbor);
                        let (no_zero, min_nz) = score(&new_path);
                        next_beam.push((new_path, no_zero, min_nz));
                    }
                }
            }

            // Keep only top-k in beam
            next_beam.sort_by(|a, b| b.1.cmp(&a.1).then(b.2.cmp(&a.2)));
            next_beam.truncate(k);
            beam = next_beam;
        }

        complete.sort_by(|a, b| b.1.cmp(&a.1).then(b.2.cmp(&a.2)));
        complete.retain(|(path, _, _)| {
            let last_node = *path.last().unwrap();
            self.graph[last_node].index == last
        });
        complete.truncate(k);
        complete
            .into_iter()
            .map(|(path, _, _)| HaplotypePath(path))
            .unique()
            .collect()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.graph.node_count() == 0
    }

    pub(crate) fn edge_reads(&self, path: &[NodeIndex]) -> Vec<HashMap<String, u32>> {
        path.windows(2)
            .map(|w| {
                let edge = self
                    .graph
                    .find_edge(w[0], w[1])
                    .or_else(|| self.graph.find_edge(w[1], w[0]))
                    .expect("Consecutive nodes in a valid path must share an edge.");
                self.graph
                    .edge_weight(edge)
                    .expect("Edge exists in graph but has no weight — graph is malformed.")
                    .supporting_reads
                    .clone()
            })
            .collect()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct EventProbs(pub(crate) HashMap<String, LogProb>);

impl EventProbs {
    fn from_record(record: &Record, tags: &Vec<String>) -> Self {
        let mut probs = HashMap::new();
        for tag in tags {
            let prob = record.info(tag.as_bytes()).float().unwrap().unwrap()[0];
            probs.insert(tag.to_string(), LogProb::from(PHREDProb(prob as f64)));
        }
        EventProbs(probs)
    }

    pub(crate) fn all_nan(&self) -> bool {
        self.0.values().all(|v| v.is_nan())
    }

    pub(crate) fn is_valid(&self) -> bool {
        self.all_nan() || self.0.values().all(|v| !v.is_nan())
    }

    /// Returns the probability of the variant being present by calculating 1 - (PROB_ABSENT + PROB_ARTIFACT)
    pub(crate) fn prob_present(&self) -> Result<LogProb> {
        Ok(self
            .0
            .get("PROB_ABSENT")
            .unwrap()
            .ln_add_exp(*self.0.get("PROB_ARTIFACT").unwrap())
            .cap_numerical_overshoot(NUMERICAL_EPSILON)
            .ln_one_minus_exp())
    }

    pub(crate) fn prob(&self, event: &str) -> Result<LogProb> {
        Ok(*self
            .0
            .get(format!("PROB_{event}").as_str())
            .ok_or_else(|| anyhow::anyhow!("Event '{}' not found", event))?)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct Edge {
    pub(crate) supporting_reads: HashMap<String, u32>,
}

// test node distance
#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::node::{Node, NodeType};
    use bio::stats::Prob;
    use petgraph::{Directed, Graph};
    use rust_htslib::bcf::{Read, Reader};
    use std::fs;

    #[test]
    fn test_event_probs_from_record() {
        let mut reader = Reader::from_path("tests/resources/calls.bcf").unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let tags = vec![
            "PROB_ABSENT".to_string(),
            "PROB_PRESENT".to_string(),
            "PROB_ARTIFACT".to_string(),
        ];
        let event_probs = EventProbs::from_record(&record, &tags);
        assert_eq!(event_probs.0.len(), 3);
    }

    #[test]
    fn all_nan_returns_true_when_all_values_are_nan() {
        let event_probs = EventProbs(HashMap::from([
            ("PROB_1".to_string(), LogProb::from(Prob(f64::NAN))),
            ("PROB_2".to_string(), LogProb::from(Prob(f64::NAN))),
        ]));
        assert!(event_probs.all_nan());
    }

    #[test]
    fn all_nan_returns_false_when_not_all_values_are_nan() {
        let event_probs = EventProbs(HashMap::from([
            ("PROB_1".to_string(), LogProb::from(Prob(f64::NAN))),
            ("PROB_2".to_string(), LogProb::from(Prob(0.5))),
        ]));
        assert!(!event_probs.all_nan());
    }

    #[test]
    fn is_valid_returns_true_when_all_values_are_nan() {
        let event_probs = EventProbs(HashMap::from([
            ("PROB_1".to_string(), LogProb::from(Prob(f64::NAN))),
            ("PROB_2".to_string(), LogProb::from(Prob(f64::NAN))),
        ]));
        assert!(event_probs.is_valid());
    }

    #[test]
    fn is_valid_returns_true_when_all_values_are_finite() {
        let event_probs = EventProbs(HashMap::from([
            ("PROB_1".to_string(), LogProb::from(Prob(0.5))),
            ("PROB_2".to_string(), LogProb::from(Prob(1.0))),
        ]));
        assert!(event_probs.is_valid());
    }

    #[test]
    fn is_valid_returns_false_when_some_values_are_nan() {
        let event_probs = EventProbs(HashMap::from([
            ("PROB_1".to_string(), LogProb::from(Prob(0.5))),
            ("PROB_2".to_string(), LogProb::from(Prob(f64::NAN))),
        ]));
        assert!(!event_probs.is_valid());
    }

    #[test]
    fn prob_present_calculates_correctly_with_valid_values() {
        let event_probs = EventProbs(HashMap::from([
            ("PROB_ABSENT".to_string(), LogProb::from(Prob(0.1))),
            ("PROB_ARTIFACT".to_string(), LogProb::from(Prob(0.2))),
        ]));
        assert!(
            event_probs.prob_present().unwrap() < LogProb::from(Prob(0.71))
                && event_probs.prob_present().unwrap() > LogProb::from(Prob(0.69))
        );
    }

    #[test]
    fn test_build_graph() {
        let calls_file = PathBuf::from("tests/resources/test_calls.vcf");
        let observations_file = PathBuf::from("tests/resources/test_observations.vcf");
        let observations = vec![ObservationFile {
            path: observations_file,
            sample: "sample".to_string(),
        }];
        let variant_graph = VariantGraph::build(
            &calls_file,
            &observations,
            "OX512233.1",
            LogProb::from(Prob(0.0)),
            0.05,
        );
        assert!(variant_graph.is_ok());
    }

    #[test]
    fn test_graph_is_empty() {
        let calls_file = PathBuf::from("tests/resources/test_calls.vcf");
        let observations_file = PathBuf::from("tests/resources/test_observations.vcf");
        let observations = vec![ObservationFile {
            path: observations_file,
            sample: "sample".to_string(),
        }];
        let variant_graph = VariantGraph::build(
            &calls_file,
            &observations,
            "not actually in file",
            LogProb::from(Prob(0.0)),
            0.05,
        );
        assert!(variant_graph.unwrap().is_empty());
    }

    #[test]
    fn test_graph_paths() {
        let calls_file = PathBuf::from("tests/resources/test_calls.vcf");
        let observations_file = PathBuf::from("tests/resources/test_observations.vcf");
        let observations = vec![ObservationFile {
            path: observations_file,
            sample: "sample".to_string(),
        }];
        let mut variant_graph = VariantGraph::build(
            &calls_file,
            &observations,
            "OX512233.1",
            LogProb::from(Prob(0.0)),
            0.00,
        )
        .unwrap();
        variant_graph.graph.add_edge(
            NodeIndex::new(0),
            NodeIndex::new(2),
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        variant_graph.graph.add_edge(
            NodeIndex::new(2),
            NodeIndex::new(5),
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        variant_graph.graph.add_edge(
            NodeIndex::new(5),
            NodeIndex::new(6),
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        let paths = variant_graph.paths();
        assert_eq!(paths.len(), 4);
    }

    pub(crate) fn setup_variant_graph_with_nodes() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let _node1 = graph.add_node(Node::new(
            NodeType::Variant,
            1,
            "G".to_string(),
            "A".to_string(),
        ));
        let _node2 = graph.add_node(Node::new(
            NodeType::Reference,
            2,
            "".to_string(),
            "".to_string(),
        ));
        let _node3 = graph.add_node(Node::new(
            NodeType::Variant,
            3,
            "C".to_string(),
            "T".to_string(),
        ));
        let _node4 = graph.add_node(Node::new(
            NodeType::Variant,
            4,
            "C".to_string(),
            "".to_string(),
        ));
        let _node5 = graph.add_node(Node::new(
            NodeType::Variant,
            8,
            "C".to_string(),
            "A".to_string(),
        ));
        let _node6 = graph.add_node(Node::new(
            NodeType::Variant,
            9,
            "A".to_string(),
            "TT".to_string(),
        ));
        VariantGraph {
            graph,
            start: 0,
            end: 2,
            target: "test".to_string(),
        }
    }

    #[test]
    fn is_variant_returns_true_for_variant_node() {
        let node_type = NodeType::Variant;
        assert!(node_type.is_variant());
    }

    #[test]
    fn is_variant_returns_false_for_reference_node() {
        let node_type = NodeType::Reference;
        assert!(!node_type.is_variant());
    }

    #[test]
    fn test_graph_serialization_and_deserialization() {
        let graph = setup_variant_graph_with_nodes();
        let temp_dir = tempfile::tempdir().unwrap();
        let file_path = temp_dir.path().join("graph.json");
        let file = fs::File::create(&file_path).unwrap();
        serde_json::to_writer(file, &graph).unwrap();
        let file = fs::File::open(&file_path).unwrap();
        let deserialized_graph: VariantGraph = serde_json::from_reader(file).unwrap();
        assert_eq!(
            graph.graph.node_count(),
            deserialized_graph.graph.node_count()
        );
    }

    fn setup_graph_with_edges() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let _node0 = graph.add_node(Node::new(
            NodeType::Reference,
            1,
            "C".to_string(),
            "A".to_string(),
        ));
        let node1 = graph.add_node(Node::new(
            NodeType::Variant,
            1,
            "C".to_string(),
            "A".to_string(),
        ));
        let node2 = graph.add_node(Node::new(
            NodeType::Variant,
            2,
            "T".to_string(),
            "T".to_string(),
        ));
        let node3 = graph.add_node(Node::new(
            NodeType::Variant,
            3,
            "T".to_string(),
            "G".to_string(),
        ));
        graph.add_edge(
            node1,
            node2,
            Edge {
                supporting_reads: HashMap::from([("s1".to_string(), 3u32)]),
            },
        );
        graph.add_edge(
            node2,
            node3,
            Edge {
                supporting_reads: HashMap::from([("s1".to_string(), 5u32)]),
            },
        );
        VariantGraph {
            graph,
            start: 0,
            end: 3,
            target: "test".to_string(),
        }
    }

    #[test]
    fn edge_reads_returns_reads_for_forward_path() {
        let g = setup_graph_with_edges();
        let idx = g.graph.node_indices().skip(1).collect_vec();
        assert_eq!(g.edge_reads(&idx).len(), 2);
    }

    #[test]
    fn edge_reads_returns_reads_for_reverse_path() {
        let g = setup_graph_with_edges();
        let idx = g.graph.node_indices().skip(1).rev().collect_vec();
        let reads = g.edge_reads(&idx);
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].get("s1"), Some(&5));
    }

    /// A graph with a reference and a variant node per variant. `edges` links consecutive
    /// variants as `(index, from_variant, to_variant, reads)`.
    pub(crate) fn phased_graph(
        variants: &[(i64, &str, &str)],
        edges: &[(usize, bool, bool, u32)],
    ) -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let nodes: Vec<[NodeIndex; 2]> = variants
            .iter()
            .enumerate()
            .map(|(index, (pos, reference_allele, alternative_allele))| {
                [NodeType::Reference, NodeType::Variant].map(|node_type| {
                    graph.add_node(Node {
                        node_type,
                        reference_allele: reference_allele.to_string(),
                        alternative_allele: alternative_allele.to_string(),
                        vaf: HashMap::from([("s1".to_string(), 0.5)]),
                        probs: EventProbs(HashMap::new()),
                        pos: *pos,
                        index: index as u32,
                    })
                })
            })
            .collect();
        for &(index, from_variant, to_variant, reads) in edges {
            graph.add_edge(
                nodes[index][from_variant as usize],
                nodes[index + 1][to_variant as usize],
                Edge {
                    supporting_reads: HashMap::from([("s1".to_string(), reads)]),
                },
            );
        }
        VariantGraph {
            graph,
            start: 0,
            end: 0,
            target: "test".to_string(),
        }
    }

    #[test]
    fn phase_blocks_split_without_reads() {
        let snvs = (0..4).map(|pos| (pos, "A", "C")).collect_vec();
        let mut edges = vec![(0, true, true, 5), (0, false, false, 3)];
        edges.extend([(2, true, false, 2), (2, false, true, 4)]);
        let graph = phased_graph(&snvs, &edges);

        let blocks = graph.phase_blocks();

        assert_eq!(blocks, vec![0..=1, 2..=3]);
        assert_eq!(graph.top_k_paths(&blocks[0], 10).len(), 2);
        assert_eq!(graph.top_k_paths(&blocks[0], 1).len(), 1);
    }
}
