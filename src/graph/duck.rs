use crate::cli::HgvsNotation;
use crate::graph::node::{Node, NodeType};
use crate::graph::score::{AnnotationInput, HaplotypeScore, ScoreRecord};
use crate::graph::transcript::Transcript;
use crate::graph::{Edge, VariantGraph};
use anyhow::Result;
use duckdb::{params, Connection};
use petgraph::Graph;
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::str::FromStr;

type NodeRow = (usize, String, String, String, String, String, i64, u32);

pub(crate) fn write_graphs(graphs: HashMap<String, VariantGraph>, path: &Path) -> Result<()> {
    let mut db = Connection::open(path)?;
    let transaction = db.transaction()?;
    transaction.execute("CREATE TABLE graphs (target STRING PRIMARY KEY, start_position INTEGER, end_position INTEGER)", [])?;
    transaction.execute(
        "CREATE TABLE nodes (target STRING, node_index INTEGER, node_type STRING, reference_allele STRING, alternative_allele STRING, vaf STRING, probs STRING, pos INTEGER, index INTEGER)",
        [],
    )?;
    transaction.execute(
        "CREATE TABLE edges (target STRING, from_node INTEGER, to_node INTEGER, supporting_reads STRING)",
        [],
    )?;

    let mut insert_graph = transaction.prepare("INSERT INTO graphs VALUES (?, ?, ?)")?;
    let mut insert_node =
        transaction.prepare("INSERT INTO nodes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)")?;
    let mut insert_edge = transaction.prepare("INSERT INTO edges VALUES (?, ?, ?, ?)")?;

    for (target, graph) in graphs {
        insert_graph.execute(params![target, graph.start, graph.end])?;
        for node_index in graph.graph.node_indices() {
            let node = graph.graph.node_weight(node_index).unwrap();
            insert_node.execute(params![
                target,
                node_index.index() as i64,
                node.node_type.to_string(),
                node.reference_allele,
                node.alternative_allele,
                json5::to_string(&node.vaf)?,
                json5::to_string(&node.probs)?,
                node.pos,
                node.index as i64,
            ])?;
        }
        for edge in graph.graph.edge_indices() {
            let (from_node, to_node) = graph.graph.edge_endpoints(edge).unwrap();
            let edge = graph.graph.edge_weight(edge).unwrap();
            insert_edge.execute(params![
                target,
                from_node.index() as i64,
                to_node.index() as i64,
                serde_json::to_string(&edge.supporting_reads)?,
            ])?;
        }
    }
    transaction.execute(
        "CREATE INDEX idx_nodes_target_pos ON nodes (target, pos)",
        [],
    )?;
    transaction.execute("CREATE INDEX idx_edges_target ON edges (target)", [])?;
    transaction.execute("CREATE INDEX idx_graphs_target ON graphs (target)", [])?;
    transaction.commit()?;
    db.close().unwrap();
    Ok(())
}

pub(crate) fn feature_graph(
    path: PathBuf,
    target: String,
    start: u64,
    end: u64,
) -> Result<VariantGraph> {
    let db = Connection::open(path)?;
    let mut graph = Graph::<Node, Edge, petgraph::Directed>::new();
    let mut stmt = db.prepare(
    "SELECT node_index, node_type, reference_allele, alternative_allele, vaf, probs, pos, index FROM nodes WHERE target = ? AND pos >= ? AND pos <= ? ORDER BY node_index"
    )?;
    let nodes: Vec<NodeRow> = stmt
        .query_map(
            params![target.to_string(), start as i64, end as i64],
            |row| {
                Ok((
                    row.get(0)?,
                    row.get(1)?,
                    row.get(2)?,
                    row.get(3)?,
                    row.get(4)?,
                    row.get(5)?,
                    row.get(6)?,
                    row.get(7)?,
                ))
            },
        )?
        .map(Result::unwrap)
        .collect();
    // Edges refer to nodes by the index they had in the full graph.
    let mut node_indices = HashMap::new();
    for (node_index, node_type, reference_allele, alternative_allele, vaf, probs, pos, index) in
        nodes
    {
        let node = Node {
            node_type: NodeType::from_str(&node_type)?,
            reference_allele,
            alternative_allele,
            vaf: json5::from_str(&vaf)?,
            probs: json5::from_str(&probs)?,
            pos,
            index,
        };
        node_indices.insert(node_index, graph.add_node(node));
    }
    let mut stmt = db.prepare(
        "SELECT edges.from_node, edges.to_node, edges.supporting_reads FROM edges \
         JOIN nodes source ON source.target = edges.target AND source.node_index = edges.from_node \
         JOIN nodes sink ON sink.target = edges.target AND sink.node_index = edges.to_node \
         WHERE edges.target = ? AND source.pos BETWEEN ? AND ? AND sink.pos BETWEEN ? AND ?",
    )?;
    let edges: Vec<(usize, usize, String)> = stmt
        .query_map(
            params![
                target.to_string(),
                start as i64,
                end as i64,
                start as i64,
                end as i64
            ],
            |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)),
        )?
        .collect::<Result<Vec<_>, _>>()?;
    for (from_node, to_node, supporting_reads) in edges {
        graph.add_edge(
            node_indices[&from_node],
            node_indices[&to_node],
            Edge {
                supporting_reads: serde_json::from_str(&supporting_reads)?,
            },
        );
    }
    Ok(VariantGraph {
        graph,
        start: start as i64,
        end: end as i64,
        target: target.clone(),
    })
}

/// Retrieves all variants from a given serialised graph file. Returns a HashMap of targets and all positions of variants in the graph.
pub(crate) fn variants_on_graph(path: &PathBuf) -> Result<HashMap<String, BTreeSet<i64>>> {
    let db = Connection::open(path)?;
    let mut stmt = db.prepare("SELECT target, pos FROM nodes")?;
    let variants: Vec<(String, i64)> = stmt
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))?
        .map(Result::unwrap)
        .collect();
    let mut variant_map = HashMap::new();
    for (target, pos) in variants {
        let target_variants = variant_map.entry(target).or_insert(BTreeSet::new());
        target_variants.insert(pos);
    }
    Ok(variant_map)
}

/// The distinct variants of a target in a serialised graph file, sorted by position.
pub(crate) fn variant_nodes(path: &Path, target: &str) -> Result<Vec<Node>> {
    let db = Connection::open(path)?;
    let mut stmt = db.prepare(
        "SELECT DISTINCT pos, reference_allele, alternative_allele FROM nodes \
         WHERE target = ? AND node_type = ? ORDER BY ALL",
    )?;
    Ok(stmt
        .query_map(params![target, NodeType::Variant.to_string()], |row| {
            Ok(Node::new(
                NodeType::Variant,
                row.get(0)?,
                row.get(1)?,
                row.get(2)?,
            ))
        })?
        .collect::<Result<_, _>>()?)
}

pub(crate) fn create_scores(output_path: &Path) -> Result<()> {
    let db = Connection::open(output_path)?;
    db.execute(
        "CREATE TABLE scores (transcript String, score FLOAT, frequencies String, hgvsc String, hgvsg String, hgvsg_full String, supporting_reads String, annotation String, protein String, consequence String)",
        [],
    )?;
    db.close().unwrap();
    Ok(())
}

pub(crate) fn write_scores(
    path: &Path,
    scores: Vec<HaplotypeScore>,
    transcript: Transcript,
) -> Result<()> {
    let mut db = Connection::open(path)?;
    let transaction = db.transaction()?;
    let mut stmt =
        transaction.prepare("INSERT INTO scores VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")?;
    let transcript_name = transcript.name();
    for (score, frequencies, supporting_reads, annotation) in scores {
        stmt.execute(params![
            transcript_name.as_str(),
            score.score(),
            json5::to_string(&frequencies)?,
            score.hgvsc,
            score.hgvsg,
            score.hgvsg_full,
            json5::to_string(&supporting_reads)?,
            json5::to_string(&annotation)?,
            score.altered_protein.to_string(),
            score.consequence.to_string(),
        ])?;
    }
    transaction.commit()?;
    db.close().unwrap();
    Ok(())
}

pub(crate) fn read_scores(
    path: &Path,
    notation: HgvsNotation,
) -> Result<HashMap<String, Vec<ScoreRecord>>> {
    let db = Connection::open(path)?;
    let mut scores = HashMap::new();
    let mut stmt = db.prepare(
        "SELECT transcript, score, frequencies, hgvsc, hgvsg, hgvsg_full, supporting_reads, annotation, protein FROM scores",
    )?;
    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        let transcript = row.get(0)?;
        let score = row.get(1)?;
        let frequencies: String = row.get(2)?;
        let hgvsc: String = row.get(3)?;
        let hgvsg: String = row.get(4)?;
        let hgvsg_full: String = row.get(5)?;
        let supporting_reads: String = row.get(6)?;
        let annotation: String = row.get(7)?;
        let protein: String = row.get(8)?;
        let haplotype = match notation {
            HgvsNotation::Hgvsc => hgvsc,
            HgvsNotation::Hgvsg => hgvsg,
            HgvsNotation::HgvsgFull => hgvsg_full,
        };
        scores.entry(transcript).or_insert(Vec::new()).push((
            score,
            json5::from_str(&frequencies)?,
            haplotype,
            json5::from_str(&supporting_reads)?,
            json5::from_str(&annotation)?,
            protein,
        ));
    }
    db.close().unwrap();
    Ok(scores)
}

/// Reads every `(transcript, haplotype)` row required to annotate a VCF/BCF.
pub(crate) fn read_annotation_rows(path: &Path) -> Result<Vec<AnnotationInput>> {
    let db = Connection::open(path)?;
    let mut rows = Vec::new();
    let mut stmt = db.prepare(
        "SELECT transcript, score, consequence, hgvsc, hgvsg, annotation, frequencies FROM scores",
    )?;
    let mut result = stmt.query([])?;
    while let Some(row) = result.next()? {
        let annotation: String = row.get(5)?;
        let frequencies: String = row.get(6)?;
        rows.push(AnnotationInput {
            transcript: row.get(0)?,
            score: row.get(1)?,
            consequence: row.get(2)?,
            hgvsc: row.get(3)?,
            hgvsg: row.get(4)?,
            annotation: json5::from_str(&annotation)?,
            frequencies: json5::from_str(&frequencies)?,
        });
    }
    db.close().unwrap();
    Ok(rows)
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::annotation::Annotation;
    use crate::graph::node::{Node, NodeType};
    use crate::graph::paths::Cds;
    use crate::graph::score::{Consequence, EffectScore};
    use crate::graph::Edge;
    use crate::translation::amino_acids::{AminoAcid, Protein};
    use crate::translation::distance::DistanceMetric;
    use bio::bio_types::strand::Strand;
    use itertools::Itertools;

    use petgraph::{Directed, Graph};

    pub(crate) fn setup_graph() -> VariantGraph {
        let mut graph = Graph::<Node, Edge, Directed>::new();
        let node1 = graph.add_node(Node::new(
            NodeType::Variant,
            1,
            "T".to_string(),
            "A".to_string(),
        ));
        let node2 = graph.add_node(Node::new(
            NodeType::Reference,
            2,
            "".to_string(),
            "".to_string(),
        ));
        let node3 = graph.add_node(Node::new(
            NodeType::Variant,
            3,
            "A".to_string(),
            "T".to_string(),
        ));
        let node4 = graph.add_node(Node::new(
            NodeType::Variant,
            4,
            "T".to_string(),
            "".to_string(),
        ));
        let _node5 = graph.add_node(Node::new(
            NodeType::Variant,
            8,
            "C".to_string(),
            "A".to_string(),
        ));
        let node6 = graph.add_node(Node::new(
            NodeType::Variant,
            9,
            "A".to_string(),
            "TT".to_string(),
        ));
        let _edge1 = graph.add_edge(
            node1,
            node2,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        let _edge2 = graph.add_edge(
            node2,
            node3,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        let _edge3 = graph.add_edge(
            node3,
            node4,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        let _edge4 = graph.add_edge(
            node4,
            node6,
            Edge {
                supporting_reads: HashMap::new(),
            },
        );
        VariantGraph {
            graph,
            start: 0,
            end: 10,
            target: "test".to_string(),
        }
    }

    #[test]
    fn test_write_graphs() {
        let mut graphs = HashMap::new();
        graphs.insert("graph1".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("graphs.duckdb");
        write_graphs(graphs, output_path.as_path()).unwrap();
        let db = Connection::open(output_path.as_path()).unwrap();
        let mut stmt = db.prepare("SELECT target FROM graphs").unwrap();
        let targets: Vec<String> = stmt
            .query_map([], |row| row.get(0))
            .unwrap()
            .map(Result::unwrap)
            .collect();
        assert_eq!(targets, vec!["graph1"]);
        let mut stmt = db
            .prepare("SELECT node_index, node_type, vaf, probs, pos, index FROM nodes")
            .unwrap();
        let nodes: Vec<(i64, String, String, String, i64, i64)> = stmt
            .query_map([], |row| {
                Ok((
                    row.get(0)?,
                    row.get(1)?,
                    row.get(2)?,
                    row.get(3)?,
                    row.get(4)?,
                    row.get(5)?,
                ))
            })
            .unwrap()
            .map(Result::unwrap)
            .collect();
        assert_eq!(nodes.len(), 6);
        let mut stmt = db
            .prepare("SELECT from_node, to_node, supporting_reads FROM edges")
            .unwrap();
        let edges: Vec<(i64, i64, String)> = stmt
            .query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))
            .unwrap()
            .map(Result::unwrap)
            .collect();
        assert_eq!(edges.len(), 4);
        db.close().unwrap();
    }

    #[test]
    fn test_write_graph_with_multiple_targets() {
        let mut graphs = HashMap::new();
        graphs.insert("graph1".to_string(), setup_graph());
        graphs.insert("graph2".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("graphs.duckdb");
        write_graphs(graphs, output_path.as_path()).unwrap();
        let db = Connection::open(output_path.as_path()).unwrap();
        let mut stmt = db.prepare("SELECT target FROM graphs").unwrap();
        let targets: Vec<String> = stmt
            .query_map([], |row| row.get(0))
            .unwrap()
            .map(Result::unwrap)
            .sorted_by(|a: &String, b| a.cmp(b))
            .collect();
        assert_eq!(targets, vec!["graph1", "graph2"]);
    }

    #[test]
    fn test_feature_graphs() {
        let mut graphs = HashMap::new();
        graphs.insert("graph1".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("graphs.duckdb");
        write_graphs(graphs, output_path.as_path()).unwrap();
        let graph = feature_graph(output_path, "graph1".to_string(), 1, 3).unwrap();
        assert_eq!(graph.graph.node_count(), 3);
        assert_eq!(graph.graph.edge_count(), 2);
        assert_eq!(graph.start, 1);
        assert_eq!(graph.end, 3);
    }

    #[test]
    fn test_feature_graphs_2() {
        let mut graphs = HashMap::new();
        graphs.insert("graph1".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("graphs.duckdb");
        write_graphs(graphs, output_path.as_path()).unwrap();
        let graph = feature_graph(output_path, "graph1".to_string(), 3, 10).unwrap();
        assert_eq!(graph.graph.node_count(), 4);
        assert_eq!(graph.graph.edge_count(), 2);
        dbg!(&graph);
        assert_eq!(graph.start, 3);
        assert_eq!(graph.end, 10);
    }

    fn write(graphs: Vec<(&str, VariantGraph)>) -> (tempfile::TempDir, PathBuf) {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("graphs.duckdb");
        let graphs = graphs
            .into_iter()
            .map(|(target, graph)| (target.to_string(), graph))
            .collect();
        write_graphs(graphs, &path).unwrap();
        (temp_dir, path)
    }

    /// Nodes keyed by position, type and alleles, and edges keyed by their endpoints.
    type Content = (
        Vec<(i64, String, String, String, u32)>,
        Vec<(
            (i64, String, String),
            (i64, String, String),
            Vec<(String, u32)>,
        )>,
    );

    fn content(graph: &VariantGraph) -> Content {
        let key = |node: &Node| {
            (
                node.pos,
                node.node_type.to_string(),
                node.alternative_allele.clone(),
            )
        };
        let nodes = graph
            .graph
            .node_weights()
            .map(|node| {
                (
                    node.pos,
                    node.node_type.to_string(),
                    node.reference_allele.clone(),
                    node.alternative_allele.clone(),
                    node.index,
                )
            })
            .sorted()
            .collect();
        let edges = graph
            .graph
            .edge_indices()
            .map(|edge| {
                let (from, to) = graph.graph.edge_endpoints(edge).unwrap();
                let reads = graph.graph[edge]
                    .supporting_reads
                    .iter()
                    .map(|(sample, reads)| (sample.clone(), *reads))
                    .sorted()
                    .collect();
                (key(&graph.graph[from]), key(&graph.graph[to]), reads)
            })
            .sorted()
            .collect();
        (nodes, edges)
    }

    #[test]
    fn feature_graph_drops_edges_leaving_the_range() {
        let (_dir, path) = write(vec![("graph1", setup_graph())]);
        let graph = feature_graph(path, "graph1".to_string(), 3, 8).unwrap();
        assert_eq!(graph.graph.node_count(), 3);
        assert_eq!(graph.graph.edge_count(), 1);
    }

    #[test]
    fn feature_graph_of_the_whole_range_restores_the_built_graph() {
        let observations = vec![crate::cli::ObservationFile {
            path: PathBuf::from("tests/resources/test_observations.vcf"),
            sample: "sample".to_string(),
        }];
        let built = VariantGraph::build(
            &PathBuf::from("tests/resources/test_calls.vcf"),
            &observations,
            "OX512233.1",
            bio::stats::LogProb::from(bio::stats::Prob(0.0)),
            0.0,
        )
        .unwrap();
        assert!(built.graph.edge_count() > 0);
        let expected = content(&built);
        let (_dir, path) = write(vec![("OX512233.1", built)]);
        let graph = feature_graph(path, "OX512233.1".to_string(), 0, 30_000).unwrap();
        assert_eq!(content(&graph), expected);
    }

    #[test]
    fn feature_graph_ignores_edges_of_other_targets() {
        let (_dir, path) = write(vec![("graph1", setup_graph()), ("graph2", setup_graph())]);
        let graph = feature_graph(path, "graph1".to_string(), 0, 10).unwrap();
        assert_eq!(graph.graph.node_count(), 6);
        assert_eq!(graph.graph.edge_count(), 4);
        assert_eq!(content(&graph), content(&setup_graph()));
    }

    #[test]
    fn feature_graph_keeps_the_stored_node_order() {
        let (_dir, path) = write(vec![("graph1", setup_graph())]);
        let graph = feature_graph(path, "graph1".to_string(), 0, 10).unwrap();
        let positions = graph
            .graph
            .node_weights()
            .map(|node| node.pos)
            .collect_vec();
        assert_eq!(positions, vec![1, 2, 3, 4, 8, 9]);
    }

    #[test]
    fn feature_graph_of_a_range_without_variants_is_empty() {
        let (_dir, path) = write(vec![("graph1", setup_graph())]);
        let graph = feature_graph(path, "graph1".to_string(), 5, 7).unwrap();
        assert_eq!(graph.graph.node_count(), 0);
        assert_eq!(graph.graph.edge_count(), 0);
    }

    #[test]
    fn test_variants_on_graph() {
        let mut graphs = HashMap::new();
        graphs.insert("chr1".to_string(), setup_graph());
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("graphs.duckdb");
        write_graphs(graphs, output_path.as_path()).unwrap();
        let result = variants_on_graph(&output_path).unwrap();
        let chr1_variants = result.get("chr1").unwrap();
        assert_eq!(chr1_variants, &vec![1, 2, 3, 4, 8, 9].into_iter().collect());
    }

    #[test]
    fn test_scores() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("scores.duckdb");
        create_scores(output_path.as_path()).unwrap();
        assert!(output_path.exists());
        let transcript = Transcript::new(
            "some feature".to_string(),
            "chr1".to_string(),
            Strand::Forward,
            vec![Cds::new(0, 100, 0)],
        );
        let p1 = Protein::new(vec![AminoAcid::Phenylalanine, AminoAcid::Leucine]);
        let p2 = Protein::new(vec![AminoAcid::Isoleucine, AminoAcid::Leucine]);
        let effect_score = EffectScore {
            original_protein: p1,
            altered_protein: p2,
            distance_metric: DistanceMetric::Epstein,
            realign: false,
            consequence: Consequence::Missense,
            hgvsc: "c.[100A>G;105C>T]".to_string(),
            hgvsg: "g.[100A>G;105C>T]".to_string(),
            hgvsg_full: "g.[100A>G;105C>T]".to_string(),
        };
        let frequencies = HashMap::from([("A".to_string(), 0.1), ("C".to_string(), 0.2)]);
        let supporting_reads = vec![HashMap::from([("A".to_string(), 10), ("C".to_string(), 5)])];
        let annotion = Annotation {
            revel_score: Some(0.8),
            spliceai_score: Some(0.7),
            alphamissense_score: Some(0.6),
            gnomad_frequencies: HashMap::new(),
        };
        let scores = vec![(effect_score, frequencies, supporting_reads, annotion)];
        write_scores(output_path.as_path(), scores, transcript).unwrap();
        let scores = read_scores(output_path.as_path(), HgvsNotation::Hgvsc).unwrap();
        assert_eq!(scores.len(), 1);
        assert!((scores.get("chr1:some feature").unwrap()[0].0 - 0.03999999910593033).abs() < 1e-6);
        assert!(
            (scores.get("chr1:some feature").unwrap()[0]
                .1
                .get("A")
                .unwrap()
                - 0.1)
                .abs()
                < 1e-6
        );
        assert!(
            (scores.get("chr1:some feature").unwrap()[0]
                .1
                .get("C")
                .unwrap()
                - 0.2)
                .abs()
                < 1e-6
        );
        assert_eq!(
            scores.get("chr1:some feature").unwrap()[0].3,
            vec![HashMap::from([
                ("A".to_string(), 10u32),
                ("C".to_string(), 5u32)
            ])]
        );
        assert_eq!(scores.get("chr1:some feature").unwrap()[0].5, "IL");
    }
}
