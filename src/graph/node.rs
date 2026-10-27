use crate::graph::hgvs::hgvsc;
use crate::graph::transcript::Transcript;
use crate::graph::EventProbs;
use anyhow::anyhow;
use rust_htslib::bcf::Record;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::Display;
use std::str::FromStr;
use varlociraptor::variants::evidence::observations::read_observation::ProcessedReadObservation;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
#[allow(dead_code)] // TODO: Remove this attribute when graph is properly serialized
pub(crate) enum NodeType {
    Variant,
    Reference,
}

impl FromStr for NodeType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Variant" => Ok(NodeType::Variant),
            "Reference" => Ok(NodeType::Reference),
            _ => Err(anyhow!("Invalid node type: {}", s)),
        }
    }
}

impl Display for NodeType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NodeType::Variant => write!(f, "Variant"),
            NodeType::Reference => write!(f, "Reference"),
        }
    }
}

impl NodeType {
    pub(crate) fn is_variant(&self) -> bool {
        match self {
            NodeType::Variant => true,
            NodeType::Reference => false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(dead_code)] // TODO: Remove this attribute when graph is properly serialized
pub(crate) struct Node {
    pub(crate) node_type: NodeType,
    pub(crate) reference_allele: String,
    pub(crate) alternative_allele: String,
    pub(crate) vaf: HashMap<String, f32>,
    pub(crate) probs: EventProbs,
    pub(crate) pos: i64,
    pub(crate) index: u32,
}

impl std::fmt::Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} at position {}", self.node_type, self.pos)
    }
}

impl Node {
    pub(crate) fn from_records(
        calls_record: &Record,
        _observations_record: &HashMap<&&String, Vec<ProcessedReadObservation>>,
        event_probs: &EventProbs,
        node_type: NodeType,
        samples: &[String],
        index: u32,
    ) -> Self {
        let vafs = calls_record.format(b"AF").float().unwrap();
        let alleles = calls_record.alleles();
        let ref_allele = String::from_utf8(alleles[0].to_vec()).unwrap();
        let alt_allele = String::from_utf8(alleles[1].to_vec()).unwrap();
        Node {
            node_type,
            reference_allele: ref_allele,
            alternative_allele: alt_allele,
            vaf: samples
                .iter()
                .zip(vafs.iter())
                .map(|(s, v)| (s.to_string(), v[0]))
                .collect(),
            probs: event_probs.to_owned(),
            pos: calls_record.pos(),
            index,
        }
    }

    pub(crate) fn new(
        node_type: NodeType,
        pos: i64,
        reference_allele: String,
        alternative_allele: String,
    ) -> Self {
        Node {
            node_type,
            reference_allele,
            alternative_allele,
            vaf: Default::default(),
            probs: EventProbs(HashMap::new()),
            pos,
            index: 0,
        }
    }

    /// Returns the frameshift caused by the variant at this node.
    pub(crate) fn frameshift(&self) -> i64 {
        match &self.node_type {
            NodeType::Variant => {
                self.alternative_allele.len() as i64 - self.reference_allele.len() as i64
            }
            NodeType::Reference => 0,
        }
    }

    /// Returns the maximum VAF across all samples
    pub(crate) fn max_vaf(&self) -> f32 {
        *self
            .vaf
            .values()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(&0.0)
    }

    pub(crate) fn hgvs_notation(&self, transcript: &Transcript) -> String {
        hgvsc(
            transcript,
            self.pos as u64,
            &self.reference_allele,
            &self.alternative_allele,
        )
        .unwrap_or_default()
    }

    /// Renders the variant as a single genomic HGVS token, e.g. `467C>T`, using a
    /// one-based position. This is the token grammar the `annotate` command parses
    /// back out of a haplotype's stored `hgvsg`.
    pub(crate) fn hgvsg_token(&self) -> String {
        format!(
            "{}{}>{}",
            self.pos + 1,
            self.reference_allele,
            self.alternative_allele
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::paths::Cds;
    use bio::bio_types::strand::Strand;

    #[test]
    fn test_frameshift() {
        let node = Node::new(NodeType::Variant, 2, "A".to_string(), "A".to_string());
        let frameshift = node.frameshift();
        assert_eq!(frameshift, 0);
        let node = Node::new(NodeType::Variant, 2, "A".to_string(), "AT".to_string());
        let frameshift = node.frameshift();
        assert_eq!(frameshift, 1);
        let node = Node::new(NodeType::Variant, 2, "A".to_string(), "".to_string());
        let frameshift = node.frameshift();
        assert_eq!(frameshift, -1);
    }

    #[test]
    fn from_str_parses_variant_node_type() {
        let result = NodeType::from_str("Variant").unwrap();
        assert_eq!(result, NodeType::Variant);
    }

    #[test]
    fn from_str_parses_reference_node_type() {
        let result = NodeType::from_str("Reference").unwrap();
        assert_eq!(result, NodeType::Reference);
    }

    #[test]
    fn from_str_returns_error_for_invalid_node_type() {
        let result = NodeType::from_str("Invalid");
        assert!(result.is_err());
    }

    #[test]
    fn from_str_returns_error_for_empty_string() {
        let result = NodeType::from_str("");
        assert!(result.is_err());
    }

    #[test]
    fn test_node_display() {
        let var_node = Node {
            node_type: NodeType::Variant,
            reference_allele: "C".to_string(),
            alternative_allele: "A".to_string(),
            vaf: Default::default(),
            probs: EventProbs(Default::default()),
            pos: 42,
            index: 0,
        };
        let ref_node = Node {
            node_type: NodeType::Reference,
            reference_allele: "".to_string(),
            alternative_allele: "".to_string(),
            vaf: Default::default(),
            probs: EventProbs(Default::default()),
            pos: 99,
            index: 1,
        };

        assert_eq!(format!("{}", var_node), "Variant at position 42");
        assert_eq!(format!("{}", ref_node), "Reference at position 99");
    }

    #[test]
    fn test_node_max_vaf() {
        let mut vafs = HashMap::new();
        vafs.insert("A".to_string(), 0.05);
        vafs.insert("C".to_string(), 0.1);
        vafs.insert("G".to_string(), 0.02);
        vafs.insert("T".to_string(), 0.03);
        let var_node = Node {
            node_type: NodeType::Variant,
            reference_allele: "C".to_string(),
            alternative_allele: "A".to_string(),
            vaf: vafs.clone(),
            probs: EventProbs(Default::default()),
            pos: 42,
            index: 0,
        };
        assert_eq!(var_node.max_vaf(), 0.1);
    }

    #[test]
    fn test_node_hgvs() {
        let var_node = Node {
            node_type: NodeType::Variant,
            reference_allele: "C".to_string(),
            alternative_allele: "A".to_string(),
            vaf: Default::default(),
            probs: EventProbs(Default::default()),
            pos: 42,
            index: 0,
        };
        let transcript = Transcript {
            feature: "ENST00000367770.8".to_string(),
            target: "chr1".to_string(),
            strand: Strand::Forward,
            coding_sequences: vec![Cds::new(0, 100, 2), Cds::new(200, 300, 1)],
        };
        assert_eq!(var_node.hgvs_notation(&transcript), "43C>A");
    }

    #[test]
    fn test_node_hgvsg_token() {
        let node = Node::new(NodeType::Variant, 466, "C".to_string(), "T".to_string());
        assert_eq!(node.hgvsg_token(), "467C>T");
        let deletion = Node::new(NodeType::Variant, 4, "AT".to_string(), "A".to_string());
        assert_eq!(deletion.hgvsg_token(), "5AT>A");
    }
}
