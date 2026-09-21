use crate::annotation::Annotation;
use crate::graph::node::Node;
use crate::graph::transcript::Transcript;
use crate::translation::amino_acids::AminoAcid;
use crate::translation::amino_acids::Protein;
use crate::translation::distance::DistanceMetric;
use anyhow::Result;
use bio::alignment::pairwise::Aligner;
use bio::alignment::AlignmentOperation::*;
use bio::bio_types::strand::Strand;
use clap::ValueEnum;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Coarse, Sequence-Ontology-style classification of a haplotype's effect on the
/// protein. Derived best-effort from the reference and altered protein sequences;
/// the graded [`EffectScore::score`] remains the primary product.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Consequence {
    Synonymous,
    Missense,
    StopGained,
    StopLost,
    StartLost,
    Frameshift,
    ProteinAltering,
}

impl Consequence {
    fn classify(
        original: &Protein,
        altered: &Protein,
        net_frameshift: i64,
        has_coding_indel: bool,
        start_lost: bool,
    ) -> Consequence {
        if start_lost {
            return Consequence::StartLost;
        }
        if net_frameshift % 3 != 0 {
            return Consequence::Frameshift;
        }
        // An in-frame change that introduces a stop upstream of the natural one is a
        // stop gain, whether or not it also rearranges residues (in-frame indel).
        if has_premature_stop(altered) && !has_premature_stop(original) {
            return Consequence::StopGained;
        }
        if net_frameshift != 0 {
            return Consequence::ProteinAltering;
        }
        if original == altered {
            return Consequence::Synonymous;
        }
        // Both proteins have the same length here, so a lost stop is the only
        // remaining categorical change; anything else is a plain substitution.
        for (reference_aa, altered_aa) in original
            .amino_acids()
            .into_iter()
            .zip(altered.amino_acids())
        {
            if reference_aa.is_stop() && !altered_aa.is_stop() {
                return Consequence::StopLost;
            }
        }
        // Offsetting coding indels keep the length unchanged but rewrite the
        // residues between them, which is structural rather than a substitution.
        if has_coding_indel {
            return Consequence::ProteinAltering;
        }
        Consequence::Missense
    }
}

/// True if the protein is truncated by a stop codon upstream of its last residue.
fn has_premature_stop(protein: &Protein) -> bool {
    protein
        .amino_acids()
        .iter()
        .rev()
        .skip(1)
        .any(AminoAcid::is_stop)
}

impl fmt::Display for Consequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let term = match self {
            Consequence::Synonymous => "synonymous_variant",
            Consequence::Missense => "missense_variant",
            Consequence::StopGained => "stop_gained",
            Consequence::StopLost => "stop_lost",
            Consequence::StartLost => "start_lost",
            Consequence::Frameshift => "frameshift_variant",
            Consequence::ProteinAltering => "protein_altering_variant",
        };
        write!(f, "{term}")
    }
}

#[derive(Debug, Clone)]
pub struct EffectScore {
    pub original_protein: Protein,
    pub altered_protein: Protein,
    pub distance_metric: DistanceMetric,
    pub realign: bool,
    pub consequence: Consequence,
    pub hgvsc: String,
    pub hgvsg: String,
    pub hgvsg_full: String,
}

impl EffectScore {
    pub(crate) fn from_haplotype(
        reference: &HashMap<String, Vec<u8>>,
        transcript: &Transcript,
        haplotype: &[Node],
        original_protein: Protein,
        distance_metric: DistanceMetric,
        mut realign: bool,
    ) -> Result<Self> {
        let mut altered_protein = Protein::from_haplotype(reference, transcript, haplotype)?;
        let start_lost = altered_protein.start_lost(&original_protein);
        if start_lost {
            altered_protein.apply_start_lost();
            realign = true;
        }
        let net_frameshift: i64 = haplotype.iter().map(Node::frameshift).sum();
        let has_coding_indel = haplotype
            .iter()
            .map(Node::frameshift)
            .any(|shift| shift != 0);
        let consequence = Consequence::classify(
            &original_protein,
            &altered_protein,
            net_frameshift,
            has_coding_indel,
            start_lost,
        );
        let mut variants: Vec<_> = haplotype
            .iter()
            .filter(|n| n.node_type.is_variant())
            .collect();
        let hgvsg = format!(
            "g.[{}]",
            variants
                .iter()
                .map(|n| n.hgvsg_token())
                .collect::<Vec<_>>()
                .join(";")
        );
        let hgvsg_full = format!(
            "g.[{}]",
            haplotype
                .iter()
                .map(|n| {
                    if n.node_type.is_variant() {
                        n.hgvsg_token()
                    } else {
                        format!("{}=", n.pos + 1)
                    }
                })
                .collect::<Vec<_>>()
                .join(";")
        );
        if transcript.strand == Strand::Reverse {
            variants.reverse();
        }
        let hgvsc = format!(
            "c.[{}]",
            variants
                .iter()
                .map(|n| n.hgvs_notation(transcript))
                .collect::<Vec<_>>()
                .join(";")
        );
        Ok(Self {
            original_protein,
            altered_protein,
            distance_metric,
            realign,
            consequence,
            hgvsc,
            hgvsg,
            hgvsg_full,
        })
    }

    pub fn score(&self) -> f64 {
        // Compare proteins:
        // If equal return 0
        // If nothing is translated (e.g. a lost start codon without a downstream start) return the maximum penalty
        // If no frameshift occured in the changed protein, simply compare amino acid by amino acid based on the distance metric and divide by the length of the protein
        // If a frameshift occurs, re-align the proteins and compare amino acid by amino acid based on the distance metric and divide by the length of the protein
        if self.original_protein == self.altered_protein {
            0.0
        } else if self.altered_protein.amino_acids().is_empty() {
            1.0
        } else if !self.realign {
            // We can assume both proteins have the same length
            let mut total = 0.0;
            for (i, (ref_aa, alt_aa)) in self
                .original_protein
                .amino_acids()
                .into_iter()
                .zip(self.altered_protein.amino_acids())
                .enumerate()
            {
                if alt_aa.is_stop() && !ref_aa.is_stop() {
                    total += (self.altered_protein.amino_acids().len() - i - 1) as f64;
                }
                total += self.distance_metric.compute(&ref_aa, &alt_aa);
            }
            total / self.altered_protein.amino_acids().len() as f64
        } else {
            // Realign using semiglobal
            let prot_x = self.original_protein.as_bytes();
            let prot_y = self.altered_protein.as_bytes();

            let score_fn = |a: u8, b: u8| {
                // Convert used distance matrix to a negative cost with integer scaling as expected by bio
                let d = self
                    .distance_metric
                    .compute(&AminoAcid::from(a), &AminoAcid::from(b));
                (-d * 10.0).round() as i32
            };

            let gap_open = -8;
            let gap_extend = -1;

            let mut aligner =
                Aligner::with_capacity(prot_x.len(), prot_y.len(), gap_open, gap_extend, &score_fn);

            let alignment = aligner.semiglobal(&prot_x, &prot_y);

            let ops = alignment.operations;

            let mut total = 0.0;

            let mut i = alignment.xstart;
            let mut j = alignment.ystart;

            for op in ops {
                match op {
                    Match => {
                        i += 1;
                        j += 1;
                    }
                    Subst => {
                        let aa_x_opt = prot_x.get(i).map(|&b| AminoAcid::from(b));
                        let aa_y_opt = prot_y.get(j).map(|&b| AminoAcid::from(b));

                        match (aa_x_opt, aa_y_opt) {
                            (Some(aa_x), Some(aa_y)) => {
                                if aa_y.is_stop() && !aa_x.is_stop() {
                                    let remaining = (prot_x.len() - i).max(prot_y.len() - j);
                                    total += remaining as f64;
                                    break;
                                }
                                total += self.distance_metric.compute(&aa_x, &aa_y);
                            }
                            _ => {
                                total += 1.0;
                            }
                        }

                        i += 1;
                        j += 1;
                    }
                    Del | Ins => {
                        total += 1.0;

                        if let Del = op {
                            i += 1;
                        } else {
                            j += 1;
                        }
                    }
                    _ => {}
                }
            }
            total / self.altered_protein.amino_acids().len() as f64
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, ValueEnum)]
pub enum HaplotypeMetric {
    Product,
    GeometricMean,
    Minimum,
}

pub type HaplotypeFrequency = HashMap<String, f32>;

pub(crate) type HaplotypeScore = (
    EffectScore,
    HaplotypeFrequency,
    Vec<HashMap<String, u32>>,
    Annotation,
);

pub(crate) type ScoreRecord = (
    f64,
    HaplotypeFrequency,
    String,
    Vec<HashMap<String, u32>>,
    Annotation,
    String,
);

/// A single `(transcript, haplotype)` row of `scores.duckdb` reduced to the fields
/// the `annotate` command needs to write predictions back into a VCF/BCF.
pub(crate) struct AnnotationInput {
    pub(crate) transcript: String,
    pub(crate) score: f64,
    pub(crate) consequence: String,
    pub(crate) hgvsc: String,
    pub(crate) hgvsg: String,
    pub(crate) annotation: Annotation,
    pub(crate) frequencies: HaplotypeFrequency,
}

impl HaplotypeMetric {
    pub fn calculate(&self, haplotype: &[Node], samples: &HashSet<String>) -> HaplotypeFrequency {
        let mut metrics = HashMap::new();
        for sample in samples {
            let vafs = haplotype
                .iter()
                .filter(|n| n.node_type.is_variant())
                .map(|n| *n.vaf.get(sample).unwrap_or(&0.0))
                .collect::<Vec<_>>();
            let result = match self {
                HaplotypeMetric::Product => vafs.iter().product(),
                HaplotypeMetric::GeometricMean => {
                    let product: f32 = vafs.iter().product();
                    product.powf(1.0 / vafs.len() as f32)
                }
                HaplotypeMetric::Minimum => vafs.iter().cloned().fold(1.0, f32::min),
            };
            metrics.insert(sample.to_string(), result);
        }
        metrics
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::node::NodeType;
    use crate::graph::EventProbs;
    use crate::translation::amino_acids::AminoAcid;

    #[test]
    fn test_product_metric() {
        let nodes = vec![
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.5)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 0,
                index: 0,
            },
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.25)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 1,
                index: 1,
            },
            Node {
                node_type: NodeType::Reference,
                reference_allele: "".to_string(),
                alternative_allele: "".to_string(),
                vaf: [("S1".to_string(), 0.75)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 2,
                index: 2,
            },
        ];
        let metric = HaplotypeMetric::Product;
        let result = metric.calculate(&nodes, &HashSet::from(["S1".to_string()]));
        assert_eq!(result.get("S1").unwrap(), &(0.5 * 0.25));
    }

    #[test]
    fn test_geometric_mean_metric() {
        let nodes = vec![
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.5)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 0,
                index: 0,
            },
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.25)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 1,
                index: 1,
            },
        ];
        let metric = HaplotypeMetric::GeometricMean;
        let result = metric.calculate(&nodes, &HashSet::from(["S1".to_string()]));
        let expected = (0.5_f32 * 0.25_f32).powf(1.0 / 2.0);
        assert!((result.get("S1").unwrap() - expected).abs() < 1e-6);
    }

    #[test]
    fn test_all_metrics_return_one_with_only_reference_nodes() {
        let nodes = vec![Node {
            node_type: NodeType::Reference,
            reference_allele: "".to_string(),
            alternative_allele: "".to_string(),
            vaf: [("S1".to_string(), 0.5)].into(),
            probs: EventProbs(HashMap::new()),
            pos: 0,
            index: 0,
        }];
        let metric = HaplotypeMetric::GeometricMean;
        let result = metric.calculate(&nodes, &HashSet::from(["S1".to_string()]));
        assert_eq!(result.get("S1").unwrap(), &1.0);
        let metric = HaplotypeMetric::Product;
        let result = metric.calculate(&nodes, &HashSet::from(["S1".to_string()]));
        assert_eq!(result.get("S1").unwrap(), &1.0);
        let metric = HaplotypeMetric::Minimum;
        let result = metric.calculate(&nodes, &HashSet::from(["S1".to_string()]));
        assert_eq!(result.get("S1").unwrap(), &1.0);
    }

    #[test]
    fn reference_haplotype_scores_one_for_every_sample() {
        let haplotype: Vec<Node> = Vec::new();
        let samples = HashSet::from(["S1".to_string(), "S2".to_string()]);
        for metric in [
            HaplotypeMetric::Minimum,
            HaplotypeMetric::Product,
            HaplotypeMetric::GeometricMean,
        ] {
            let result = metric.calculate(&haplotype, &samples);
            assert_eq!(result.get("S1"), Some(&1.0));
            assert_eq!(result.get("S2"), Some(&1.0));
        }
    }

    #[test]
    fn test_minimum_metric() {
        let nodes = vec![
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.5)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 0,
                index: 0,
            },
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.25)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 1,
                index: 1,
            },
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "A".to_string(),
                vaf: [("S1".to_string(), 0.75)].into(),
                probs: EventProbs(HashMap::new()),
                pos: 2,
                index: 2,
            },
        ];
        let metric = HaplotypeMetric::Minimum;
        let result = metric.calculate(&nodes, &HashSet::from(["S1".to_string()]));
        assert_eq!(*result.get("S1").unwrap(), 0.25);
    }

    #[test]
    fn test_score_equal_proteins() {
        let p1 = Protein::new(vec![AminoAcid::Phenylalanine, AminoAcid::Leucine]);
        let p2 = Protein::new(vec![AminoAcid::Phenylalanine, AminoAcid::Leucine]);
        let score = EffectScore {
            original_protein: p1,
            altered_protein: p2,
            distance_metric: DistanceMetric::Epstein,
            realign: false,
            consequence: Consequence::Missense,
            hgvsc: "c.[100A>G;105C>T]".to_string(),
            hgvsg: "g.[100A>G;105C>T]".to_string(),
            hgvsg_full: "g.[100A>G;105C>T]".to_string(),
        };
        assert!(score.score().abs() < 1e-6);
    }

    #[test]
    fn test_score_different_proteins() {
        let p1 = Protein::new(vec![AminoAcid::Phenylalanine, AminoAcid::Leucine]);
        let p2 = Protein::new(vec![AminoAcid::Phenylalanine, AminoAcid::Valine]);
        let score = EffectScore {
            original_protein: p1,
            altered_protein: p2,
            distance_metric: DistanceMetric::Epstein,
            realign: false,
            consequence: Consequence::Missense,
            hgvsc: "c.[100A>G;105C>T]".to_string(),
            hgvsg: "g.[100A>G;105C>T]".to_string(),
            hgvsg_full: "g.[100A>G;105C>T]".to_string(),
        };
        // Distance between Leucine and Valine is 0.03 / 2 since len = 2
        assert!((score.score() - 0.015).abs() < 1e-6)
    }

    #[test]
    fn test_score_different_proteins_with_stop() {
        let p1 = Protein::new(vec![
            AminoAcid::Phenylalanine,
            AminoAcid::Leucine,
            AminoAcid::Valine,
        ]);
        let p2 = Protein::new(vec![
            AminoAcid::Phenylalanine,
            AminoAcid::Stop,
            AminoAcid::Valine,
        ]);
        let score = EffectScore {
            original_protein: p1,
            altered_protein: p2,
            distance_metric: DistanceMetric::Epstein,
            realign: false,
            consequence: Consequence::Missense,
            hgvsc: "c.[100A>G;105C>T]".to_string(),
            hgvsg: "g.[100A>G;105C>T]".to_string(),
            hgvsg_full: "g.[100A>G;105C>T]".to_string(),
        };
        assert!((score.score() - (2.0 / 3.0)).abs() < 1e-6)
    }

    #[test]
    fn test_score_different_proteins_realign() {
        let p1 = Protein::new(vec![
            AminoAcid::Phenylalanine,
            AminoAcid::Leucine,
            AminoAcid::Valine,
        ]);
        let p2 = Protein::new(vec![AminoAcid::Leucine, AminoAcid::Valine]);
        let score = EffectScore {
            original_protein: p1,
            altered_protein: p2,
            distance_metric: DistanceMetric::Epstein,
            realign: true,
            consequence: Consequence::Frameshift,
            hgvsc: "c.[100A>G;105C>T]".to_string(),
            hgvsg: "g.[100A>G;105C>T]".to_string(),
            hgvsg_full: "g.[100A>G;105C>T]".to_string(),
        };
        assert!((score.score() - 0.5).abs() < 1e-6)
    }

    #[test]
    fn test_score_different_proteins_realign_with_substitution() {
        let p1 = Protein::new(vec![
            AminoAcid::Phenylalanine,
            AminoAcid::Phenylalanine,
            AminoAcid::Leucine,
            AminoAcid::Phenylalanine,
            AminoAcid::Phenylalanine,
            AminoAcid::Phenylalanine,
        ]);
        let p2 = Protein::new(vec![
            AminoAcid::Phenylalanine,
            AminoAcid::Valine,
            AminoAcid::Phenylalanine,
            AminoAcid::Phenylalanine,
            AminoAcid::Phenylalanine,
        ]);
        let score = EffectScore {
            original_protein: p1,
            altered_protein: p2,
            distance_metric: DistanceMetric::Epstein,
            realign: true,
            consequence: Consequence::Frameshift,
            hgvsc: "c.[100A>G;105C>T]".to_string(),
            hgvsg: "g.[100A>G;105C>T]".to_string(),
            hgvsg_full: "g.[100A>G;105C>T]".to_string(),
        };
        assert!((score.score() - 0.2).abs() < 1e-6)
    }

    #[test]
    fn classify_identical_proteins_is_synonymous() {
        let protein = Protein::new(vec![AminoAcid::Methionine, AminoAcid::Leucine]);
        let consequence = Consequence::classify(&protein, &protein.clone(), 0, false, false);
        assert_eq!(consequence, Consequence::Synonymous);
    }

    #[test]
    fn classify_single_substitution_is_missense() {
        let original = Protein::new(vec![AminoAcid::Methionine, AminoAcid::Leucine]);
        let altered = Protein::new(vec![AminoAcid::Methionine, AminoAcid::Valine]);
        assert_eq!(
            Consequence::classify(&original, &altered, 0, false, false),
            Consequence::Missense
        );
    }

    #[test]
    fn classify_premature_stop_is_stop_gained() {
        let original = Protein::new(vec![
            AminoAcid::Methionine,
            AminoAcid::Leucine,
            AminoAcid::Valine,
        ]);
        let altered = Protein::new(vec![
            AminoAcid::Methionine,
            AminoAcid::Stop,
            AminoAcid::Valine,
        ]);
        assert_eq!(
            Consequence::classify(&original, &altered, 0, false, false),
            Consequence::StopGained
        );
    }

    #[test]
    fn classify_inframe_indel_that_introduces_a_stop_is_stop_gained() {
        let original = Protein::new(vec![
            AminoAcid::Methionine,
            AminoAcid::Lysine,
            AminoAcid::Leucine,
            AminoAcid::Valine,
            AminoAcid::Stop,
        ]);
        let altered = Protein::new(vec![
            AminoAcid::Methionine,
            AminoAcid::Stop,
            AminoAcid::Valine,
            AminoAcid::Stop,
        ]);
        assert_eq!(
            Consequence::classify(&original, &altered, -3, true, false),
            Consequence::StopGained
        );
    }

    #[test]
    fn classify_inframe_deletion_without_new_stop_is_protein_altering() {
        let original = Protein::new(vec![
            AminoAcid::Methionine,
            AminoAcid::Lysine,
            AminoAcid::Leucine,
            AminoAcid::Valine,
            AminoAcid::Stop,
        ]);
        let altered = Protein::new(vec![
            AminoAcid::Methionine,
            AminoAcid::Leucine,
            AminoAcid::Valine,
            AminoAcid::Stop,
        ]);
        assert_eq!(
            Consequence::classify(&original, &altered, -3, true, false),
            Consequence::ProteinAltering
        );
    }

    #[test]
    fn classify_offsetting_coding_indels_is_protein_altering() {
        let original = Protein::new(vec![
            AminoAcid::Methionine,
            AminoAcid::Lysine,
            AminoAcid::Leucine,
            AminoAcid::Valine,
        ]);
        let altered = Protein::new(vec![
            AminoAcid::Methionine,
            AminoAcid::Threonine,
            AminoAcid::Proline,
            AminoAcid::Valine,
        ]);
        // An insertion and a deletion of equal length cancel to a net frameshift
        // of zero while still rewriting the residues between them.
        assert_eq!(
            Consequence::classify(&original, &altered, 0, true, false),
            Consequence::ProteinAltering
        );
    }

    #[test]
    fn classify_lost_stop_is_stop_lost() {
        let original = Protein::new(vec![AminoAcid::Methionine, AminoAcid::Stop]);
        let altered = Protein::new(vec![AminoAcid::Methionine, AminoAcid::Valine]);
        assert_eq!(
            Consequence::classify(&original, &altered, 0, false, false),
            Consequence::StopLost
        );
    }

    #[test]
    fn classify_uses_frameshift_and_start_lost_before_residue_comparison() {
        let original = Protein::new(vec![AminoAcid::Methionine, AminoAcid::Leucine]);
        let altered = Protein::new(vec![AminoAcid::Methionine, AminoAcid::Valine]);
        assert_eq!(
            Consequence::classify(&original, &altered, 1, true, false),
            Consequence::Frameshift
        );
        assert_eq!(
            Consequence::classify(&original, &altered, 3, true, false),
            Consequence::ProteinAltering
        );
        assert_eq!(
            Consequence::classify(&original, &altered, 0, false, true),
            Consequence::StartLost
        );
    }

    #[test]
    fn consequence_renders_sequence_ontology_terms() {
        assert_eq!(Consequence::Missense.to_string(), "missense_variant");
        assert_eq!(Consequence::StopGained.to_string(), "stop_gained");
        assert_eq!(Consequence::Frameshift.to_string(), "frameshift_variant");
    }
}
