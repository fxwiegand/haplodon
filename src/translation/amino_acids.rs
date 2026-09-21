use crate::graph::node::Node;
use crate::graph::transcript::Transcript;
use crate::translation::dna_to_amino_acids;
use anyhow::Result;
use bio::alphabets::dna;
use bio::bio_types::strand::Strand;
use std::collections::HashMap;
use std::fmt::Display;

/// A protein consisting of a sequence of amino acids
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Protein {
    pub(crate) sequence: Vec<AminoAcid>,
}

impl Protein {
    /// Creates a new protein from a sequence of amino acids
    pub(crate) fn new(sequence: Vec<AminoAcid>) -> Protein {
        Protein { sequence }
    }

    pub(crate) fn amino_acids(&self) -> Vec<AminoAcid> {
        self.sequence.clone()
    }

    pub(crate) fn from_transcript(
        reference: &HashMap<String, Vec<u8>>,
        transcript: &Transcript,
    ) -> Result<Protein> {
        let target_ref = reference
            .get(&transcript.target)
            .ok_or_else(|| anyhow::anyhow!("Missing reference for {}", transcript.target))?;

        let cds_seq: Vec<u8> = transcript
            .cds()
            .flat_map(|cds| {
                let region = &target_ref[cds.start as usize..=cds.end as usize];
                match transcript.strand {
                    Strand::Reverse => dna::revcomp(region),
                    _ => region.to_vec(),
                }
            })
            .collect();
        Ok(Protein::new(dna_to_amino_acids(&cds_seq)?))
    }

    pub(crate) fn from_haplotype(
        reference: &HashMap<String, Vec<u8>>,
        transcript: &Transcript,
        haplotype: &[Node],
    ) -> Result<Protein> {
        let target_ref = reference
            .get(&transcript.target)
            .ok_or_else(|| anyhow::anyhow!("Missing reference for {}", transcript.target))?;

        let cds_seq: Vec<u8> = transcript
            .cds()
            .flat_map(|cds| {
                let mut region = target_ref[cds.start as usize..=cds.end as usize].to_vec();
                let mut offset: isize = 0;

                let mut nodes: Vec<&Node> = haplotype
                    .iter()
                    .filter(|n| {
                        n.pos >= cds.start as i64
                            && n.pos <= cds.end as i64
                            && n.node_type.is_variant()
                    })
                    .collect();
                nodes.sort_by_key(|n| n.pos);

                for node in nodes {
                    let pos = ((node.pos - cds.start as i64) as isize + offset) as usize;
                    // Check if variant exceeds the boundaries of CDS, if so only use the proportion of the variant within the CDS
                    // This might affect splicing
                    let remaining_bases = region.len().saturating_sub(pos);
                    if remaining_bases == 0 {
                        continue; // nothing left in CDS
                    }
                    let ref_len = node.reference_allele.len().min(remaining_bases);
                    let alt_bytes = node.alternative_allele.as_bytes();
                    region.splice(pos..pos + ref_len, alt_bytes.iter().copied());
                    offset += alt_bytes.len() as isize - ref_len as isize;
                }

                match transcript.strand {
                    Strand::Reverse => dna::revcomp(&region),
                    _ => region.to_vec(),
                }
            })
            .collect();

        Ok(Protein::new(dna_to_amino_acids(&cds_seq)?))
    }

    /// Returns true if the protein has lost the start codon compared to the original protein sequence
    pub(crate) fn start_lost(&self, original: &Protein) -> bool {
        matches!(original.sequence.first(), Some(AminoAcid::Methionine))
            && !matches!(self.sequence.first(), Some(AminoAcid::Methionine))
    }

    /// Applies the start lost effect to the protein sequence by removing all amino acids before the next methionine if present, otherwise retrunss an empty sequence
    pub(crate) fn apply_start_lost(&mut self) {
        let reinitiated = match self
            .sequence
            .iter()
            .skip(1)
            .position(|amino_acid| *amino_acid == AminoAcid::Methionine)
        {
            Some(offset) => self.sequence[offset + 1..].to_vec(),
            None => Vec::new(),
        };
        self.sequence = reinitiated;
    }

    pub fn as_bytes(&self) -> Vec<u8> {
        self.amino_acids()
            .iter()
            .map(|aa| aa.short_abbreviation() as u8)
            .collect()
    }
}

impl Display for Protein {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", String::from_utf8(self.as_bytes()).unwrap())
    }
}

/// Amino acids and stop
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum AminoAcid {
    Alanine,
    Arginine,
    Asparagine,
    AsparticAcid,
    Cysteine,
    GlutamicAcid,
    Glutamine,
    Glycine,
    Histidine,
    Isoleucine,
    Leucine,
    Lysine,
    Methionine,
    Phenylalanine,
    Proline,
    Serine,
    Threonine,
    Tryptophan,
    Tyrosine,
    Valine,
    Stop,
}

impl AminoAcid {
    /// Converts a codon to an amino acid or stop codon. Returns an error if the codon is invalid.
    pub(crate) fn from_codon(codon: &[u8]) -> Result<AminoAcid> {
        match codon {
            b"GCU" | b"GCC" | b"GCA" | b"GCG" => Ok(AminoAcid::Alanine),
            b"CGU" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => Ok(AminoAcid::Arginine),
            b"AAU" | b"AAC" => Ok(AminoAcid::Asparagine),
            b"GAU" | b"GAC" => Ok(AminoAcid::AsparticAcid),
            b"UGU" | b"UGC" => Ok(AminoAcid::Cysteine),
            b"GAA" | b"GAG" => Ok(AminoAcid::GlutamicAcid),
            b"CAA" | b"CAG" => Ok(AminoAcid::Glutamine),
            b"GGU" | b"GGC" | b"GGA" | b"GGG" => Ok(AminoAcid::Glycine),
            b"CAU" | b"CAC" => Ok(AminoAcid::Histidine),
            b"AUU" | b"AUC" | b"AUA" => Ok(AminoAcid::Isoleucine),
            b"UUA" | b"UUG" | b"CUU" | b"CUC" | b"CUA" | b"CUG" => Ok(AminoAcid::Leucine),
            b"AAA" | b"AAG" => Ok(AminoAcid::Lysine),
            b"AUG" => Ok(AminoAcid::Methionine),
            b"UUU" | b"UUC" => Ok(AminoAcid::Phenylalanine),
            b"CCU" | b"CCC" | b"CCA" | b"CCG" => Ok(AminoAcid::Proline),
            b"UCU" | b"UCC" | b"UCA" | b"UCG" | b"AGU" | b"AGC" => Ok(AminoAcid::Serine),
            b"ACU" | b"ACC" | b"ACA" | b"ACG" => Ok(AminoAcid::Threonine),
            b"UGG" => Ok(AminoAcid::Tryptophan),
            b"UAU" | b"UAC" => Ok(AminoAcid::Tyrosine),
            b"GUU" | b"GUC" | b"GUA" | b"GUG" => Ok(AminoAcid::Valine),
            b"UAA" | b"UAG" | b"UGA" => Ok(AminoAcid::Stop),
            _ => anyhow::bail!("Invalid codon"),
        }
    }

    /// Returns the 3 letter abbreviation of the amino acid
    fn abbreviation(&self) -> &str {
        match self {
            AminoAcid::Alanine => "Ala",
            AminoAcid::Arginine => "Arg",
            AminoAcid::Asparagine => "Asn",
            AminoAcid::AsparticAcid => "Asp",
            AminoAcid::Cysteine => "Cys",
            AminoAcid::GlutamicAcid => "Glu",
            AminoAcid::Glutamine => "Gln",
            AminoAcid::Glycine => "Gly",
            AminoAcid::Histidine => "His",
            AminoAcid::Isoleucine => "Ile",
            AminoAcid::Leucine => "Leu",
            AminoAcid::Lysine => "Lys",
            AminoAcid::Methionine => "Met",
            AminoAcid::Phenylalanine => "Phe",
            AminoAcid::Proline => "Pro",
            AminoAcid::Serine => "Ser",
            AminoAcid::Threonine => "Thr",
            AminoAcid::Tryptophan => "Trp",
            AminoAcid::Tyrosine => "Tyr",
            AminoAcid::Valine => "Val",
            AminoAcid::Stop => "Stop",
        }
    }

    pub(crate) fn short_abbreviation(&self) -> char {
        match self {
            AminoAcid::Alanine => 'A',
            AminoAcid::Arginine => 'R',
            AminoAcid::Asparagine => 'N',
            AminoAcid::AsparticAcid => 'D',
            AminoAcid::Cysteine => 'C',
            AminoAcid::GlutamicAcid => 'E',
            AminoAcid::Glutamine => 'Q',
            AminoAcid::Glycine => 'G',
            AminoAcid::Histidine => 'H',
            AminoAcid::Isoleucine => 'I',
            AminoAcid::Leucine => 'L',
            AminoAcid::Lysine => 'K',
            AminoAcid::Methionine => 'M',
            AminoAcid::Phenylalanine => 'F',
            AminoAcid::Proline => 'P',
            AminoAcid::Serine => 'S',
            AminoAcid::Threonine => 'T',
            AminoAcid::Tryptophan => 'W',
            AminoAcid::Tyrosine => 'Y',
            AminoAcid::Valine => 'V',
            AminoAcid::Stop => 'X',
        }
    }

    /// Returns true if the amino acid is a stop codon
    pub(crate) fn is_stop(&self) -> bool {
        matches!(self, AminoAcid::Stop)
    }
}

impl Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.abbreviation())
    }
}

impl From<u8> for AminoAcid {
    fn from(byte: u8) -> Self {
        match byte {
            b'A' => AminoAcid::Alanine,
            b'R' => AminoAcid::Arginine,
            b'N' => AminoAcid::Asparagine,
            b'D' => AminoAcid::AsparticAcid,
            b'C' => AminoAcid::Cysteine,
            b'Q' => AminoAcid::Glutamine,
            b'E' => AminoAcid::GlutamicAcid,
            b'G' => AminoAcid::Glycine,
            b'H' => AminoAcid::Histidine,
            b'I' => AminoAcid::Isoleucine,
            b'L' => AminoAcid::Leucine,
            b'K' => AminoAcid::Lysine,
            b'M' => AminoAcid::Methionine,
            b'F' => AminoAcid::Phenylalanine,
            b'P' => AminoAcid::Proline,
            b'S' => AminoAcid::Serine,
            b'T' => AminoAcid::Threonine,
            b'W' => AminoAcid::Tryptophan,
            b'Y' => AminoAcid::Tyrosine,
            b'V' => AminoAcid::Valine,
            b'X' => AminoAcid::Stop,
            _ => panic!("Invalid amino acid byte"),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::graph::node::NodeType;
    use crate::graph::paths::Cds;
    use crate::graph::EventProbs;

    use super::*;

    #[test]
    fn test_valid_codons() {
        assert!(matches!(
            AminoAcid::from_codon(b"GCU"),
            Ok(AminoAcid::Alanine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"CGU"),
            Ok(AminoAcid::Arginine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"AAU"),
            Ok(AminoAcid::Asparagine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"GAU"),
            Ok(AminoAcid::AsparticAcid)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UGU"),
            Ok(AminoAcid::Cysteine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"GAA"),
            Ok(AminoAcid::GlutamicAcid)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"CAA"),
            Ok(AminoAcid::Glutamine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"GGU"),
            Ok(AminoAcid::Glycine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"CAU"),
            Ok(AminoAcid::Histidine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"AUU"),
            Ok(AminoAcid::Isoleucine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UUA"),
            Ok(AminoAcid::Leucine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"AAA"),
            Ok(AminoAcid::Lysine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"AUG"),
            Ok(AminoAcid::Methionine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UUU"),
            Ok(AminoAcid::Phenylalanine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"CCU"),
            Ok(AminoAcid::Proline)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UCU"),
            Ok(AminoAcid::Serine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"ACU"),
            Ok(AminoAcid::Threonine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UGG"),
            Ok(AminoAcid::Tryptophan)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"UAU"),
            Ok(AminoAcid::Tyrosine)
        ));
        assert!(matches!(
            AminoAcid::from_codon(b"GUU"),
            Ok(AminoAcid::Valine)
        ));
        assert!(matches!(AminoAcid::from_codon(b"UAA"), Ok(AminoAcid::Stop)));
        assert!(matches!(AminoAcid::from_codon(b"UAG"), Ok(AminoAcid::Stop)));
        assert!(matches!(AminoAcid::from_codon(b"UGA"), Ok(AminoAcid::Stop)));
    }

    #[test]
    fn test_invalid_codons() {
        assert!(AminoAcid::from_codon(b"XXX").is_err());
        assert!(AminoAcid::from_codon(b"XYZ").is_err());
        assert!(AminoAcid::from_codon(b"UGX").is_err());
        assert!(AminoAcid::from_codon(b"ACX").is_err());
    }

    #[test]
    fn test_abbreviations() {
        assert_eq!(AminoAcid::Alanine.abbreviation(), "Ala");
        assert_eq!(AminoAcid::Arginine.abbreviation(), "Arg");
        assert_eq!(AminoAcid::Asparagine.abbreviation(), "Asn");
        assert_eq!(AminoAcid::AsparticAcid.abbreviation(), "Asp");
        assert_eq!(AminoAcid::Cysteine.abbreviation(), "Cys");
        assert_eq!(AminoAcid::GlutamicAcid.abbreviation(), "Glu");
        assert_eq!(AminoAcid::Glutamine.abbreviation(), "Gln");
        assert_eq!(AminoAcid::Glycine.abbreviation(), "Gly");
        assert_eq!(AminoAcid::Histidine.abbreviation(), "His");
        assert_eq!(AminoAcid::Isoleucine.abbreviation(), "Ile");
        assert_eq!(AminoAcid::Leucine.abbreviation(), "Leu");
        assert_eq!(AminoAcid::Lysine.abbreviation(), "Lys");
        assert_eq!(AminoAcid::Methionine.abbreviation(), "Met");
        assert_eq!(AminoAcid::Phenylalanine.abbreviation(), "Phe");
        assert_eq!(AminoAcid::Proline.abbreviation(), "Pro");
        assert_eq!(AminoAcid::Serine.abbreviation(), "Ser");
        assert_eq!(AminoAcid::Threonine.abbreviation(), "Thr");
        assert_eq!(AminoAcid::Tryptophan.abbreviation(), "Trp");
        assert_eq!(AminoAcid::Tyrosine.abbreviation(), "Tyr");
        assert_eq!(AminoAcid::Valine.abbreviation(), "Val");
        assert_eq!(AminoAcid::Stop.abbreviation(), "Stop");
    }

    #[test]
    fn test_short_abbreviations() {
        assert_eq!(AminoAcid::Alanine.short_abbreviation(), 'A');
        assert_eq!(AminoAcid::Arginine.short_abbreviation(), 'R');
        assert_eq!(AminoAcid::Asparagine.short_abbreviation(), 'N');
        assert_eq!(AminoAcid::AsparticAcid.short_abbreviation(), 'D');
        assert_eq!(AminoAcid::Cysteine.short_abbreviation(), 'C');
        assert_eq!(AminoAcid::GlutamicAcid.short_abbreviation(), 'E');
        assert_eq!(AminoAcid::Glutamine.short_abbreviation(), 'Q');
        assert_eq!(AminoAcid::Glycine.short_abbreviation(), 'G');
        assert_eq!(AminoAcid::Histidine.short_abbreviation(), 'H');
        assert_eq!(AminoAcid::Isoleucine.short_abbreviation(), 'I');
        assert_eq!(AminoAcid::Leucine.short_abbreviation(), 'L');
        assert_eq!(AminoAcid::Lysine.short_abbreviation(), 'K');
        assert_eq!(AminoAcid::Methionine.short_abbreviation(), 'M');
        assert_eq!(AminoAcid::Phenylalanine.short_abbreviation(), 'F');
        assert_eq!(AminoAcid::Proline.short_abbreviation(), 'P');
        assert_eq!(AminoAcid::Serine.short_abbreviation(), 'S');
        assert_eq!(AminoAcid::Threonine.short_abbreviation(), 'T');
        assert_eq!(AminoAcid::Tryptophan.short_abbreviation(), 'W');
        assert_eq!(AminoAcid::Tyrosine.short_abbreviation(), 'Y');
        assert_eq!(AminoAcid::Valine.short_abbreviation(), 'V');
        assert_eq!(AminoAcid::Stop.short_abbreviation(), 'X');
    }

    #[test]
    fn test_protein_to_string() {
        let protein = Protein::new(vec![
            AminoAcid::Alanine,
            AminoAcid::Arginine,
            AminoAcid::Asparagine,
            AminoAcid::AsparticAcid,
            AminoAcid::Cysteine,
            AminoAcid::GlutamicAcid,
            AminoAcid::Glutamine,
            AminoAcid::Stop,
        ]);
        assert_eq!(protein.to_string(), "ARNDCEQX");
    }

    #[test]
    #[allow(clippy::vec_init_then_push)]
    fn test_protein_from_transcript() {
        let mut reference = HashMap::new();
        reference.insert(
            "chr1".to_string(),
            "AAAAAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTTAAATTT"
                .bytes()
                .collect(),
        );
        let mut cds = Vec::new();
        cds.push(Cds {
            start: 3,
            end: 9,
            phase: 0,
        });
        cds.push(Cds {
            start: 12,
            end: 18,
            phase: 0,
        });
        let transcript = Transcript {
            feature: "ENSG007".to_string(),
            target: "chr1".to_string(),
            strand: Strand::Forward,
            coding_sequences: cds,
        };
        let protein = Protein::from_transcript(&reference, &transcript);
        assert!(protein.is_ok());
        let expected_protein = Protein {
            sequence: vec![
                AminoAcid::Lysine,
                AminoAcid::Phenylalanine,
                AminoAcid::Isoleucine,
                AminoAcid::Stop,
            ],
        };
        assert_eq!(protein.unwrap(), expected_protein);
    }

    #[test]
    fn test_protein_from_haplotype() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAACCC".to_vec());
        let cds = vec![Cds {
            start: 0,
            end: 5,
            phase: 0,
        }];

        let transcript = Transcript {
            feature: "ENST_TEST".to_string(),
            target: "chr1".to_string(),
            strand: Strand::Forward,
            coding_sequences: cds,
        };

        let haplotype = vec![Node {
            node_type: NodeType::Variant,
            reference_allele: "A".to_string(),
            alternative_allele: "T".to_string(),
            vaf: HashMap::new(),
            probs: EventProbs(HashMap::new()),
            pos: 2,
            index: 0,
        }];

        let protein = Protein::from_haplotype(&reference, &transcript, &haplotype);
        assert!(protein.is_ok());

        let expected = Protein {
            sequence: vec![AminoAcid::Asparagine, AminoAcid::Proline],
        };

        assert_eq!(protein.unwrap(), expected);
    }

    #[test]
    fn from_haplotype_applies_full_insertion_allele() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAACCC".to_vec());
        let transcript = Transcript {
            feature: "ENST_INS".to_string(),
            target: "chr1".to_string(),
            strand: Strand::Forward,
            coding_sequences: vec![Cds {
                start: 0,
                end: 5,
                phase: 0,
            }],
        };
        let haplotype = vec![Node {
            node_type: NodeType::Variant,
            reference_allele: "A".to_string(),
            alternative_allele: "ATT".to_string(),
            vaf: HashMap::new(),
            probs: EventProbs(HashMap::new()),
            pos: 2,
            index: 0,
        }];

        let protein = Protein::from_haplotype(&reference, &transcript, &haplotype).unwrap();

        assert_eq!(
            protein,
            Protein {
                sequence: vec![AminoAcid::Lysine, AminoAcid::Phenylalanine],
            }
        );
    }

    #[test]
    fn from_haplotype_applies_reverse_strand_variants_in_genomic_order() {
        let mut reference = HashMap::new();
        reference.insert("chr1".to_string(), b"AAACCC".to_vec());
        let transcript = Transcript {
            feature: "ENST_REV".to_string(),
            target: "chr1".to_string(),
            strand: Strand::Reverse,
            coding_sequences: vec![Cds {
                start: 0,
                end: 5,
                phase: 0,
            }],
        };
        let haplotype = vec![
            Node {
                node_type: NodeType::Variant,
                reference_allele: "C".to_string(),
                alternative_allele: "CG".to_string(),
                vaf: HashMap::new(),
                probs: EventProbs(HashMap::new()),
                pos: 4,
                index: 1,
            },
            Node {
                node_type: NodeType::Variant,
                reference_allele: "A".to_string(),
                alternative_allele: "T".to_string(),
                vaf: HashMap::new(),
                probs: EventProbs(HashMap::new()),
                pos: 1,
                index: 0,
            },
        ];

        let protein = Protein::from_haplotype(&reference, &transcript, &haplotype).unwrap();

        assert_eq!(
            protein,
            Protein {
                sequence: vec![AminoAcid::Alanine, AminoAcid::Valine],
            }
        );
    }

    #[test]
    fn start_lost_compares_first_residue_against_original() {
        let original = Protein::new(vec![
            AminoAcid::Methionine,
            AminoAcid::Leucine,
            AminoAcid::Valine,
        ]);
        assert!(Protein::new(vec![
            AminoAcid::Threonine,
            AminoAcid::Leucine,
            AminoAcid::Valine
        ])
        .start_lost(&original));
        assert!(!Protein::new(vec![
            AminoAcid::Methionine,
            AminoAcid::Tryptophan,
            AminoAcid::Valine
        ])
        .start_lost(&original));
        let without_annotated_start = Protein::new(vec![AminoAcid::Leucine, AminoAcid::Valine]);
        assert!(!Protein::new(vec![AminoAcid::Threonine, AminoAcid::Valine])
            .start_lost(&without_annotated_start));
    }

    #[test]
    fn apply_start_lost_restarts_at_next_methionine() {
        let mut protein = Protein::new(vec![
            AminoAcid::Threonine,
            AminoAcid::Leucine,
            AminoAcid::Methionine,
            AminoAcid::Valine,
        ]);
        protein.apply_start_lost();
        assert_eq!(
            protein.sequence,
            vec![AminoAcid::Methionine, AminoAcid::Valine]
        );

        let mut without_downstream_start = Protein::new(vec![
            AminoAcid::Threonine,
            AminoAcid::Leucine,
            AminoAcid::Valine,
        ]);
        without_downstream_start.apply_start_lost();
        assert!(without_downstream_start.sequence.is_empty());
    }

    #[test]
    fn test_amino_acid_from_u8() {
        let cases = vec![
            (b'A', AminoAcid::Alanine),
            (b'R', AminoAcid::Arginine),
            (b'N', AminoAcid::Asparagine),
            (b'D', AminoAcid::AsparticAcid),
            (b'E', AminoAcid::GlutamicAcid),
            (b'Q', AminoAcid::Glutamine),
            (b'K', AminoAcid::Lysine),
            (b'H', AminoAcid::Histidine),
            (b'Y', AminoAcid::Tyrosine),
            (b'W', AminoAcid::Tryptophan),
            (b'F', AminoAcid::Phenylalanine),
            (b'S', AminoAcid::Serine),
            (b'T', AminoAcid::Threonine),
            (b'P', AminoAcid::Proline),
            (b'G', AminoAcid::Glycine),
            (b'C', AminoAcid::Cysteine),
            (b'M', AminoAcid::Methionine),
            (b'V', AminoAcid::Valine),
            (b'I', AminoAcid::Isoleucine),
            (b'L', AminoAcid::Leucine),
            (b'X', AminoAcid::Stop),
        ];

        for (byte, expected) in cases {
            assert_eq!(AminoAcid::from(byte), expected);
        }
    }

    #[test]
    #[should_panic]
    fn test_amino_acid_from_invalid_input() {
        let _aa = AminoAcid::from(b'Z');
    }
}
