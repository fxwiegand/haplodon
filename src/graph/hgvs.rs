use crate::graph::transcript::Transcript;
use bio::alphabets::dna;
use bio::bio_types::strand::Strand;

pub(crate) fn hgvsc(
    transcript: &Transcript,
    position: u64,
    reference_allele: &str,
    alternative_allele: &str,
) -> Option<String> {
    let cdna_pos = transcript.position_in_transcript(position as usize).ok()? + 1; // 1-based, strand-aware

    let (ref_allele, alt_allele) = if transcript.strand == Strand::Reverse {
        (
            String::from_utf8(dna::revcomp(reference_allele.as_bytes())).unwrap(),
            String::from_utf8(dna::revcomp(alternative_allele.as_bytes())).unwrap(),
        )
    } else {
        (reference_allele.to_string(), alternative_allele.to_string())
    };

    let notation = match (ref_allele.is_empty(), alt_allele.is_empty()) {
        (true, false) => format!("{cdna_pos}_{}ins{alt_allele}", cdna_pos + 1),
        (false, true) => match ref_allele.len() {
            1 => format!("{cdna_pos}del"),
            n => format!("{cdna_pos}_{}del", cdna_pos + n - 1),
        },
        (false, false) => match (ref_allele.len(), alt_allele.len()) {
            (1, 1) => format!("{cdna_pos}{ref_allele}>{alt_allele}"),
            (1, _) => format!("{cdna_pos}delins{alt_allele}"),
            (n, _) => format!("{cdna_pos}_{}delins{alt_allele}", cdna_pos + n - 1),
        },
        (true, true) => return None,
    };

    Some(notation)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{graph::paths::Cds, graph::transcript::Transcript, transcripts};
    use bio::bio_types::strand::Strand;
    use std::path::PathBuf;

    fn transcript(cdss: Vec<Cds>, strand: Strand) -> Transcript {
        Transcript {
            feature: "ENST1".into(),
            target: "chr1".into(),
            strand,
            coding_sequences: cdss,
        }
    }

    #[test]
    fn snv_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                114,
                "C",
                "G"
            ),
            Some("15C>G".into())
        );
    }

    #[test]
    fn snv_first_base_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                100,
                "A",
                "T"
            ),
            Some("1A>T".into())
        );
    }

    #[test]
    fn del_single_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                105,
                "A",
                ""
            ),
            Some("6del".into())
        );
    }

    #[test]
    fn del_multi_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                105,
                "ATG",
                ""
            ),
            Some("6_8del".into())
        );
    }

    #[test]
    fn ins_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                105,
                "",
                "ATTCAC"
            ),
            Some("6_7insATTCAC".into())
        );
    }

    #[test]
    fn multi_exon_forward() {
        let t = transcript(
            vec![Cds::new(100, 109, 0), Cds::new(200, 299, 0)],
            Strand::Forward,
        );
        assert_eq!(hgvsc(&t, 200, "C", "T"), Some("11C>T".into()));
        assert_eq!(hgvsc(&t, 204, "G", "A"), Some("15G>A".into()));
    }

    #[test]
    fn snv_reverse_allele_flip() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Reverse),
                114,
                "C",
                "G"
            ),
            Some("87G>C".into())
        );
    }

    #[test]
    fn snv_last_base_reverse() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Reverse),
                200,
                "A",
                "T"
            ),
            Some("1T>A".into())
        );
    }

    #[test]
    fn ins_reverse() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Reverse),
                195,
                "",
                "ATG"
            ),
            Some("6_7insCAT".into())
        );
    }

    #[test]
    fn outside_cds() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                50,
                "A",
                "T"
            ),
            None
        );
    }

    #[test]
    fn empty_alleles() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                105,
                "",
                ""
            ),
            None
        );
    }

    #[test]
    fn delins_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                105,
                "ATG",
                "TC"
            ),
            Some("6_8delinsTC".into())
        );
    }

    #[test]
    fn delins_single_ref_forward() {
        assert_eq!(
            hgvsc(
                &transcript(vec![Cds::new(100, 200, 0)], Strand::Forward),
                105,
                "A",
                "TT"
            ),
            Some("6delinsTT".into())
        );
    }

    #[test]
    fn snv_reverse_multi_exon() {
        let gff_file = PathBuf::from("tests/resources/ENSP00000355304.gff3");
        let graph = PathBuf::from("tests/resources/graph.duckdb");
        let transcripts = transcripts(&gff_file, &graph).unwrap();
        let t = &transcripts[0];
        assert_eq!(hgvsc(t, 102229355, "T", "C"), Some("1193A>G".into()));
        assert_eq!(hgvsc(t, 102263548, "G", "A"), Some("280C>T".into()));
    }
}
