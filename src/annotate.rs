use crate::annotation::Annotation;
use crate::graph::duck::read_annotation_rows;
use crate::graph::score::HaplotypeFrequency;
use anyhow::{anyhow, bail, Result};
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{Format, Header, Read, Reader, Record, Writer};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::path::Path;

const HAPLODON_INFO: &[u8] = b"##INFO=<ID=HAPLODON,Number=.,Type=String,Description=\"haplodon haplotype-aware variant effect. One entry per transcript and haplotype; scores are haplotype-level and shared across a haplotype's records, joined via haplotype_id (per-sample frequencies in FORMAT/HF). Reserved characters are percent-encoded per VCF 4.3. Format: allele|consequence|feature|hgvsc|hgvsg|haplotype_id|sequence_impact|general_impact|structure_impact|splice_impact|acmg\">";

const HAPLODON_MAX_INFO: &[u8] = b"##INFO=<ID=HAPLODON_MAX,Number=A,Type=Float,Description=\"Maximum haplodon haplotype score for this ALT across all transcripts and haplotypes (0 to 1).\">";

const HF_FORMAT: &[u8] = b"##FORMAT=<ID=HF,Number=.,Type=String,Description=\"Per-sample haplotype frequencies as comma-separated haplotype_id|frequency tokens; haplotype_id matches INFO/HAPLODON.\">";

const VERSION_LINE: &str = concat!("##haplodonVersion=", env!("CARGO_PKG_VERSION"));

/// A variant node of a haplotype, recovered from a stored `hgvsg` token.
#[derive(Debug, PartialEq)]
struct Variant {
    /// Zero-based genomic position, matching `Record::pos`.
    pos: i64,
    reference_allele: String,
    alternative_allele: String,
}

/// A scored `(transcript, haplotype)` ready to write into its records.
struct Haplotype {
    feature: String,
    consequence: String,
    hgvsc: String,
    hgvsg: String,
    id: String,
    score: f64,
    annotation: Annotation,
    frequencies: HaplotypeFrequency,
}

/// Identifies a variant record across the calls file and the scores.
type VariantKey = (String, i64, String, String);

/// Writes the haplotype scores from `scores` back into `calls` as VCF/BCF
/// annotations, using the same variant records that `build` consumed.
pub(crate) fn annotate(calls: &Path, scores: &Path, output: &Path) -> Result<()> {
    let mut haplotypes: Vec<Haplotype> = Vec::new();
    let mut membership: HashMap<VariantKey, Vec<usize>> = HashMap::new();
    for row in read_annotation_rows(scores)? {
        let (target, feature) = row.transcript.split_once(':').ok_or_else(|| {
            anyhow!(
                "Malformed transcript identifier in scores: {}",
                row.transcript
            )
        })?;
        // GFF3 gives CDS features an `ID` like `CDS:ENSP...`; report the bare protein id.
        let feature = feature.strip_prefix("CDS:").unwrap_or(feature);
        let tokens = hgvsg_tokens(&row.hgvsg)?;
        if tokens.is_empty() {
            continue;
        }
        let id = haplotype_id(&row.transcript, &tokens);
        let index = haplotypes.len();
        for token in &tokens {
            let variant = parse_variant(token)?;
            membership
                .entry((
                    target.to_string(),
                    variant.pos,
                    variant.reference_allele,
                    variant.alternative_allele,
                ))
                .or_default()
                .push(index);
        }
        haplotypes.push(Haplotype {
            feature: feature.to_string(),
            consequence: row.consequence,
            hgvsc: encode(&row.hgvsc),
            hgvsg: encode(&row.hgvsg),
            id,
            score: row.score,
            annotation: row.annotation,
            frequencies: row.frequencies,
        });
    }

    let mut reader = Reader::from_path(calls)?;
    let input_header = reader.header().clone();
    for tag in ["HAPLODON", "HAPLODON_MAX", "HF"] {
        if input_header.info_type(tag.as_bytes()).is_ok()
            || input_header.format_type(tag.as_bytes()).is_ok()
        {
            bail!(
                "calls file {} already declares {tag}; annotate expects the original calls file from build",
                calls.display()
            );
        }
    }
    let samples: Vec<String> = input_header
        .samples()
        .iter()
        .map(|sample| String::from_utf8_lossy(sample).into_owned())
        .collect();

    let mut header = Header::from_template(&input_header);
    header.push_record(HAPLODON_INFO);
    header.push_record(HAPLODON_MAX_INFO);
    header.push_record(HF_FORMAT);
    header.push_record(VERSION_LINE.as_bytes());
    header.push_record(
        format!(
            "##haplodonAnnotateCommand=annotate --calls {} --scores {} --output {}",
            calls.display(),
            scores.display(),
            output.display()
        )
        .as_bytes(),
    );

    let (uncompressed, format) = output_format(output);
    let mut writer = Writer::from_path(output, &header, uncompressed, format)?;

    let mut matched = 0usize;
    for result in reader.records() {
        let mut record = result?;
        let target = match record.rid() {
            Some(rid) => Some(String::from_utf8(input_header.rid2name(rid)?.to_vec())?),
            None => None,
        };
        writer.translate(&mut record);
        if let Some(target) = target {
            if annotate_record(&mut record, &target, &haplotypes, &membership, &samples)? {
                matched += 1;
            }
        }
        writer.write(&record)?;
    }

    if !haplotypes.is_empty() && matched == 0 {
        bail!(
            "no record in {} matched any of the {} scored haplotypes; is this the calls file that was used by build?",
            calls.display(),
            haplotypes.len()
        );
    }
    Ok(())
}

/// Adds the haplodon annotations to a single record, leaving records that are not
/// part of any scored haplotype untouched.
fn annotate_record(
    record: &mut Record,
    target: &str,
    haplotypes: &[Haplotype],
    membership: &HashMap<VariantKey, Vec<usize>>,
    samples: &[String],
) -> Result<bool> {
    let position = record.pos();
    let (reference, alternatives) = {
        let alleles = record.alleles();
        let reference = String::from_utf8_lossy(alleles[0]).into_owned();
        let alternatives: Vec<String> = alleles[1..]
            .iter()
            .map(|allele| String::from_utf8_lossy(allele).into_owned())
            .collect();
        (reference, alternatives)
    };

    let mut entries: Vec<String> = Vec::new();
    let mut maxima: Vec<f32> = Vec::with_capacity(alternatives.len());
    let mut frequencies: Vec<Vec<String>> = vec![Vec::new(); samples.len()];
    let mut annotated = false;

    for alternative in &alternatives {
        let key = (
            target.to_string(),
            position,
            reference.clone(),
            alternative.clone(),
        );
        let mut maximum = f32::missing();
        if let Some(indices) = membership.get(&key) {
            let mut indices = indices.clone();
            indices.sort_by(|&a, &b| {
                haplotypes[b]
                    .score
                    .total_cmp(&haplotypes[a].score)
                    .then_with(|| haplotypes[a].id.cmp(&haplotypes[b].id))
            });
            for &index in &indices {
                let haplotype = &haplotypes[index];
                entries.push(format_entry(haplotype, alternative));
                let score = haplotype.score as f32;
                if maximum.is_missing() || score > maximum {
                    maximum = score;
                }
                for (sample, tokens) in samples.iter().zip(frequencies.iter_mut()) {
                    // Only carriers (a non-zero haplotype frequency) are listed per sample.
                    if let Some(&frequency) = haplotype.frequencies.get(sample) {
                        if frequency > 0.0 {
                            tokens.push(format!(
                                "{}|{}",
                                haplotype.id,
                                format_float(frequency as f64)
                            ));
                        }
                    }
                }
                annotated = true;
            }
        }
        maxima.push(maximum);
    }

    if !annotated {
        return Ok(false);
    }

    let info = entries.join(",");
    record.push_info_string(b"HAPLODON", &[info.as_bytes()])?;
    record.push_info_float(b"HAPLODON_MAX", &maxima)?;
    if !samples.is_empty() {
        let values: Vec<Vec<u8>> = frequencies
            .into_iter()
            .map(|mut tokens| {
                if tokens.is_empty() {
                    b".".to_vec()
                } else {
                    tokens.sort();
                    tokens.join(",").into_bytes()
                }
            })
            .collect();
        record.push_format_string(b"HF", &values)?;
    }
    Ok(true)
}

/// Renders one `HAPLODON` entry for the given ALT.
fn format_entry(haplotype: &Haplotype, allele: &str) -> String {
    let annotation = &haplotype.annotation;
    [
        encode(allele),
        encode(&haplotype.consequence),
        encode(&haplotype.feature),
        haplotype.hgvsc.clone(),
        haplotype.hgvsg.clone(),
        haplotype.id.clone(),
        format_float(haplotype.score),
        format_optional(annotation.revel_score),
        format_optional(annotation.alphamissense_score),
        format_optional(annotation.spliceai_score),
        format_optional(annotation.acmg_score),
    ]
    .join("|")
}

/// Derives a stable, order-independent identifier for a haplotype from its
/// transcript and its sorted set of `hgvsg` tokens.
fn haplotype_id(transcript: &str, tokens: &[&str]) -> String {
    let mut tokens = tokens.to_vec();
    tokens.sort_unstable();
    let mut hasher = Sha256::new();
    hasher.update(transcript.as_bytes());
    hasher.update(b"\x00");
    hasher.update(tokens.join(";").as_bytes());
    hasher
        .finalize()
        .iter()
        .take(6)
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

/// Splits a composite `hgvsg` string such as `g.[467C>T;1211G>A]` into its tokens.
fn hgvsg_tokens(hgvsg: &str) -> Result<Vec<&str>> {
    let inner = hgvsg
        .strip_prefix("g.[")
        .and_then(|rest| rest.strip_suffix(']'))
        .ok_or_else(|| anyhow!("Unexpected hgvsg notation: {hgvsg}"))?;
    if inner.is_empty() {
        Ok(Vec::new())
    } else {
        Ok(inner.split(';').collect())
    }
}

/// Parses a single `hgvsg` token such as `467C>T` (a one-based `ref>alt` change).
fn parse_variant(token: &str) -> Result<Variant> {
    let (change, alternative) = token
        .split_once('>')
        .ok_or_else(|| anyhow!("Unexpected hgvsg token: {token}"))?;
    let split = change
        .find(|character: char| !character.is_ascii_digit())
        .ok_or_else(|| anyhow!("Unexpected hgvsg token: {token}"))?;
    let (position, reference) = change.split_at(split);
    Ok(Variant {
        pos: position.parse::<i64>()? - 1,
        reference_allele: reference.to_string(),
        alternative_allele: alternative.to_string(),
    })
}

/// Percent-encodes the characters that VCF 4.3 reserves within a subfield.
fn encode(value: &str) -> String {
    let mut encoded = String::with_capacity(value.len());
    for &byte in value.as_bytes() {
        if needs_encoding(byte) {
            encoded.push_str(&format!("%{byte:02X}"));
        } else {
            encoded.push(byte as char);
        }
    }
    encoded
}

fn needs_encoding(byte: u8) -> bool {
    matches!(byte, b'%' | b';' | b',' | b'=' | b'|' | b':') || byte <= b' ' || byte > b'~'
}

fn format_optional(value: Option<f64>) -> String {
    value.map(format_float).unwrap_or_default()
}

/// Formats a score for a VCF field, dropping trailing zeros; non-finite values become empty.
fn format_float(value: f64) -> String {
    if !value.is_finite() {
        return String::new();
    }
    let formatted = format!("{value:.4}");
    let trimmed = formatted.trim_end_matches('0').trim_end_matches('.');
    trimmed.to_string()
}

/// Chooses the output format from the file extension.
fn output_format(path: &Path) -> (bool, Format) {
    let name = path.to_string_lossy().to_lowercase();
    if name.ends_with(".bcf") {
        (false, Format::Bcf)
    } else if name.ends_with(".vcf.gz") || name.ends_with(".vcf.bgz") {
        (false, Format::Vcf)
    } else {
        (true, Format::Vcf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::Annotation;
    use crate::graph::duck::{create_scores, write_scores};
    use crate::graph::paths::Cds;
    use crate::graph::score::{Consequence, EffectScore};
    use crate::graph::transcript::Transcript;
    use crate::translation::amino_acids::{AminoAcid, Protein};
    use crate::translation::distance::DistanceMetric;
    use bio::bio_types::strand::Strand;
    use std::collections::HashMap;
    use std::io::Write;

    #[test]
    fn encode_escapes_reserved_characters_only() {
        assert_eq!(encode("g.[1C>T;2G>A]"), "g.[1C>T%3B2G>A]");
        assert_eq!(encode("a,b=c|d:e f"), "a%2Cb%3Dc%7Cd%3Ae%20f");
        assert_eq!(encode("100%"), "100%25");
        assert_eq!(encode("missense_variant"), "missense_variant");
    }

    #[test]
    fn hgvsg_tokens_split_and_parse_substitutions_and_indels() {
        assert_eq!(hgvsg_tokens("g.[]").unwrap(), Vec::<&str>::new());
        assert_eq!(
            hgvsg_tokens("g.[467C>T;1211G>A]").unwrap(),
            vec!["467C>T", "1211G>A"]
        );
        assert_eq!(
            parse_variant("467C>T").unwrap(),
            Variant {
                pos: 466,
                reference_allele: "C".to_string(),
                alternative_allele: "T".to_string(),
            }
        );
        assert_eq!(
            parse_variant("5AT>A").unwrap(),
            Variant {
                pos: 4,
                reference_allele: "AT".to_string(),
                alternative_allele: "A".to_string(),
            }
        );
    }

    #[test]
    fn hgvsg_tokens_and_parse_variant_reject_malformed_notation() {
        assert!(hgvsg_tokens("c.[1C>T]").is_err());
        assert!(parse_variant("1CT").is_err());
    }

    #[test]
    fn haplotype_id_is_stable_and_independent_of_token_order() {
        let forward = ["467C>T", "1211G>A"];
        let reverse = ["1211G>A", "467C>T"];
        let id = haplotype_id("chr1:ENSP0", &forward);
        assert_eq!(id.len(), 12);
        assert_eq!(id, haplotype_id("chr1:ENSP0", &reverse));
        assert_ne!(id, haplotype_id("chr1:ENSP1", &forward));
    }

    #[test]
    fn format_float_trims_trailing_zeros_and_hides_non_finite() {
        assert_eq!(format_float(0.42), "0.42");
        assert_eq!(format_float(1.0), "1");
        assert_eq!(format_float(0.0), "0");
        assert_eq!(format_float(f64::NAN), "");
    }

    #[test]
    fn output_format_follows_the_extension() {
        assert!(matches!(
            output_format(Path::new("a.bcf")),
            (false, Format::Bcf)
        ));
        assert!(matches!(
            output_format(Path::new("a.vcf.gz")),
            (false, Format::Vcf)
        ));
        assert!(matches!(
            output_format(Path::new("a.vcf")),
            (true, Format::Vcf)
        ));
    }

    fn write_scores_fixture(path: &Path) {
        create_scores(path).unwrap();
        let transcript = Transcript::new(
            "ENSP00000TEST".to_string(),
            "chr1".to_string(),
            Strand::Forward,
            vec![Cds::new(0, 100, 0)],
        );
        let effect_score = EffectScore {
            original_protein: Protein::new(vec![AminoAcid::Methionine, AminoAcid::Leucine]),
            altered_protein: Protein::new(vec![AminoAcid::Methionine, AminoAcid::Valine]),
            distance_metric: DistanceMetric::Grantham,
            realign: false,
            consequence: Consequence::Missense,
            hgvsc: "c.[100C>T]".to_string(),
            hgvsg: "g.[7675089C>T]".to_string(),
            hgvsg_full: "g.[7675089C>T]".to_string(),
        };
        let frequencies = HashMap::from([("S1".to_string(), 0.42)]);
        let annotation = Annotation {
            revel_score: Some(0.31),
            acmg_score: None,
            spliceai_score: None,
            alphamissense_score: None,
        };
        write_scores(
            path,
            vec![(effect_score, frequencies, vec![], annotation)],
            transcript,
        )
        .unwrap();
    }

    fn write_calls_vcf(path: &Path, records: &str) {
        let mut file = std::fs::File::create(path).unwrap();
        write!(
            file,
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr1>\n\
             ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
             {records}"
        )
        .unwrap();
    }

    #[test]
    fn annotate_writes_haplotype_fields_into_matching_records() {
        let dir = tempfile::tempdir().unwrap();
        let calls = dir.path().join("calls.vcf");
        let scores = dir.path().join("scores.duckdb");
        let output = dir.path().join("annotated.vcf");

        write_calls_vcf(
            &calls,
            "chr1\t7675089\t.\tC\tT\t60\tPASS\tDP=133\tGT\t0/1\t0/1\n\
             chr1\t7680000\t.\tA\tG\t60\tPASS\tDP=90\tGT\t0/1\t0/1\n",
        );
        write_scores_fixture(&scores);
        annotate(&calls, &scores, &output).unwrap();

        let mut reader = Reader::from_path(&output).unwrap();
        let mut records = reader.records();

        let annotated = records.next().unwrap().unwrap();
        let haplodon = {
            let info = annotated.info(b"HAPLODON");
            let buffer = info.string().unwrap().unwrap();
            String::from_utf8(buffer[0].to_vec()).unwrap()
        };
        assert!(haplodon.starts_with("T|missense_variant|ENSP00000TEST|c.[100C>T]|g.[7675089C>T]|"));
        assert!(haplodon.contains("|0.31|||"));
        let maximum = annotated.info(b"HAPLODON_MAX").float().unwrap().unwrap()[0];
        let score: f32 = haplodon.split('|').nth(6).unwrap().parse().unwrap();
        assert!(maximum > 0.0 && maximum < 1.0);
        assert!((maximum - score).abs() < 1e-4);
        let frequencies = {
            let format = annotated.format(b"HF");
            let buffer = format.string().unwrap();
            (
                String::from_utf8(buffer[0].to_vec()).unwrap(),
                String::from_utf8(buffer[1].to_vec()).unwrap(),
            )
        };
        // S1 carries the haplotype; S2 has no frequency, so it is dropped to a missing value.
        assert!(frequencies.0.ends_with("|0.42"));
        assert_eq!(frequencies.1, ".");

        let untouched = records.next().unwrap().unwrap();
        assert!(untouched.info(b"HAPLODON").string().unwrap().is_none());
    }

    #[test]
    fn annotate_orders_entries_by_descending_score() {
        let dir = tempfile::tempdir().unwrap();
        let calls = dir.path().join("calls.vcf");
        let scores = dir.path().join("scores.duckdb");
        let output = dir.path().join("out.vcf");

        write_calls_vcf(
            &calls,
            "chr1\t7675089\t.\tC\tT\t60\tPASS\tDP=99\tGT\t0/1\t0/1\n",
        );

        // Two transcripts carry the same variant with different divergence, so the
        // record lists both haplotypes, most divergent first.
        create_scores(&scores).unwrap();
        for (feature, altered, consequence) in [
            ("ENSP00000ZERO", AminoAcid::Leucine, Consequence::Synonymous),
            ("ENSP00000HIGH", AminoAcid::Valine, Consequence::Missense),
        ] {
            let effect_score = EffectScore {
                original_protein: Protein::new(vec![AminoAcid::Methionine, AminoAcid::Leucine]),
                altered_protein: Protein::new(vec![AminoAcid::Methionine, altered]),
                distance_metric: DistanceMetric::Grantham,
                realign: false,
                consequence,
                hgvsc: "c.[100C>T]".to_string(),
                hgvsg: "g.[7675089C>T]".to_string(),
                hgvsg_full: "g.[7675089C>T]".to_string(),
            };
            let annotation = Annotation {
                revel_score: None,
                acmg_score: None,
                spliceai_score: None,
                alphamissense_score: None,
            };
            let transcript = Transcript::new(
                feature.to_string(),
                "chr1".to_string(),
                Strand::Forward,
                vec![Cds::new(0, 100, 0)],
            );
            write_scores(
                &scores,
                vec![(
                    effect_score,
                    HashMap::from([("S1".to_string(), 0.5)]),
                    vec![],
                    annotation,
                )],
                transcript,
            )
            .unwrap();
        }

        annotate(&calls, &scores, &output).unwrap();
        let mut reader = Reader::from_path(&output).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let info = record.info(b"HAPLODON");
        let buffer = info.string().unwrap().unwrap();
        let features: Vec<String> = buffer
            .iter()
            .map(|entry| {
                String::from_utf8_lossy(entry)
                    .split('|')
                    .nth(2)
                    .unwrap()
                    .to_string()
            })
            .collect();
        assert_eq!(features, ["ENSP00000HIGH", "ENSP00000ZERO"]);
    }

    #[test]
    fn annotate_fails_when_no_record_matches_a_scored_haplotype() {
        let dir = tempfile::tempdir().unwrap();
        let calls = dir.path().join("calls.vcf");
        let scores = dir.path().join("scores.duckdb");
        let output = dir.path().join("out.vcf");

        write_calls_vcf(
            &calls,
            "chr1\t4242\t.\tA\tG\t60\tPASS\tDP=10\tGT\t0/1\t0/1\n",
        );
        write_scores_fixture(&scores);

        let error = annotate(&calls, &scores, &output).unwrap_err();
        assert!(error.to_string().contains("matched any"));
    }

    #[test]
    fn annotate_fails_when_calls_file_is_already_annotated() {
        let dir = tempfile::tempdir().unwrap();
        let calls = dir.path().join("calls.vcf");
        let scores = dir.path().join("scores.duckdb");
        let output = dir.path().join("out.vcf");

        let mut file = std::fs::File::create(&calls).unwrap();
        write!(
            file,
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr1>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             ##FORMAT=<ID=HF,Number=.,Type=String,Description=\"stale\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
             chr1\t7675089\t.\tC\tT\t60\tPASS\t.\tGT:HF\t0/1:stale\n"
        )
        .unwrap();
        drop(file);
        write_scores_fixture(&scores);

        let error = annotate(&calls, &scores, &output).unwrap_err();
        assert!(error.to_string().contains("already declares HF"));
    }
}
