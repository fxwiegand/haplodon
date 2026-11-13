//! Labelled missense and synonymous SNVs from ClinVar, selected like in Pejaver et al. 2022.

use anyhow::{bail, Context, Result};
use genebears::Variant;
use itertools::Itertools;
use rust_htslib::bcf::header::HeaderRecord;
use rust_htslib::bcf::{self, Read};
use std::collections::HashSet;
use std::path::Path;

const PATHOGENIC: [&str; 3] = [
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
];
const BENIGN: [&str; 3] = ["Benign", "Likely_benign", "Benign/Likely_benign"];
/// Review statuses with at least one star.
const REVIEWED: [&str; 4] = [
    "criteria_provided,_single_submitter",
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "reviewed_by_expert_panel",
    "practice_guideline",
];
/// Consequences on any transcript that can make a variant pathogenic other than through its group.
const OTHER_MECHANISMS: [&str; 7] = [
    "nonsense",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "initiator_codon_variant",
    "start_lost",
    "stop_lost",
];

#[derive(Debug)]
pub(crate) struct Labelled {
    pub(crate) variant: Variant,
    /// Missense on some transcript, otherwise synonymous.
    pub(crate) missense: bool,
    pub(crate) pathogenic: bool,
}

/// The release from the header and the SNVs with 1+ stars, labelled pathogenic or benign, in genes
/// with at least one pathogenic record of any kind. Variants listed twice with different labels are
/// dropped.
pub(crate) fn read(path: &Path) -> Result<(String, Vec<Labelled>)> {
    let mut reader = bcf::Reader::from_path(path)?;
    let records = reader.header().header_records();
    let header = |tag: &str| {
        records.iter().find_map(|record| match record {
            HeaderRecord::Generic { key, value } if key == tag => Some(value.as_str()),
            _ => None,
        })
    };
    // The scores are looked up for hg38. The curves work for GRCh37 data as well.
    if header("reference") != Some("GRCh38") {
        bail!("The ClinVar VCF must be the one for GRCh38");
    }
    let release = header("fileDate")
        .context("ClinVar VCF has no fileDate")?
        .replace('-', "");

    let mut disease_genes = HashSet::new();
    let mut candidates = Vec::new();
    for record in reader.records() {
        let record = record?;
        let significance = info(&record, b"CLNSIG")?;
        let pathogenic = PATHOGENIC.contains(&significance.as_str());
        let genes = info(&record, b"GENEINFO")?
            .split('|')
            .filter_map(|gene| gene.split(':').next())
            .filter(|gene| !gene.is_empty())
            .map(String::from)
            .collect_vec();
        if pathogenic {
            disease_genes.extend(genes.iter().cloned());
        }
        if !pathogenic && !BENIGN.contains(&significance.as_str())
            || !REVIEWED.contains(&info(&record, b"CLNREVSTAT")?.as_str())
        {
            continue;
        }
        let alleles = record.alleles();
        let [reference, alternative] = alleles.as_slice() else {
            continue;
        };
        let snv = |allele: &[u8]| matches!(allele, [b'A' | b'C' | b'G' | b'T']);
        if !snv(reference) || !snv(alternative) {
            continue;
        }
        let consequences = info(&record, b"MC")?;
        let terms = consequences
            .split(',')
            .filter_map(|consequence| consequence.split_once('|'))
            .map(|(_, term)| term)
            .collect_vec();
        let missense = terms.contains(&"missense_variant");
        if terms.iter().any(|term| OTHER_MECHANISMS.contains(term))
            || !missense && !terms.contains(&"synonymous_variant")
        {
            continue;
        }
        let rid = record.rid().context("ClinVar record without chromosome")?;
        let variant = Variant::new(
            std::str::from_utf8(record.header().rid2name(rid)?)?,
            record.pos() as u64 + 1,
            std::str::from_utf8(reference)?,
            std::str::from_utf8(alternative)?,
        );
        candidates.push((
            genes,
            Labelled {
                variant,
                missense,
                pathogenic,
            },
        ));
    }

    let variants = candidates
        .into_iter()
        .filter(|(genes, _)| genes.iter().any(|gene| disease_genes.contains(gene)))
        .map(|(_, labelled)| labelled)
        .into_group_map_by(|labelled| labelled.variant.clone())
        .into_values()
        .filter(|duplicates| duplicates.iter().map(|l| l.pathogenic).all_equal())
        .filter_map(|duplicates| duplicates.into_iter().next())
        .collect();
    Ok((release, variants))
}

/// The INFO value, joined again at the commas rust-htslib splits at.
fn info(record: &bcf::Record, tag: &[u8]) -> Result<String> {
    Ok(record
        .info(tag)
        .string()?
        .map(|values| values.iter().map(|v| String::from_utf8_lossy(v)).join(","))
        .unwrap_or_default())
}
