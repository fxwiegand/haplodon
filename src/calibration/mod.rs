//! Turns raw variant scores into probabilities of pathogenicity, following Pejaver et al. 2022
//! (https://doi.org/10.1016/j.ajhg.2022.10.013): ClinVar variants binned by score give a
//! likelihood ratio per bin, which Bayes' rule turns into a probability with the prior of the
//! variant group.

mod clinvar;

use crate::annotation::{gnomad_af, Scores, VariantAnnotation, HUB_DATABASES};
use crate::graph::transcript::read_transcripts;
use anyhow::{bail, Context, Result};
use clap::ValueEnum;
use itertools::Itertools;
use log::info;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;

/// Variants at least this frequent in gnomAD get no probability, since the curves are for rare
/// ones.
const MAX_AF: f64 = 0.01;

/// Minimum variants per bin before merging. Only the small edge group needs it, since the merging
/// and `MIN_LABEL` already make the other bins larger.
const MIN_BIN: f64 = 50.0;

/// Minimum of each label per bin, so that the bins at the ends stay stable.
const MIN_LABEL: f64 = 10.0;

const BUILT_IN: &str = include_str!("../../resources/calibration.tsv");

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, ValueEnum)]
#[serde(rename_all = "lowercase")]
#[value(rename_all = "lower")]
pub(crate) enum Score {
    Revel,
    AlphaMissense,
    SpliceAi,
    Cadd,
}

impl Score {
    /// The groups the score is learned on.
    fn groups(self) -> &'static [Group] {
        match self {
            Score::Revel | Score::AlphaMissense => &[Group::Missense],
            Score::SpliceAi => &[Group::SynonymousEdge, Group::SynonymousOther],
            Score::Cadd => &[
                Group::Missense,
                Group::SynonymousEdge,
                Group::SynonymousOther,
            ],
        }
    }

    fn value(self, annotation: &VariantAnnotation) -> Option<f64> {
        let value = match self {
            Score::Revel => annotation.revel_score,
            Score::AlphaMissense => annotation.alphamissense_score,
            Score::SpliceAi => annotation.spliceai_max_score,
            Score::Cadd => annotation.cadd_score,
        };
        value.filter(|value| value.is_finite())
    }

    /// The GeneBe Hub database of the score.
    fn database(self) -> &'static str {
        match self {
            Score::Revel => HUB_DATABASES[0].0,
            Score::AlphaMissense => HUB_DATABASES[1].0,
            Score::SpliceAi => HUB_DATABASES[2].0,
            Score::Cadd => HUB_DATABASES[5].0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum Group {
    Missense,
    /// Synonymous at the last base of an exon, the exonic base of a splice donor.
    SynonymousEdge,
    SynonymousOther,
}

impl Group {
    pub(crate) fn new(missense: bool, edge: bool) -> Group {
        match (missense, edge) {
            (true, _) => Group::Missense,
            (false, true) => Group::SynonymousEdge,
            (false, false) => Group::SynonymousOther,
        }
    }

    /// The share of pathogenic variants among the rare ones of the group. Missense is Pejaver's
    /// estimate, the synonymous ones come from gnomAD v4.1 and ClinVar, scaled to match it
    /// (koesterlab/haplodon-bench, workflow/scripts/calibration/priors.py).
    fn prior(self) -> f64 {
        match self {
            Group::Missense => 0.0441,
            Group::SynonymousEdge => 0.149,
            Group::SynonymousOther => 0.000736,
        }
    }
}

/// A bin of a curve, from `lower` up to the next bin.
#[derive(Debug, Serialize, Deserialize)]
struct Bin {
    score: Score,
    group: Group,
    /// The ClinVar release the bin was learned on.
    clinvar: String,
    /// The GeneBe Hub database or the API the scores came from.
    source: String,
    prior: f64,
    lower: f64,
    /// Labelled variants, including a pseudo-variant split by the shares of the whole group.
    pathogenic: f64,
    benign: f64,
    probability: f64,
}

/// The curves by score and group, each sorted by score.
#[derive(Debug)]
pub(crate) struct Calibration(BTreeMap<(Score, Group), Vec<Bin>>);

impl Calibration {
    /// The table at the path, or the one built into haplodon.
    pub(crate) fn load(path: Option<&Path>) -> Result<Self> {
        let table = path
            .map(|path| {
                fs::read_to_string(path).with_context(|| format!("Cannot read {}", path.display()))
            })
            .transpose()?;
        let calibration = Self::parse(table.as_deref().unwrap_or(BUILT_IN))?;
        for score in Score::value_variants() {
            for group in score.groups() {
                if !calibration.0.contains_key(&(*score, *group)) {
                    bail!("The calibration table has no curve for {score:?} on {group:?} variants");
                }
            }
        }
        Ok(calibration)
    }

    fn parse(table: &str) -> Result<Self> {
        let mut curves: BTreeMap<_, Vec<Bin>> = BTreeMap::new();
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(table.as_bytes());
        for bin in reader.deserialize() {
            let bin: Bin = bin?;
            curves.entry((bin.score, bin.group)).or_default().push(bin);
        }
        for bins in curves.values_mut() {
            bins.sort_by(|a, b| a.lower.total_cmp(&b.lower));
        }
        Ok(Calibration(curves))
    }

    fn write(&self, path: &Path) -> Result<()> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .with_context(|| format!("Cannot write {}", path.display()))?;
        for bin in self.0.values().flatten() {
            writer.serialize(bin)?;
        }
        writer.flush()?;
        Ok(())
    }

    /// The probability of the bin the value falls into. Values below the first bin take the first.
    fn probability(&self, score: Score, group: Group, value: f64) -> Option<f64> {
        let bins = self.0.get(&(score, group))?;
        let above = bins.partition_point(|bin| bin.lower <= value);
        bins.get(above.saturating_sub(1)).map(|bin| bin.probability)
    }

    /// The scores of a variant of the group as probabilities, without the ones that have no curve
    /// for the group. Common variants get none.
    pub(crate) fn apply(
        &self,
        annotation: &VariantAnnotation,
        group: Option<Group>,
        edge: bool,
    ) -> VariantAnnotation {
        let group = group.filter(|_| rare(annotation));
        let probability = |score: Score| {
            let group = group?;
            // SpliceAI takes the synonymous curve of the same position for missense variants as
            // well, since splicing doesn't depend on the amino acid.
            let curve = match score {
                Score::SpliceAi => Group::new(false, edge),
                _ => group,
            };
            self.probability(score, curve, score.value(annotation)?)
        };
        VariantAnnotation {
            revel_score: probability(Score::Revel),
            alphamissense_score: probability(Score::AlphaMissense),
            spliceai_max_score: probability(Score::SpliceAi),
            cadd_score: probability(Score::Cadd),
            ..annotation.clone()
        }
    }

    /// The ClinVar releases the curves were learned on.
    pub(crate) fn sources(&self) -> String {
        self.0
            .values()
            .flatten()
            .map(|bin| format!("ClinVar {}", bin.clinvar))
            .unique()
            .join(",")
    }
}

/// Learns the curves of the scores on the ClinVar VCF and replaces their rows in the table.
pub(crate) fn calibrate(
    clinvar: &Path,
    features: &Path,
    scores: &[Score],
    source: &Scores,
    table: &Path,
) -> Result<()> {
    if scores.contains(&Score::Cadd) && !matches!(source, Scores::Hub(_)) {
        bail!(
            "CADD scores are only in the GeneBe Hub databases, download them with `haplodon cache`"
        );
    }
    // A new table starts from the built-in one, so that it has all curves. Writing it right away
    // fails before the long annotation if the path is wrong.
    let mut calibration = Calibration::load(Some(table).filter(|table| table.exists()))?;
    calibration.write(table)?;
    calibration
        .0
        .retain(|(score, _), _| !scores.contains(score));
    let chromosome = |name: &str| name.trim_start_matches("chr").to_string();
    // 1-based like the ClinVar variants.
    let edges: HashSet<(String, u64)> = read_transcripts(features)?
        .iter()
        .flat_map(|transcript| {
            let chromosome = chromosome(&transcript.target);
            transcript
                .exon_edges()
                .into_iter()
                .map(move |pos| (chromosome.clone(), pos + 1))
        })
        .collect();
    let (release, labelled) = clinvar::read(clinvar)?;
    let (labels, variants): (Vec<_>, Vec<_>) = labelled
        .into_iter()
        .map(|labelled| {
            let variant = labelled.variant;
            let edge = edges.contains(&(chromosome(&variant.chr), variant.pos));
            (
                (Group::new(labelled.missense, edge), labelled.pathogenic),
                variant,
            )
        })
        .unzip();
    if !labels
        .iter()
        .any(|(group, _)| *group == Group::SynonymousEdge)
    {
        bail!(
            "No ClinVar variant is at an exon edge of the features, do the chromosome names match?"
        );
    }
    info!("Annotating {} ClinVar SNVs", variants.len());
    let annotations = source
        .annotate(&variants)
        .context("GeneBe failed to annotate the ClinVar SNVs")?;
    let labelled = labels
        .into_iter()
        .zip(annotations)
        .filter_map(|(label, annotation)| Some((label, annotation.filter(rare)?)))
        .collect_vec();

    for &score in scores {
        let source = match source {
            Scores::Hub(_) => score.database().to_string(),
            api => api.sources(),
        };
        for &group in score.groups() {
            let values = labelled
                .iter()
                .filter(|((variant_group, _), _)| *variant_group == group)
                .filter_map(|((_, pathogenic), annotation)| {
                    Some((score.value(annotation)?, *pathogenic))
                })
                .collect_vec();
            info!(
                "Calibrating {score:?} on {} {group:?} variants",
                values.len()
            );
            let bins = curve(values, group)
                .with_context(|| format!("Cannot calibrate {score:?} on {group:?} variants"))?
                .into_iter()
                .map(|(counts, probability)| Bin {
                    score,
                    group,
                    clinvar: release.clone(),
                    source: source.clone(),
                    prior: group.prior(),
                    lower: counts.lower,
                    pathogenic: rounded(counts.pathogenic),
                    benign: rounded(counts.benign),
                    probability: rounded(probability),
                })
                .collect();
            calibration.0.insert((score, group), bins);
        }
    }
    calibration.write(table)
}

/// Below `MAX_AF` in gnomAD or not in it.
fn rare(annotation: &VariantAnnotation) -> bool {
    gnomad_af(annotation).is_none_or(|af| af < MAX_AF)
}

/// Six significant digits, so the table only changes with the data.
fn rounded(value: f64) -> f64 {
    format!("{value:.5e}").parse().unwrap_or(value)
}

/// Labelled variants of a bin.
#[derive(Debug, Clone, Copy, Default)]
struct Counts {
    lower: f64,
    pathogenic: f64,
    benign: f64,
}

impl Counts {
    fn total(&self) -> f64 {
        self.pathogenic + self.benign
    }

    fn merge(self, next: Counts) -> Counts {
        Counts {
            lower: self.lower,
            pathogenic: self.pathogenic + next.pathogenic,
            benign: self.benign + next.benign,
        }
    }
}

/// Bins of the labelled scores and their probabilities. The bins start with at least `MIN_BIN`
/// variants and `MIN_LABEL` of each label, never splitting equal scores, and are merged until the
/// pathogenic share no longer goes down while the score goes up. Which bins merge doesn't depend on
/// the prior.
fn curve(mut labelled: Vec<(f64, bool)>, group: Group) -> Result<Vec<(Counts, f64)>> {
    labelled.sort_by(|a, b| a.0.total_cmp(&b.0));
    let small = |bin: &Counts| bin.total() < MIN_BIN || bin.pathogenic.min(bin.benign) < MIN_LABEL;
    let mut bins: Vec<Counts> = Vec::new();
    let mut previous = None;
    for (score, pathogenic) in labelled {
        let variant = Counts {
            lower: score,
            pathogenic: f64::from(u8::from(pathogenic)),
            benign: f64::from(u8::from(!pathogenic)),
        };
        match bins.last_mut() {
            Some(bin) if small(bin) || previous == Some(score) => *bin = bin.merge(variant),
            _ => bins.push(variant),
        }
        previous = Some(score);
    }
    // A last bin that is too small joins the one before.
    if let Some(last) = bins.pop_if(|last| small(last)) {
        match bins.last_mut() {
            Some(before) => *before = before.merge(last),
            None => bins.push(last),
        }
    }

    let sum = |bins: &[Counts]| bins.iter().copied().fold(Counts::default(), Counts::merge);
    let all = sum(&bins);
    if all.pathogenic == 0.0 || all.benign == 0.0 {
        bail!(
            "{} pathogenic and {} benign variants",
            all.pathogenic,
            all.benign
        );
    }
    // One pseudo-variant per bin, split by the shares of the group, so that no bin gets a
    // probability of exactly 0 or 1. Adding it can break the order again.
    let pseudo = all.pathogenic / all.total();
    let bins = pool(bins, |bin| bin.pathogenic / bin.total());
    let bins = pool(bins, |bin| (bin.pathogenic + pseudo) / (bin.total() + 1.0))
        .into_iter()
        .map(|bin| Counts {
            pathogenic: bin.pathogenic + pseudo,
            benign: bin.benign + 1.0 - pseudo,
            ..bin
        })
        .collect_vec();

    let all = sum(&bins);
    Ok(bins
        .into_iter()
        .map(|bin| {
            let in_pathogenic = group.prior() * bin.pathogenic / all.pathogenic;
            let in_benign = (1.0 - group.prior()) * bin.benign / all.benign;
            (bin, in_pathogenic / (in_pathogenic + in_benign))
        })
        .collect())
}

/// Merges neighbouring bins wherever the share goes down (pool adjacent violators).
fn pool(bins: Vec<Counts>, share: impl Fn(&Counts) -> f64) -> Vec<Counts> {
    let mut pooled: Vec<Counts> = Vec::new();
    for mut bin in bins {
        while let Some(before) = pooled.pop_if(|before| share(before) > share(&bin)) {
            bin = before.merge(bin);
        }
        pooled.push(bin);
    }
    pooled
}

#[cfg(test)]
pub(crate) fn table(rows: &str) -> Calibration {
    let header = "score\tgroup\tclinvar\tsource\tprior\tlower\tpathogenic\tbenign\tprobability\n";
    Calibration::parse(&format!("{header}{rows}")).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn built_in_table() {
        // Loading checks that all curves are there.
        let calibration = Calibration::load(None).unwrap();
        for bins in calibration.0.values() {
            assert!(bins
                .iter()
                .tuple_windows()
                .all(|(a, b)| a.probability <= b.probability));
            // A new prior needs a new calibration.
            assert!(bins.iter().all(|bin| bin.prior == bin.group.prior()));
        }
        // A new CADD database needs a new calibration.
        assert!(calibration
            .0
            .values()
            .flatten()
            .all(|bin| bin.score != Score::Cadd || bin.source == Score::Cadd.database()));
    }

    #[test]
    fn curve_merges_bins() {
        // 120 variants each at 0, 1, 2 and 3. The pathogenic share drops from score 0 to 1, so
        // they merge, and 3 has no benign variants, so it joins 2.
        let labelled = (0..480)
            .map(|i| {
                let pathogenic = [i % 2 == 0, i % 4 == 0, i % 4 != 0, true][i / 120];
                ((i / 120) as f64, pathogenic)
            })
            .collect_vec();

        let bins = curve(labelled, Group::Missense).unwrap();

        let lower = bins.iter().map(|(counts, _)| counts.lower).collect_vec();
        assert_eq!(lower, [0.0, 2.0]);
        assert!(bins[0].1 < bins[1].1);
    }

    #[test]
    fn lookup() {
        let calibration = table(
            "revel\tmissense\t1\tt\t0\t0.2\t1\t9\t0.1\n\
             revel\tmissense\t1\tt\t0\t0.5\t5\t5\t0.5\n",
        );
        // Below, at and above the bins.
        for (value, probability) in [(0.1, 0.1), (0.5, 0.5), (0.9, 0.5)] {
            assert_eq!(
                calibration.probability(Score::Revel, Group::Missense, value),
                Some(probability)
            );
        }
    }
}
