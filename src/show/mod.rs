use crate::graph::score::ScoreRecord;
use anyhow::anyhow;
use anyhow::Result;
use csv::WriterBuilder;
use itertools::Itertools;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Writes given per-haplotype scores into a single TSV file with the columns:
/// `transcript`, `score`, `haplotype`, annotation score columns (`revel_score`,
/// `acmg_score`, `spliceai_score`, `alphamissense_score`) and two columns
/// per sample containing the frequency and supporting read information. Each row
/// corresponds to one (score, per-sample) pair for a given transcript.
///
/// # Arguments
///
/// * `output_path` - A reference to a `PathBuf` that holds the path of the output TSV file.
/// * `scores` - A reference to a `HashMap` that holds the scores, likelihoods, and haplotypes for each transcript.
/// * `report_protein` - Whether to include the alternative protein sequence as an additional column.
pub(crate) fn render_scores(
    output_path: &PathBuf,
    scores: &HashMap<String, Vec<ScoreRecord>>,
    report_protein: bool,
) -> Result<()> {
    let mut wtr = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(Path::new(output_path))?;

    let samples: Vec<String> = scores
        .values()
        .flat_map(|scores| {
            scores
                .iter()
                .flat_map(|(_, sample_scores, _, _, _, _)| sample_scores.keys().cloned())
        })
        .unique()
        .collect();

    let mut headers = vec![
        "transcript".to_string(),
        "score".to_string(),
        "haplotype".to_string(),
        "revel_score".to_string(),
        "acmg_score".to_string(),
        "spliceai_score".to_string(),
        "alphamissense_score".to_string(),
    ];
    if report_protein {
        headers.push("protein".to_string());
    }
    headers.extend(samples.iter().map(|s| format!("{}:frequency", s)));
    headers.extend(samples.iter().map(|s| format!("{}:supporting_reads", s)));
    wtr.write_record(headers)?;

    for (transcript, hap_scores) in scores {
        for (score_val, sample_scores, haplotype, supporting_reads, annotation, protein) in
            hap_scores
        {
            let mut row = vec![
                transcript.to_string(),
                score_val.to_string(),
                haplotype.to_string(),
                annotation
                    .revel_score
                    .map_or("".to_string(), |s| s.to_string()),
                annotation
                    .acmg_score
                    .map_or("".to_string(), |s| s.to_string()),
                annotation
                    .spliceai_score
                    .map_or("".to_string(), |s| s.to_string()),
                annotation
                    .alphamissense_score
                    .map_or("".to_string(), |s| s.to_string()),
            ];
            if report_protein {
                row.push(protein.to_string());
            }
            for sample in &samples {
                let val = sample_scores.get(sample).ok_or_else(|| {
                    anyhow!(
                        "No score found for sample '{}' in transcript '{}'",
                        sample,
                        transcript
                    )
                })?;
                row.push(val.to_string());
            }
            for sample in &samples {
                let reads = supporting_reads
                    .iter()
                    .map(|m| m.get(sample).unwrap_or(&0))
                    .cloned()
                    .join(";");
                row.push(reads);
            }
            wtr.write_record(row)?;
        }
    }

    wtr.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::Annotation;

    #[test]
    fn test_render_scores() {
        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.keep().join("scores.tsv");
        let annotion = Annotation {
            revel_score: Some(0.8),
            acmg_score: Some(0.9),
            spliceai_score: Some(0.7),
            alphamissense_score: Some(0.6),
            gnomad_frequencies: HashMap::new(),
        };
        let scores = HashMap::from([(
            "chr1:some feature".to_string(),
            vec![(
                0.04f64,
                HashMap::from([("A".to_string(), 0.1f32), ("C".to_string(), 0.2f32)]),
                "c.[100A>G;105C>T]".to_string(),
                vec![HashMap::from([
                    ("A".to_string(), 10u32),
                    ("C".to_string(), 5u32),
                ])],
                annotion,
                "FL".to_string(),
            )],
        )]);
        render_scores(&output_path, &scores, true).unwrap();
        let content = std::fs::read_to_string(&output_path).unwrap();
        assert!(output_path.exists());
        assert!(content.contains("transcript"));
        assert!(content.contains("chr1:some feature"));
        assert!(content.contains("10"));
        assert!(content.contains("protein"));
        assert!(content.contains("FL"));
    }
}
