use crate::annotation::{annotate_coding_variants, Scores, HUB_DATABASES};
use crate::cli::{Command, GeneBeCredentials, Haplodon};
use crate::graph::duck::{create_scores, read_scores, write_graphs, write_scores};
use crate::graph::peptide::write_peptides;
use crate::graph::transcript::transcripts;
use crate::graph::VariantGraph;
use crate::show::render_scores;
use crate::utils::bcf::get_targets;
use crate::utils::create_output_dir;
use anyhow::{bail, Result};
use bio::stats::{LogProb, Prob};
use clap::Parser;
use env_logger::Env;
use genebears::{ClientConfig, GeneBears, Genome, Hub, Store};
use log::{info, warn};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::sync::{Arc, Mutex};

mod annotate;
mod annotation;
mod cli;
mod graph;
mod show;
mod transcription;
mod translation;
mod utils;

fn main() -> Result<()> {
    let args = Haplodon::parse();
    let log_level = if args.verbose { "info" } else { "off" };
    env_logger::Builder::from_env(Env::default().default_filter_or(log_level)).init();
    if let Some(n) = args.threads {
        ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to build thread pool");
    }
    args.command.run()
}

impl Command {
    fn run(&self) -> Result<()> {
        match self {
            Command::Build {
                calls,
                observations,
                min_prob_present,
                min_vaf,
                output,
            } => {
                let targets = get_targets(calls)?;
                let graphs = targets
                    .into_par_iter()
                    .filter_map(|target| {
                        info!("Building graph for target {target}");
                        let variant_graph = VariantGraph::build(
                            calls,
                            observations,
                            &target,
                            LogProb::from(Prob(*min_prob_present)),
                            *min_vaf,
                        )
                        .ok()?;
                        info!(
                            "Finished building graph for target {} with {} nodes",
                            target,
                            variant_graph.graph.node_count()
                        );
                        if !variant_graph.is_empty() {
                            Some((target, variant_graph))
                        } else {
                            None
                        }
                    })
                    .collect();
                write_graphs(graphs, output)?;
            }
            Command::Process {
                features,
                reference,
                graph,
                distance_metric,
                haplotype_metric,
                output,
                max_haplotypes_per_block,
                genebe_cache,
                genebe_email,
                genebe_api_key,
                genome_build,
                genebe_store,
            } => {
                info!("Reading reference genome from {reference:?}");
                let reference_genome = utils::fasta::read_reference(reference);
                // Scores are looked up by position, so the reference has to be of the genome build.
                let chr1 = reference_genome.get("1").or(reference_genome.get("chr1"));
                let length = match genome_build {
                    Genome::Hg38 => 248_956_422,
                    Genome::Hg19 => 249_250_621,
                    Genome::T2t => 248_387_328,
                };
                if chr1.is_some_and(|chr1| chr1.len() != length) {
                    bail!(
                        "The reference is not {}, set --genome-build",
                        genome_build.as_str()
                    );
                }
                let write_lock = Arc::new(Mutex::new(()));
                let credentials = match (genebe_email, genebe_api_key) {
                    (Some(email), Some(api_key)) => Some(GeneBeCredentials {
                        email: email.clone(),
                        api_key: api_key.clone(),
                    }),
                    _ => None,
                };
                let mut config = ClientConfig::default();
                if let Some(cache) = genebe_cache {
                    config = config.with_cache(cache);
                }
                if let Some(credentials) = &credentials {
                    config.email = Some(credentials.email.clone());
                    config.api_key = Some(credentials.api_key.clone());
                }
                let scores = Scores::new(
                    *genome_build,
                    Store::new(genebe_store),
                    GeneBears::new(config)?,
                )?;
                let transcripts = transcripts(features, graph)?;
                let annotations = annotate_coding_variants(graph, &transcripts, &scores)?;
                create_scores(output)?;
                let total = transcripts.len();
                transcripts.into_par_iter().enumerate().try_for_each(
                    |(i, transcript)| -> anyhow::Result<()> {
                        info!(
                            "Processing transcript {} ({}/{})",
                            transcript.name(),
                            i + 1,
                            total
                        );
                        if !transcript.is_translatable(&reference_genome) {
                            if let Some(base) =
                                transcript.untranslatable_reference_base(&reference_genome)
                            {
                                warn!(
                                    "Skipping transcript {}: reference coding sequence contains non-ACGT base '{}'",
                                    transcript.name(),
                                    base as char
                                );
                            }
                            return Ok(());
                        }
                        let (blocks, samples) =
                            transcript.block_haplotypes(graph, *max_haplotypes_per_block)?;
                        let scores = transcript.score_haplotypes(
                            &blocks,
                            &samples,
                            &reference_genome,
                            *haplotype_metric,
                            *distance_metric,
                            &annotations[&transcript.target],
                        )?;
                        info!(
                            "Writing scores for {} different haplotypes for transcript {}",
                            scores.len(),
                            transcript.name()
                        );
                        let _lock = write_lock.lock().unwrap();
                        write_scores(output, scores, transcript)?;
                        Ok(())
                    },
                )?;
            }
            Command::Cache {
                genebe_email,
                genebe_api_key,
                genebe_store,
            } => {
                let hub = Hub::new(genebe_email, genebe_api_key)?;
                let store = Store::new(genebe_store);
                let rt = tokio::runtime::Runtime::new()?;
                for (id, _, _) in HUB_DATABASES {
                    rt.block_on(hub.pull(&id.parse()?, &store))?;
                }
            }
            Command::Peptides {
                features,
                reference,
                graph,
                interval,
                sample,
                output,
                events,
                min_event_prob,
                background_events,
                max_background_event_prob,
                max_cds_length: _max_cds_length,
            } => {
                info!("Reading reference genome from {reference:?}");
                let reference_genome = utils::fasta::read_reference(reference);
                let mut peptides = Vec::new();
                for transcript in transcripts(features, graph)? {
                    info!("Processing transcript {}", transcript.name());
                    let transcript_peptides = transcript.peptides(
                        graph,
                        &reference_genome,
                        interval.clone(),
                        sample,
                        events,
                        LogProb::from(Prob(*min_event_prob)),
                        background_events,
                        LogProb::from(Prob(*max_background_event_prob)),
                    )?;
                    peptides.extend(transcript_peptides);
                }
                info!("Writing peptides to {output:?}");
                write_peptides(peptides, output)?;
            }
            Command::Plot {
                input,
                notation,
                report_protein,
                output,
            } => {
                if let Some(parent) = output.parent() {
                    create_output_dir(parent)?;
                }
                let scores = read_scores(input, *notation)?;
                render_scores(output, &scores, *report_protein)?;
            }
            Command::Annotate {
                calls,
                scores,
                output,
            } => {
                if let Some(parent) = output.parent() {
                    create_output_dir(parent)?;
                }
                annotate::annotate(calls, scores, output)?;
            }
        }
        Ok(())
    }
}
