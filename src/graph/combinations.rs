//! Most and least severe combinations of haplotypes across phase blocks. Blocks only interact
//! through the frame and whether an earlier stop codon ended translation, so larger cases use a
//! dynamic program over that state. It is exact unless coding indels or a lost start make the
//! protein be realigned.

use anyhow::Result;
use itertools::{iproduct, Itertools};
use std::collections::{BTreeMap, HashMap};

const MAX_EXHAUSTIVE: usize = 64;

/// Frame offset 0, 1 or 2 when entering a block, or [`STOPPED`].
pub(crate) type State = usize;
pub(crate) const STOPPED: State = 3;
const STATES: usize = 4;

/// One haplotype per block, `None` for blocks not chosen yet, which are scored at the reference.
pub(crate) type Combination = Vec<Option<usize>>;
type ByState = [Option<Combination>; STATES];

/// The dynamic program assumes that, given the state a block is entered in, a more severe part for
/// the block never makes the whole combination less severe, as for sums, maxima and probabilistic
/// ORs.
pub(crate) trait Scoring {
    type Impact: Copy + 'static;
    const IMPACTS: &'static [Self::Impact];

    /// State after `block` for the haplotypes chosen in `prefix`. May only depend on the haplotype
    /// of `block` and the state before it.
    fn state_after(&mut self, prefix: &[Option<usize>], block: usize) -> Result<State>;

    fn impact(
        &mut self,
        combination: &[Option<usize>],
        impact: Self::Impact,
    ) -> Result<Option<f64>>;

    /// Whether the impact only counts in the original frame. Its minimum then has to enter the
    /// block in frame if possible, otherwise an upstream frameshift makes every haplotype harmless.
    fn needs_frame(&self, impact: Self::Impact) -> bool;
}

/// For every block haplotype and impact, the combinations containing it with the highest and the
/// lowest value. `sizes[b]` is the number of haplotypes of block `b`, in transcript order.
pub(crate) fn extreme_combinations<S: Scoring>(
    sizes: &[usize],
    scoring: &mut S,
) -> Result<Vec<Combination>> {
    let exhaustive = sizes
        .iter()
        .try_fold(1usize, |combinations, &size| combinations.checked_mul(size))
        .is_some_and(|combinations| combinations <= MAX_EXHAUSTIVE);
    let mut search = Search::new(sizes, scoring);
    let mut combinations = Vec::new();
    for impact in 0..S::IMPACTS.len() {
        search.score_candidates(impact, exhaustive)?;
        for maximise in [true, false] {
            combinations.extend(search.select(impact, maximise)?);
        }
    }
    Ok(combinations.into_iter().unique().collect())
}

struct Search<'a, S> {
    sizes: &'a [usize],
    scoring: &'a mut S,
    // Cache of state_after by block, haplotype and state before.
    transitions: HashMap<(usize, Option<usize>, State), State>,
    values: Vec<BTreeMap<Combination, Option<f64>>>,
}

impl<'a, S: Scoring> Search<'a, S> {
    fn new(sizes: &'a [usize], scoring: &'a mut S) -> Self {
        Search {
            sizes,
            scoring,
            transitions: HashMap::new(),
            values: vec![BTreeMap::new(); S::IMPACTS.len()],
        }
    }

    fn score_candidates(&mut self, impact: usize, exhaustive: bool) -> Result<()> {
        if !exhaustive {
            return self.dp(impact);
        }
        let all = self
            .sizes
            .iter()
            .map(|&size| (0..size).map(Some))
            .multi_cartesian_product()
            .collect_vec();
        for combination in all {
            self.value(&combination, impact)?;
        }
        Ok(())
    }

    /// Scores every block haplotype between the least severe prefix and suffix, and between the
    /// most severe ones, for each state.
    fn dp(&mut self, impact: usize) -> Result<()> {
        let least = self.prefixes(impact, false)?;
        let most = self.prefixes(impact, true)?;
        // Suffixes are ranked behind the least severe prefix, a severe one could mask them.
        let least_after = self.suffixes(&least, impact, false)?;
        let most_after = self.suffixes(&least, impact, true)?;
        for (prefixes, suffixes) in [(&least, &least_after), (&most, &most_after)] {
            for (block, state) in iproduct!(0..self.sizes.len(), 0..STATES) {
                for haplotype in 0..self.sizes[block] {
                    if let Some(combination) =
                        self.join(prefixes, suffixes, block, state, haplotype)?
                    {
                        self.value(&combination, impact)?;
                    }
                }
            }
        }
        Ok(())
    }

    /// Best prefix per block and the state it enters the block in.
    fn prefixes(&mut self, impact: usize, maximise: bool) -> Result<Vec<ByState>> {
        let blocks = self.sizes.len();
        let mut prefixes: Vec<ByState> = vec![Default::default(); blocks + 1];
        prefixes[0][0] = Some(vec![None; blocks]);
        for (block, state) in iproduct!(0..blocks, 0..STATES) {
            let Some(prefix) = prefixes[block][state].clone() else {
                continue;
            };
            for haplotype in 0..self.sizes[block] {
                let mut candidate = prefix.clone();
                candidate[block] = Some(haplotype);
                let next = self.transition(&candidate, block, state)?;
                self.keep_extreme(&mut prefixes[block + 1][next], candidate, impact, maximise)?;
            }
        }
        Ok(prefixes)
    }

    /// Best suffix per block and entry state, ranked behind `prefixes`.
    fn suffixes(
        &mut self,
        prefixes: &[ByState],
        impact: usize,
        maximise: bool,
    ) -> Result<Vec<ByState>> {
        let blocks = self.sizes.len();
        let mut suffixes: Vec<ByState> = vec![Default::default(); blocks + 1];
        suffixes[blocks] = std::array::from_fn(|_| Some(vec![None; blocks]));
        for (block, state) in iproduct!((0..blocks).rev(), 0..STATES) {
            let mut best = None;
            for haplotype in 0..self.sizes[block] {
                if let Some(candidate) = self.join(prefixes, &suffixes, block, state, haplotype)? {
                    self.keep_extreme(&mut best, candidate, impact, maximise)?;
                }
            }
            suffixes[block][state] = best.map(|mut combination| {
                combination[..block].fill(None);
                combination
            });
        }
        Ok(suffixes)
    }

    /// Prefix, haplotype and the matching suffix, if both exist.
    fn join(
        &mut self,
        prefixes: &[ByState],
        suffixes: &[ByState],
        block: usize,
        state: State,
        haplotype: usize,
    ) -> Result<Option<Combination>> {
        let Some(prefix) = prefixes[block][state].as_ref() else {
            return Ok(None);
        };
        let mut head = prefix.clone();
        head[block] = Some(haplotype);
        let next = self.transition(&head, block, state)?;
        let Some(suffix) = suffixes[block + 1][next].as_ref() else {
            return Ok(None);
        };
        Ok(Some(
            head.iter()
                .zip(suffix)
                .map(|(before, after)| before.or(*after))
                .collect(),
        ))
    }

    /// For every block haplotype, the most (or least) severe scored combination containing it.
    fn select(&mut self, impact: usize, maximise: bool) -> Result<Vec<Combination>> {
        let needs_frame = !maximise && self.scoring.needs_frame(S::IMPACTS[impact]);
        let mut scored = Vec::new();
        for (combination, value) in self.values[impact].clone() {
            if combination.iter().all(Option::is_some) {
                // Entering a block out of frame ranks higher, so the minimum avoids it if possible.
                let shifted = if needs_frame {
                    self.shifted(&combination)?
                } else {
                    vec![false; combination.len()]
                };
                scored.push((combination, value, shifted));
            }
        }

        let mut extremes = Vec::new();
        for (block, &size) in self.sizes.iter().enumerate() {
            for haplotype in 0..size {
                let extreme = scored
                    .iter()
                    .filter(|(combination, ..)| combination[block] == Some(haplotype))
                    .reduce(|old, new| {
                        let (_, old_value, old_shifted) = old;
                        let (_, new_value, new_shifted) = new;
                        let replace = more_extreme(
                            (new_shifted[block], new_value),
                            (old_shifted[block], old_value),
                            maximise,
                        );
                        if replace {
                            new
                        } else {
                            old
                        }
                    });
                extremes.extend(extreme.map(|(combination, ..)| combination.clone()));
            }
        }
        Ok(extremes)
    }

    fn transition(
        &mut self,
        prefix: &[Option<usize>],
        block: usize,
        before: State,
    ) -> Result<State> {
        if before == STOPPED {
            return Ok(STOPPED);
        }
        let key = (block, prefix[block], before);
        if let Some(&state) = self.transitions.get(&key) {
            return Ok(state);
        }
        let state = self.scoring.state_after(prefix, block)?;
        self.transitions.insert(key, state);
        Ok(state)
    }

    /// Whether a complete combination enters each block out of frame or stopped.
    fn shifted(&mut self, combination: &[Option<usize>]) -> Result<Vec<bool>> {
        let mut prefix = vec![None; combination.len()];
        let mut state = 0;
        let mut shifted = Vec::with_capacity(combination.len());
        for (block, &haplotype) in combination.iter().enumerate() {
            shifted.push(state != 0);
            prefix[block] = haplotype;
            state = self.transition(&prefix, block, state)?;
        }
        Ok(shifted)
    }

    /// Puts the candidate into `slot` if it is empty or the candidate is more extreme.
    fn keep_extreme(
        &mut self,
        slot: &mut Option<Combination>,
        candidate: Combination,
        impact: usize,
        maximise: bool,
    ) -> Result<()> {
        let new = self.value(&candidate, impact)?;
        let replace = match slot {
            Some(old) => more_extreme(new, self.value(old, impact)?, maximise),
            None => true,
        };
        if replace {
            *slot = Some(candidate);
        }
        Ok(())
    }

    fn value(&mut self, combination: &[Option<usize>], impact: usize) -> Result<Option<f64>> {
        if let Some(&value) = self.values[impact].get(combination) {
            return Ok(value);
        }
        let value = self.scoring.impact(combination, S::IMPACTS[impact])?;
        self.values[impact].insert(combination.to_vec(), value);
        Ok(value)
    }
}

/// Strictly higher, or lower when minimising. Missing values rank lowest.
fn more_extreme<T: PartialOrd>(new: T, old: T, maximise: bool) -> bool {
    if maximise {
        new > old
    } else {
        new < old
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    // Per block and haplotype the frameshift (3 stops translation) and the effect in each frame.
    struct Toy(Vec<Vec<(usize, [f64; 3])>>);

    impl Toy {
        fn translate(&self, combination: &[Option<usize>]) -> (State, Vec<f64>) {
            let mut state = 0;
            let mut effects = Vec::new();
            for (haplotypes, haplotype) in self.0.iter().zip(combination) {
                if let (Some(haplotype), true) = (haplotype, state != STOPPED) {
                    let (shift, effect) = haplotypes[*haplotype];
                    effects.push(effect[state]);
                    state = if shift == 3 {
                        STOPPED
                    } else {
                        (state + shift) % 3
                    };
                }
            }
            (state, effects)
        }
    }

    // Maximum and probabilistic OR of the effects, the latter only counting in frame.
    impl Scoring for Toy {
        type Impact = usize;
        const IMPACTS: &'static [usize] = &[0, 1];

        fn state_after(&mut self, prefix: &[Option<usize>], block: usize) -> Result<State> {
            Ok(self.translate(&prefix[..=block]).0)
        }

        fn impact(&mut self, combination: &[Option<usize>], impact: usize) -> Result<Option<f64>> {
            let (_, effects) = self.translate(combination);
            Ok(match impact {
                0 => effects.into_iter().reduce(f64::max),
                _ => Some(1.0 - effects.iter().map(|effect| 1.0 - effect).product::<f64>()),
            })
        }

        fn needs_frame(&self, impact: usize) -> bool {
            impact == 1
        }
    }

    #[test]
    fn dp_matches_brute_force() {
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..200 {
            let mut blocks = Vec::new();
            for _ in 0..rng.gen_range(1..7) {
                let mut haplotypes = Vec::new();
                for _ in 0..rng.gen_range(1..4) {
                    let effect = std::array::from_fn(|_| rng.gen_range(0..5) as f64 / 4.0);
                    haplotypes.push((rng.gen_range(0..4), effect));
                }
                blocks.push(haplotypes);
            }
            let sizes = blocks.iter().map(Vec::len).collect_vec();
            let mut toy = Toy(blocks);
            for impact in 0..Toy::IMPACTS.len() {
                for maximise in [true, false] {
                    let mut extremes = |exhaustive| {
                        let mut search = Search::new(&sizes, &mut toy);
                        search.score_candidates(impact, exhaustive).unwrap();
                        let combinations = search.select(impact, maximise).unwrap();
                        combinations
                            .iter()
                            .map(|combination| search.value(combination, impact).unwrap())
                            .collect_vec()
                    };
                    let dp = extremes(false);
                    assert_eq!(dp, extremes(true), "{sizes:?} {impact} {maximise}");
                }
            }
        }
    }

    #[test]
    fn minimum_stays_in_frame() {
        // Block 0 can shift the frame, behind which block 1 has no effect.
        let mut toy = Toy(vec![
            vec![(0, [0.0; 3]), (1, [0.0; 3])],
            vec![(0, [0.5, 0.0, 0.0])],
        ]);
        let mut search = Search::new(&[2, 1], &mut toy);
        search.score_candidates(0, true).unwrap();
        search.score_candidates(1, true).unwrap();

        // Index 2 is the haplotype of block 1.
        assert_eq!(search.select(0, false).unwrap()[2], vec![Some(1), Some(0)]);
        assert_eq!(search.select(1, false).unwrap()[2], vec![Some(0), Some(0)]);
    }
}
