//! Hidden Markov model kernel shared by `roh` and `cnv`: a port of bcftools
//! `HMM.c`. A transition matrix is laid out `P(i|j) = m[n*i + j]` (transition
//! j → i); powers of it are precomputed so a gap of `d` positions costs one
//! matrix copy plus `d / ntprob` multiplications. Floating-point operation
//! order follows the C code so posteriors agree to the last digit.

/// Per-site adjustment of the transition matrix (`set_tprob_f`).
pub trait TransitionFn {
    fn apply(&mut self, prev_pos: u32, pos: u32, tprob: &mut [f64]);
}

impl TransitionFn for () {
    fn apply(&mut self, _prev_pos: u32, _pos: u32, _tprob: &mut [f64]) {}
}

/// Viterbi/forward state captured at a position, to resume a run after a
/// buffer flush (`hmm_snapshot`).
#[derive(Clone, Debug)]
pub struct Snapshot {
    pub snap_at_pos: u32,
    vit_prob: Vec<f64>,
    fwd_prob: Vec<f64>,
}

#[derive(Clone, Debug)]
struct State {
    snap_at_pos: u32,
    vit_prob: Vec<f64>,
    fwd_prob: Vec<f64>,
    bwd_prob: Vec<f64>,
}

pub struct Hmm<T: TransitionFn = ()> {
    nstates: usize,
    vprob: Vec<f64>,
    vprob_tmp: Vec<f64>,
    /// `[nstates * n]`: back-pointers during the pass, then the decoded state
    /// of site `i` at `vpath[i * nstates]`.
    vpath: Vec<u8>,
    bwd: Vec<f64>,
    bwd_tmp: Vec<f64>,
    /// `[nstates * (n + 1)]`: row 0 is the initial distribution, row `i + 1`
    /// the posterior of site `i` after `run_fwd_bwd`.
    fwd: Vec<f64>,
    ntprob_arr: usize,
    curr_tprob: Vec<f64>,
    tmp: Vec<f64>,
    tprob_arr: Vec<f64>,
    pub tfn: T,
    init: State,
    state: State,
    snapshot: Option<Snapshot>,
}

/// `dst = a × b` for `n × n` matrices in the `P(i|j)` layout.
fn mat_mul(n: usize, a: &[f64], b: &[f64], dst: &mut [f64]) {
    for i in 0..n {
        for j in 0..n {
            let mut val = 0.0;
            for k in 0..n {
                val += a[n * i + k] * b[n * k + j];
            }
            dst[n * i + j] = val;
        }
    }
}

impl<T: TransitionFn> Hmm<T> {
    /// `hmm_init`: `tprob` is the one-step matrix; `ntprob` powers are cached.
    pub fn new(nstates: usize, tprob: &[f64], ntprob: usize, tfn: T) -> Self {
        let blank = || State {
            snap_at_pos: 0,
            vit_prob: vec![0.0; nstates],
            fwd_prob: vec![0.0; nstates],
            bwd_prob: vec![0.0; nstates],
        };
        let mut h = Self {
            nstates,
            vprob: vec![0.0; nstates],
            vprob_tmp: vec![0.0; nstates],
            vpath: Vec::new(),
            bwd: vec![0.0; nstates],
            bwd_tmp: vec![0.0; nstates],
            fwd: Vec::new(),
            ntprob_arr: 0,
            curr_tprob: vec![0.0; nstates * nstates],
            tmp: vec![0.0; nstates * nstates],
            tprob_arr: Vec::new(),
            tfn,
            init: blank(),
            state: blank(),
            snapshot: None,
        };
        h.set_tprob(tprob, ntprob);
        h.init_states(None);
        h
    }

    pub fn nstates(&self) -> usize {
        self.nstates
    }

    /// The one-step transition matrix (`hmm_get_tprob`).
    pub fn tprob(&self) -> &[f64] {
        &self.tprob_arr[..self.nstates * self.nstates]
    }

    /// `hmm_set_tprob`: install a one-step matrix and precompute its powers.
    pub fn set_tprob(&mut self, tprob: &[f64], ntprob: usize) {
        let n = self.nstates;
        let nn = n * n;
        self.ntprob_arr = ntprob;
        let cnt = ntprob.max(1);
        if self.tprob_arr.len() < nn * cnt {
            self.tprob_arr.resize(nn * cnt, 0.0);
        }
        self.tprob_arr[..nn].copy_from_slice(&tprob[..nn]);
        for i in 1..cnt {
            let (head, tail) = self.tprob_arr.split_at_mut(i * nn);
            let (a, b) = (&head[..nn], &head[(i - 1) * nn..i * nn]);
            mat_mul(n, a, b, &mut tail[..nn]);
        }
    }

    /// `hmm_init_states`: initial distribution (uniform when `None`).
    pub fn init_states(&mut self, probs: Option<&[f64]>) {
        let n = self.nstates;
        self.init.snap_at_pos = 0;
        self.state.snap_at_pos = 0;
        match probs {
            Some(p) => {
                self.init.vit_prob.copy_from_slice(&p[..n]);
                let sum: f64 = self.init.vit_prob.iter().sum();
                for v in &mut self.init.vit_prob {
                    *v /= sum;
                }
            }
            None => {
                for v in &mut self.init.vit_prob {
                    *v = 1.0 / n as f64;
                }
            }
        }
        self.init.fwd_prob.copy_from_slice(&self.init.vit_prob);
        self.init.bwd_prob.copy_from_slice(&self.init.vit_prob);
        self.state.vit_prob.copy_from_slice(&self.init.vit_prob);
        self.state.fwd_prob.copy_from_slice(&self.init.fwd_prob);
        self.state.bwd_prob.copy_from_slice(&self.init.bwd_prob);
    }

    /// `hmm_snapshot`: ask the next runs to record their state at `pos`.
    /// Take the filled snapshot back with [`Hmm::take_snapshot`].
    pub fn request_snapshot(&mut self, prev: Option<Snapshot>, pos: u32) {
        let n = self.nstates;
        let mut s = match prev {
            Some(s) if s.vit_prob.len() == n => s,
            _ => Snapshot {
                snap_at_pos: 0,
                vit_prob: vec![0.0; n],
                fwd_prob: vec![0.0; n],
            },
        };
        s.snap_at_pos = pos;
        self.snapshot = Some(s);
    }

    pub fn take_snapshot(&mut self) -> Option<Snapshot> {
        self.snapshot.take()
    }

    /// `hmm_restore`: start the next run from a snapshot, or from the initial
    /// distribution when there is none.
    pub fn restore(&mut self, snap: Option<&Snapshot>) {
        match snap {
            Some(s) if s.snap_at_pos != 0 => {
                self.state.snap_at_pos = s.snap_at_pos;
                self.state.vit_prob.copy_from_slice(&s.vit_prob);
                self.state.fwd_prob.copy_from_slice(&s.fwd_prob);
            }
            _ => {
                self.state.snap_at_pos = 0;
                self.state.vit_prob.copy_from_slice(&self.init.vit_prob);
                self.state.fwd_prob.copy_from_slice(&self.init.fwd_prob);
            }
        }
    }

    /// `hmm_reset`: forget the snapshot and restart from the initial state.
    pub fn reset(&mut self, snap: Option<&mut Snapshot>) {
        if let Some(s) = snap {
            s.snap_at_pos = 0;
        }
        self.state.snap_at_pos = 0;
        self.state.vit_prob.copy_from_slice(&self.init.vit_prob);
        self.state.fwd_prob.copy_from_slice(&self.init.fwd_prob);
    }

    /// `_set_tprob`: `curr_tprob = M^(pos_diff + 1)`.
    fn set_curr_tprob(&mut self, pos_diff: i64) {
        assert!(pos_diff >= 0);
        let n = self.nstates;
        let nn = n * n;
        let ntp = self.ntprob_arr as i64;
        let idx = if ntp > 0 {
            (pos_diff % ntp) as usize
        } else {
            0
        };
        self.curr_tprob
            .copy_from_slice(&self.tprob_arr[idx * nn..(idx + 1) * nn]);
        if ntp > 0 {
            let blocks = pos_diff / ntp;
            let last = (self.ntprob_arr - 1) * nn;
            for _ in 0..blocks {
                mat_mul(
                    n,
                    &self.tprob_arr[last..last + nn],
                    &self.curr_tprob,
                    &mut self.tmp,
                );
                self.curr_tprob.copy_from_slice(&self.tmp);
            }
        }
    }

    fn step_tprob(&mut self, prev_pos: u32, pos: u32, pos_diff: i64) {
        self.set_curr_tprob(pos_diff);
        self.tfn.apply(prev_pos, pos, &mut self.curr_tprob);
    }

    /// `hmm_run_viterbi`: `eprobs` holds `nstates` emission probabilities per
    /// site; `sites` are the positions.
    pub fn run_viterbi(&mut self, eprobs: &[f64], sites: &[u32]) {
        let n = sites.len();
        let ns = self.nstates;
        if self.vpath.len() < n * ns {
            self.vpath.resize(n * ns, 0);
        }
        self.vprob.copy_from_slice(&self.state.vit_prob);
        let mut prev_pos = if self.state.snap_at_pos != 0 {
            self.state.snap_at_pos
        } else {
            sites[0]
        };
        for i in 0..n {
            let pos_diff = if sites[i] == prev_pos {
                0
            } else {
                sites[i] as i64 - prev_pos as i64 - 1
            };
            self.step_tprob(prev_pos, sites[i], pos_diff);
            prev_pos = sites[i];
            let eprob = &eprobs[i * ns..(i + 1) * ns];
            let mut vnorm = 0.0;
            for j in 0..ns {
                let mut vmax = 0.0;
                let mut k_vmax = 0usize;
                for k in 0..ns {
                    let pval = self.vprob[k] * self.curr_tprob[ns * j + k];
                    if vmax < pval {
                        vmax = pval;
                        k_vmax = k;
                    }
                }
                self.vpath[i * ns + j] = k_vmax as u8;
                self.vprob_tmp[j] = vmax * eprob[j];
                vnorm += self.vprob_tmp[j];
            }
            for j in 0..ns {
                self.vprob_tmp[j] /= vnorm;
            }
            std::mem::swap(&mut self.vprob, &mut self.vprob_tmp);
            if let Some(s) = &mut self.snapshot
                && sites[i] == s.snap_at_pos {
                    s.vit_prob.copy_from_slice(&self.vprob);
                }
        }
        let mut iptr = 0usize;
        for i in 1..ns {
            if self.vprob[iptr] < self.vprob[i] {
                iptr = i;
            }
        }
        for i in (0..n).rev() {
            iptr = self.vpath[i * ns + iptr] as usize;
            self.vpath[i * ns] = iptr as u8;
        }
    }

    /// Decoded Viterbi state of site `i` after [`Hmm::run_viterbi`].
    pub fn vstate(&self, i: usize) -> u8 {
        self.vpath[i * self.nstates]
    }

    /// `hmm_run_fwd_bwd`: posteriors land in [`Hmm::fwd_bwd_prob`].
    pub fn run_fwd_bwd(&mut self, eprobs: &[f64], sites: &[u32]) {
        let n = sites.len();
        let ns = self.nstates;
        if self.fwd.len() < (n + 1) * ns {
            self.fwd.resize((n + 1) * ns, 0.0);
        }
        self.fwd[..ns].copy_from_slice(&self.state.fwd_prob);
        self.bwd.copy_from_slice(&self.state.bwd_prob);
        let mut prev_pos = if self.state.snap_at_pos != 0 {
            self.state.snap_at_pos
        } else {
            sites[0]
        };
        for i in 0..n {
            let pos_diff = if sites[i] == prev_pos {
                0
            } else {
                sites[i] as i64 - prev_pos as i64 - 1
            };
            self.step_tprob(prev_pos, sites[i], pos_diff);
            prev_pos = sites[i];
            let eprob = &eprobs[i * ns..(i + 1) * ns];
            let mut norm = 0.0;
            for j in 0..ns {
                let mut pval = 0.0;
                for k in 0..ns {
                    pval += self.fwd[i * ns + k] * self.curr_tprob[ns * j + k];
                }
                self.fwd[(i + 1) * ns + j] = pval * eprob[j];
                norm += self.fwd[(i + 1) * ns + j];
            }
            for j in 0..ns {
                self.fwd[(i + 1) * ns + j] /= norm;
            }
            if let Some(s) = &mut self.snapshot
                && sites[i] == s.snap_at_pos {
                    s.fwd_prob
                        .copy_from_slice(&self.fwd[(i + 1) * ns..(i + 2) * ns]);
                }
        }
        let mut prev_pos = sites[n - 1];
        for i in 0..n {
            let isite = n - i - 1;
            let pos_diff = if sites[isite] == prev_pos {
                0
            } else {
                prev_pos as i64 - sites[isite] as i64 - 1
            };
            self.step_tprob(sites[isite], prev_pos, pos_diff);
            prev_pos = sites[isite];
            let eprob = &eprobs[isite * ns..(isite + 1) * ns];
            let mut bwd_norm = 0.0;
            for j in 0..ns {
                let mut pval = 0.0;
                for k in 0..ns {
                    pval += self.bwd[k] * eprob[k] * self.curr_tprob[ns * k + j];
                }
                self.bwd_tmp[j] = pval;
                bwd_norm += pval;
            }
            let row = (n - i) * ns;
            let mut norm = 0.0;
            for j in 0..ns {
                self.bwd_tmp[j] /= bwd_norm;
                self.fwd[row + j] *= self.bwd_tmp[j];
                norm += self.fwd[row + j];
            }
            for j in 0..ns {
                self.fwd[row + j] /= norm;
            }
            std::mem::swap(&mut self.bwd, &mut self.bwd_tmp);
        }
    }

    /// Rows of `nstates` values: row 0 the initial distribution, row `i + 1`
    /// the posterior of site `i`.
    pub fn fwd_bwd_prob(&self) -> &[f64] {
        &self.fwd
    }

    /// `hmm_run_baum_welch`: one re-estimation of the transition matrix;
    /// the result stays valid until the next run.
    pub fn run_baum_welch(&mut self, eprobs: &[f64], sites: &[u32]) -> &[f64] {
        let n = sites.len();
        let ns = self.nstates;
        if self.fwd.len() < (n + 1) * ns {
            self.fwd.resize((n + 1) * ns, 0.0);
        }
        self.fwd[..ns].copy_from_slice(&self.state.fwd_prob);
        self.bwd.copy_from_slice(&self.state.bwd_prob);
        let mut prev_pos = if self.state.snap_at_pos != 0 {
            self.state.snap_at_pos
        } else {
            sites[0]
        };
        let mut tmp_xi = vec![0.0f64; ns * ns];
        let mut tmp_gamma = vec![0.0f64; ns];
        let mut fwd_bwd = vec![0.0f64; ns];
        for i in 0..n {
            let pos_diff = if sites[i] == prev_pos {
                0
            } else {
                sites[i] as i64 - prev_pos as i64 - 1
            };
            self.step_tprob(prev_pos, sites[i], pos_diff);
            prev_pos = sites[i];
            let eprob = &eprobs[i * ns..(i + 1) * ns];
            let mut norm = 0.0;
            for j in 0..ns {
                let mut pval = 0.0;
                for k in 0..ns {
                    pval += self.fwd[i * ns + k] * self.curr_tprob[ns * j + k];
                }
                self.fwd[(i + 1) * ns + j] = pval * eprob[j];
                norm += self.fwd[(i + 1) * ns + j];
            }
            for j in 0..ns {
                self.fwd[(i + 1) * ns + j] /= norm;
            }
        }
        let mut prev_pos = sites[n - 1];
        for i in 0..n {
            let isite = n - i - 1;
            let pos_diff = if sites[isite] == prev_pos {
                0
            } else {
                prev_pos as i64 - sites[isite] as i64 - 1
            };
            self.step_tprob(sites[isite], prev_pos, pos_diff);
            prev_pos = sites[isite];
            let eprob = &eprobs[isite * ns..(isite + 1) * ns];
            let mut bwd_norm = 0.0;
            for j in 0..ns {
                let mut pval = 0.0;
                for k in 0..ns {
                    pval += self.bwd[k] * eprob[k] * self.curr_tprob[ns * k + j];
                }
                self.bwd_tmp[j] = pval;
                bwd_norm += pval;
            }
            let row = (n - i) * ns;
            let mut norm = 0.0;
            for j in 0..ns {
                self.bwd_tmp[j] /= bwd_norm;
                fwd_bwd[j] = self.fwd[row + j] * self.bwd_tmp[j];
                norm += fwd_bwd[j];
            }
            for j in 0..ns {
                fwd_bwd[j] /= norm;
                tmp_gamma[j] += fwd_bwd[j];
            }
            for j in 0..ns {
                for k in 0..ns {
                    tmp_xi[ns * k + j] +=
                        self.fwd[row + j] * self.bwd[k] * self.tprob_arr[ns * k + j] * eprob[k]
                            / norm;
                }
            }
            self.fwd[row..row + ns].copy_from_slice(&fwd_bwd[..ns]);
            std::mem::swap(&mut self.bwd, &mut self.bwd_tmp);
        }
        for j in 0..ns {
            let mut norm = 0.0;
            for k in 0..ns {
                self.curr_tprob[ns * k + j] = tmp_xi[ns * k + j] / tmp_gamma[j];
                norm += self.curr_tprob[ns * k + j];
            }
            for k in 0..ns {
                self.curr_tprob[ns * k + j] /= norm;
            }
        }
        &self.curr_tprob
    }
}

#[cfg(test)]
#[path = "../tests/unit/hmm.rs"]
mod tests;
