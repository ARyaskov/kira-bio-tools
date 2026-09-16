//! FxHash exactly as the `fxhash` crate (rustc's hasher) computes it. The
//! 64-bit variant keys the entries of every `.ani` index on disk, so the
//! algorithm is pinned here instead of depending on an unmaintained crate,
//! and a test checks it against `fxhash` on random input.

use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hash, Hasher};

const SEED: u64 = 0x51_7c_c1_b7_27_22_0a_95;

/// The `fxhash` crate's `FxHasher64`: rotate-xor-multiply over 8-byte words.
#[derive(Clone, Copy, Debug, Default)]
pub struct FxHasher64 {
    hash: u64,
}

impl FxHasher64 {
    #[inline]
    fn add(&mut self, word: u64) {
        self.hash = (self.hash.rotate_left(5) ^ word).wrapping_mul(SEED);
    }
}

impl Hasher for FxHasher64 {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut b = bytes;
        while b.len() >= 8 {
            self.add(u64::from_le_bytes(b[..8].try_into().unwrap()));
            b = &b[8..];
        }
        if b.len() >= 4 {
            self.add(u32::from_le_bytes(b[..4].try_into().unwrap()) as u64);
            b = &b[4..];
        }
        for &x in b {
            self.add(x as u64);
        }
    }

    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(i as u64);
    }

    #[inline]
    fn write_u16(&mut self, i: u16) {
        self.add(i as u64);
    }

    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(i as u64);
    }

    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i);
    }

    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.hash
    }
}

/// The in-memory hasher; on 64-bit targets the crate's `FxHasher` is this.
pub type FxHasher = FxHasher64;
pub type FxBuildHasher = BuildHasherDefault<FxHasher64>;
pub type FxHashMap<K, V> = HashMap<K, V, FxBuildHasher>;
pub type FxHashSet<T> = HashSet<T, FxBuildHasher>;

/// The crate's `hash64`: the 64-bit hash of anything `Hash`, including the
/// length prefix `[u8]` writes.
pub fn hash64<T: Hash + ?Sized>(v: &T) -> u64 {
    let mut h = FxHasher64::default();
    v.hash(&mut h);
    h.finish()
}

#[cfg(test)]
#[path = "../tests/unit/fx.rs"]
mod tests;
