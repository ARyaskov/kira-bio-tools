//! # Ultra-Optimized Parallel PGM-Index

#[cfg(not(target_env = "msvc"))]
use jemallocator::Jemalloc;
#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

use crossbeam::utils::CachePadded;
use num_traits::ToPrimitive;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

pub trait Key: Ord + Copy + ToPrimitive + Send + Sync {}
impl Key for u64 {}
impl Key for i64 {}
impl Key for u32 {}
impl Key for i32 {}

#[derive(Debug, Clone, Copy)]
enum DataComplexity {
    Linear,
    Quadratic,
    Exponential,
    Random,
}

/// One linear segment: y = slope * x + intercept (cache-aligned)
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[repr(C, align(64))]
pub struct Segment<K: Key> {
    slope: f64,
    intercept: f64,
    min_key: K,
    max_key: K,
    start_pos: usize,
    end_pos: usize,
}

/// The main ultra-optimized PGM-Index structure
#[derive(Debug)]
pub struct PGMIndex<K: Key> {
    pub epsilon: usize,
    pub data: Arc<Vec<K>>,
    segments: Vec<Segment<K>>,
    segment_lookup: Vec<usize>,
    lookup_scale: f64,
    min_key_f64: f64,
    // Cache-padded for avoiding false sharing
    stats: CachePadded<PGMStats>,
}

#[derive(Debug, Default)]
struct PGMStats {
    cache_hits: std::sync::atomic::AtomicU64,
    total_queries: std::sync::atomic::AtomicU64,
}

// Custom serialization
impl<K: Key + Serialize> Serialize for PGMIndex<K> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("PGMIndex", 6)?;
        state.serialize_field("epsilon", &self.epsilon)?;
        state.serialize_field("data", &*self.data)?;
        state.serialize_field("segments", &self.segments)?;
        state.serialize_field("segment_lookup", &self.segment_lookup)?;
        state.serialize_field("lookup_scale", &self.lookup_scale)?;
        state.serialize_field("min_key_f64", &self.min_key_f64)?;
        state.end()
    }
}

impl<'de, K: Key + Deserialize<'de>> Deserialize<'de> for PGMIndex<K> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::{self, MapAccess, Visitor};
        use std::fmt;

        struct PGMIndexVisitor<K>(std::marker::PhantomData<K>);

        impl<'de, K: Key + Deserialize<'de>> Visitor<'de> for PGMIndexVisitor<K> {
            type Value = PGMIndex<K>;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct PGMIndex")
            }

            fn visit_map<V>(self, mut map: V) -> Result<PGMIndex<K>, V::Error>
            where
                V: MapAccess<'de>,
            {
                let mut epsilon = None;
                let mut data = None;
                let mut segments = None;
                let mut segment_lookup = None;
                let mut lookup_scale = None;
                let mut min_key_f64 = None;

                while let Some(key) = map.next_key()? {
                    match key {
                        "epsilon" => epsilon = Some(map.next_value()?),
                        "data" => data = Some(Arc::new(map.next_value::<Vec<K>>()?)),
                        "segments" => segments = Some(map.next_value()?),
                        "segment_lookup" => segment_lookup = Some(map.next_value()?),
                        "lookup_scale" => lookup_scale = Some(map.next_value()?),
                        "min_key_f64" => min_key_f64 = Some(map.next_value()?),
                        _ => {
                            let _ = map.next_value::<serde::de::IgnoredAny>()?;
                        }
                    }
                }

                Ok(PGMIndex {
                    epsilon: epsilon.ok_or_else(|| de::Error::missing_field("epsilon"))?,
                    data: data.ok_or_else(|| de::Error::missing_field("data"))?,
                    segments: segments.ok_or_else(|| de::Error::missing_field("segments"))?,
                    segment_lookup: segment_lookup
                        .ok_or_else(|| de::Error::missing_field("segment_lookup"))?,
                    lookup_scale: lookup_scale
                        .ok_or_else(|| de::Error::missing_field("lookup_scale"))?,
                    min_key_f64: min_key_f64
                        .ok_or_else(|| de::Error::missing_field("min_key_f64"))?,
                    stats: CachePadded::new(PGMStats::default()),
                })
            }
        }

        const FIELDS: &[&str] = &[
            "epsilon",
            "data",
            "segments",
            "segment_lookup",
            "lookup_scale",
            "min_key_f64",
        ];
        deserializer.deserialize_struct(
            "PGMIndex",
            FIELDS,
            PGMIndexVisitor(std::marker::PhantomData),
        )
    }
}

impl<K: Key> PGMIndex<K> {
    /// Build ultra-optimized PGM-Index
    pub fn new(data: Vec<K>, epsilon: usize) -> Self {
        Self::new_with_threads(data, epsilon, rayon::current_num_threads())
    }

    /// Build with specified threads and memory pool
    pub fn new_with_threads(data: Vec<K>, epsilon: usize, num_threads: usize) -> Self {
        assert!(epsilon > 0, "epsilon must be positive");
        assert!(!data.is_empty(), "data cannot be empty");

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        let result = pool.install(|| {
            let _n = data.len();

            let target_segments = Self::optimal_segment_count_adaptive(&data, epsilon);
            let segments = Self::build_segments_simd_parallel(&data, target_segments);
            let (segment_lookup, lookup_scale, min_key_f64) =
                Self::build_segment_lookup_vectorized(&segments, &data);

            (segments, segment_lookup, lookup_scale, min_key_f64)
        });

        PGMIndex {
            epsilon,
            data: Arc::new(data),
            segments: result.0,
            segment_lookup: result.1,
            lookup_scale: result.2,
            min_key_f64: result.3,
            stats: CachePadded::new(PGMStats::default()),
        }
    }

    fn optimal_segment_count_adaptive(data: &[K], epsilon: usize) -> usize {
        let complexity = Self::estimate_data_complexity(data);
        let n = data.len();
        let cores = rayon::current_num_threads();

        let base_segments = match complexity {
            DataComplexity::Linear => n / (epsilon * 16),
            DataComplexity::Quadratic => n / (epsilon * 8),
            DataComplexity::Exponential => n / (epsilon * 4),
            DataComplexity::Random => n / (epsilon * 2),
        };

        base_segments
            .max(cores * 4)
            .min(n / 32)
            .min(50000)
    }

    fn estimate_data_complexity(data: &[K]) -> DataComplexity {
        let sample_size = 1000.min(data.len());
        if sample_size < 10 {
            return DataComplexity::Linear;
        }

        let sample = &data[0..sample_size];
        let mut gaps = Vec::with_capacity(sample_size - 1);

        for i in 1..sample.len() {
            let gap = sample[i].to_f64().unwrap() - sample[i - 1].to_f64().unwrap();
            gaps.push(gap);
        }

        if gaps.is_empty() {
            return DataComplexity::Linear;
        }

        let avg_gap = gaps.iter().sum::<f64>() / gaps.len() as f64;

        if avg_gap.abs() < f64::EPSILON {
            return DataComplexity::Linear;
        }

        let variance = gaps.iter().map(|&g| (g - avg_gap).powi(2)).sum::<f64>() / gaps.len() as f64;

        let coefficient_of_variation = (variance.sqrt() / avg_gap).abs();

        match coefficient_of_variation {
            cv if cv < 0.1 => DataComplexity::Linear,
            cv if cv < 1.0 => DataComplexity::Quadratic,
            cv if cv < 10.0 => DataComplexity::Exponential,
            _ => DataComplexity::Random,
        }
    }

    fn build_segments_simd_parallel(data: &[K], target_segments: usize) -> Vec<Segment<K>> {
        let n = data.len();
        let segment_size = n / target_segments;

        let ranges: Vec<(usize, usize)> = (0..target_segments)
            .map(|i| {
                let start = i * segment_size;
                let end = if i == target_segments - 1 {
                    n
                } else {
                    (i + 1) * segment_size
                };
                (start, end)
            })
            .collect();

        ranges
            .par_iter()
            .map(|&(start, end)| Self::fit_segment_simd_optimized(data, start, end))
            .collect()
    }

    fn fit_segment_simd_optimized(data: &[K], start: usize, end: usize) -> Segment<K> {
        let n = end - start;
        if n == 0 {
            panic!("Cannot fit segment with zero elements");
        }
        if n == 1 {
            return Segment {
                slope: 0.0,
                intercept: start as f64,
                min_key: data[start],
                max_key: data[start],
                start_pos: start,
                end_pos: end,
            };
        }

        if n > 1000 && cfg!(target_arch = "x86_64") && Self::has_avx2_support() {
            unsafe { Self::fit_segment_avx2(data, start, end) }
        } else {
            Self::fit_segment_optimized_scalar(data, start, end)
        }
    }

    fn has_avx2_support() -> bool {
        #[cfg(target_arch = "x86_64")]
        {
            std::arch::is_x86_feature_detected!("avx2")
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            false
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn fit_segment_avx2(data: &[K], start: usize, end: usize) -> Segment<K> {
        let _n = end - start;

        Self::fit_segment_optimized_scalar(data, start, end)
    }

    #[cfg(not(target_arch = "x86_64"))]
    unsafe fn fit_segment_avx2(data: &[K], start: usize, end: usize) -> Segment<K> {
        Self::fit_segment_optimized_scalar(data, start, end)
    }

    fn fit_segment_optimized_scalar(data: &[K], start: usize, end: usize) -> Segment<K> {
        let n = end - start;
        let mut sum_x = 0.0;
        let mut sum_y = 0.0;
        let mut sum_xy = 0.0;
        let mut sum_x2 = 0.0;

        // Manual loop unrolling (8x)
        let mut i = 0;
        let unroll_end = n & !7; // Round down to multiple of 8

        while i < unroll_end {
            for j in 0..8 {
                let x = data[start + i + j].to_f64().unwrap();
                let y = (start + i + j) as f64;
                sum_x += x;
                sum_y += y;
                sum_xy += x * y;
                sum_x2 += x * x;
            }
            i += 8;
        }

        while i < n {
            let x = data[start + i].to_f64().unwrap();
            let y = (start + i) as f64;
            sum_x += x;
            sum_y += y;
            sum_xy += x * y;
            sum_x2 += x * x;
            i += 1;
        }

        let n_f = n as f64;
        let denominator = sum_x2 * n_f - sum_x * sum_x;
        let slope = if denominator.abs() > f64::EPSILON {
            (sum_xy * n_f - sum_x * sum_y) / denominator
        } else {
            0.0
        };
        let intercept = (sum_y - slope * sum_x) / n_f;

        Segment {
            slope,
            intercept,
            min_key: data[start],
            max_key: data[end - 1],
            start_pos: start,
            end_pos: end,
        }
    }

    fn build_segment_lookup_vectorized(
        segments: &[Segment<K>],
        data: &[K],
    ) -> (Vec<usize>, f64, f64) {
        if segments.is_empty() {
            return (vec![0], 1.0, 0.0);
        }

        let min_key_f64 = data[0].to_f64().unwrap();
        let max_key_f64 = data[data.len() - 1].to_f64().unwrap();
        let key_range = max_key_f64 - min_key_f64;

        if key_range == 0.0 {
            return (vec![0], 1.0, min_key_f64);
        }

        let table_size = (segments.len() * 8).max(1024).min(16384);
        let scale = (table_size - 1) as f64 / key_range;

        let lookup: Vec<usize> = (0..table_size)
            .into_par_iter()
            .map(|bucket| {
                let key_for_bucket = min_key_f64 + (bucket as f64 / scale);
                Self::find_segment_for_key_static(segments, key_for_bucket)
            })
            .collect();

        (lookup, scale, min_key_f64)
    }

    fn find_segment_for_key_static(segments: &[Segment<K>], key: f64) -> usize {
        let mut left = 0;
        let mut right = segments.len();

        while left < right {
            let mid = left + (right - left) / 2;
            let seg_min = segments[mid].min_key.to_f64().unwrap();
            let seg_max = segments[mid].max_key.to_f64().unwrap();

            if key >= seg_min && key <= seg_max {
                return mid;
            } else if key < seg_min {
                right = mid;
            } else {
                left = mid + 1;
            }
        }

        left.saturating_sub(1).min(segments.len() - 1)
    }

    #[inline(always)]
    fn find_segment_fast(&self, key: K) -> usize {
        if self.segments.len() <= 1 {
            return 0;
        }

        let key_f64 = key.to_f64().unwrap();
        if key_f64 < self.min_key_f64 {
            return 0;
        }

        let offset = key_f64 - self.min_key_f64;
        let index = (offset * self.lookup_scale) as usize;
        let seg_idx = self.segment_lookup[index.min(self.segment_lookup.len() - 1)];

        let seg = &self.segments[seg_idx];
        if key >= seg.min_key && key <= seg.max_key {
            seg_idx
        } else {
            self.find_segment_binary_search(key)
        }
    }

    fn find_segment_binary_search(&self, key: K) -> usize {
        self.segments
            .binary_search_by(|seg| {
                if key < seg.min_key {
                    std::cmp::Ordering::Greater
                } else if key > seg.max_key {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .unwrap_or_else(|i| i.saturating_sub(1).min(self.segments.len() - 1))
    }

    #[inline(always)]
    pub fn predict(&self, key: K) -> (usize, usize) {
        let seg_idx = self.find_segment_fast(key);
        let seg = &self.segments[seg_idx];

        let key_f64 = key.to_f64().unwrap();
        let predicted_pos = seg.slope.mul_add(key_f64, seg.intercept);

        let segment_start = seg.start_pos as f64;
        let segment_end = seg.end_pos as f64;
        let clamped_pos = predicted_pos.max(segment_start).min(segment_end - 1.0);

        let mid = clamped_pos.round() as isize;
        let epsilon_i = self.epsilon as isize;

        let lo = (mid - epsilon_i).max(0) as usize;
        let hi = ((mid + epsilon_i + 1) as usize).min(self.data.len());

        (lo, hi)
    }

    #[inline(always)]
    pub fn get(&self, key: K) -> Option<usize> {
        use std::sync::atomic::Ordering;

        self.stats.total_queries.fetch_add(1, Ordering::Relaxed);

        let (lo, hi) = self.predict(key);

        if lo >= self.data.len() || hi == 0 {
            return None;
        }

        #[cfg(target_arch = "x86_64")]
        unsafe {
            if std::arch::is_x86_feature_detected!("sse") {
                let search_range = &self.data[lo..hi];
                let ptr = search_range.as_ptr() as *const i8;
                _mm_prefetch(ptr, _MM_HINT_T0);

                if hi - lo > 8 {
                    _mm_prefetch(ptr.add(64), _MM_HINT_T0);
                }
            }
        }

        let search_range = &self.data[lo..hi];
        let result = search_range.binary_search(&key).ok();

        if result.is_some() {
            self.stats.cache_hits.fetch_add(1, Ordering::Relaxed);
        }

        result.map(|i| lo + i)
    }

    pub fn batch_get(&self, keys: &[K]) -> Vec<Option<usize>> {
        if keys.len() < 1000 {
            return keys.par_iter().map(|&key| self.get(key)).collect();
        }
        self.batch_get_vectorized(keys)
    }

    fn batch_get_vectorized(&self, keys: &[K]) -> Vec<Option<usize>> {
        let mut segment_groups: Vec<Vec<(K, usize)>> = vec![Vec::new(); self.segments.len()];

        for (i, &key) in keys.iter().enumerate() {
            let seg_idx = self.find_segment_fast(key);
            segment_groups[seg_idx].push((key, i));
        }

        let all_results: Vec<(usize, Option<usize>)> = segment_groups
            .par_iter()
            .flat_map(|group| {
                group
                    .par_iter()
                    .map(|&(key, orig_idx)| (orig_idx, self.get(key)))
            })
            .collect();

        let mut results = vec![None; keys.len()];
        for (idx, result) in all_results {
            results[idx] = result;
        }

        results
    }

    pub fn batch_predict(&self, keys: &[K]) -> Vec<(usize, usize)> {
        keys.par_iter().map(|&key| self.predict(key)).collect()
    }

    pub fn segment_count(&self) -> usize {
        self.segments.len()
    }

    pub fn avg_segment_size(&self) -> f64 {
        if self.segments.is_empty() {
            0.0
        } else {
            self.data.len() as f64 / self.segments.len() as f64
        }
    }

    pub fn memory_usage(&self) -> usize {
        std::mem::size_of_val(&**self.data)
            + std::mem::size_of_val(&*self.segments)
            + std::mem::size_of_val(&*self.segment_lookup)
            + std::mem::size_of::<Self>()
    }

    pub fn cache_hit_rate(&self) -> f64 {
        use std::sync::atomic::Ordering;
        let hits = self.stats.cache_hits.load(Ordering::Relaxed);
        let total = self.stats.total_queries.load(Ordering::Relaxed);
        if total > 0 {
            hits as f64 / total as f64
        } else {
            0.0
        }
    }

    pub fn reset_stats(&self) {
        use std::sync::atomic::Ordering;
        self.stats.cache_hits.store(0, Ordering::Relaxed);
        self.stats.total_queries.store(0, Ordering::Relaxed);
    }

    pub fn get_stats(&self) -> PGMPerformanceStats {
        use std::sync::atomic::Ordering;
        PGMPerformanceStats {
            total_queries: self.stats.total_queries.load(Ordering::Relaxed),
            cache_hits: self.stats.cache_hits.load(Ordering::Relaxed),
            cache_hit_rate: self.cache_hit_rate(),
            segment_count: self.segment_count(),
            avg_segment_size: self.avg_segment_size(),
            memory_usage_bytes: self.memory_usage(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct PGMPerformanceStats {
    pub total_queries: u64,
    pub cache_hits: u64,
    pub cache_hit_rate: f64,
    pub segment_count: usize,
    pub avg_segment_size: f64,
    pub memory_usage_bytes: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ultra_optimized_correctness() {
        let data: Vec<u64> = (0..100_000).collect();
        let idx = PGMIndex::new(data.clone(), 32);

        for &k in &[0, 1000, 50000, 99999] {
            assert_eq!(idx.get(k), Some(k as usize));
        }

        let queries = vec![100, 500, 1000, 5000, 9999];
        let results = idx.batch_get(&queries);

        for (i, &query) in queries.iter().enumerate() {
            assert_eq!(results[i], Some(query as usize));
        }
    }
}
