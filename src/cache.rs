//! Caching layer for ferro-hgvs operations
//!
//! This module provides LRU caches for frequently accessed data:
//! - Parsed HGVS variants
//! - Transcript lookups
//! - Normalized variants
//!
//! # Example
//!
//! ```
//! use ferro_hgvs::cache::ParseCache;
//!
//! let cache = ParseCache::new(1000);
//! let variant = cache.get_or_parse("NM_000088.3:c.459A>G").unwrap();
//! println!("Parsed: {}", variant);
//! println!("Cache stats: {:?}", cache.stats());
//! ```

use std::collections::HashMap;
use std::hash::Hash;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::RwLock;

use crate::error::FerroError;
use crate::hgvs::parser::parse_hgvs;
use crate::hgvs::variant::HgvsVariant;

/// Statistics for cache usage
#[derive(Debug, Clone, Default)]
pub struct CacheStats {
    /// Number of cache hits
    pub hits: u64,
    /// Number of cache misses
    pub misses: u64,
    /// Number of items currently in cache
    pub size: usize,
    /// Maximum cache capacity
    pub capacity: usize,
    /// Number of evictions
    pub evictions: u64,
}

impl CacheStats {
    /// Calculate hit rate as a percentage
    pub fn hit_rate(&self) -> f64 {
        let total = self.hits + self.misses;
        if total == 0 {
            0.0
        } else {
            (self.hits as f64 / total as f64) * 100.0
        }
    }

    /// Calculate miss rate as a percentage
    pub fn miss_rate(&self) -> f64 {
        100.0 - self.hit_rate()
    }
}

/// Thread-safe LRU cache implementation
pub struct LruCache<K: Hash + Eq + Clone, V: Clone> {
    /// Cache entries
    entries: RwLock<HashMap<K, (V, u64)>>,
    /// Maximum capacity
    capacity: usize,
    /// Access counter for LRU tracking
    access_counter: AtomicU64,
    /// Hit counter
    hits: AtomicU64,
    /// Miss counter
    misses: AtomicU64,
    /// Eviction counter
    evictions: AtomicU64,
}

impl<K: Hash + Eq + Clone, V: Clone> LruCache<K, V> {
    /// Create a new LRU cache with the given capacity
    pub fn new(capacity: usize) -> Self {
        Self {
            entries: RwLock::new(HashMap::with_capacity(capacity)),
            capacity,
            access_counter: AtomicU64::new(0),
            hits: AtomicU64::new(0),
            misses: AtomicU64::new(0),
            evictions: AtomicU64::new(0),
        }
    }

    /// Get a value from the cache
    ///
    /// Note: Statistics counters use `Relaxed` ordering for performance.
    /// This means counts may be slightly inconsistent under heavy concurrent
    /// access, but this is acceptable for non-critical statistics.
    pub fn get(&self, key: &K) -> Option<V> {
        let entries = self.entries.read().unwrap();
        if let Some((value, _)) = entries.get(key) {
            self.hits.fetch_add(1, Ordering::Relaxed);
            // Note: We don't update access time on read to avoid write lock
            // This is a trade-off for better concurrent read performance
            Some(value.clone())
        } else {
            self.misses.fetch_add(1, Ordering::Relaxed);
            None
        }
    }

    /// Insert a value into the cache
    pub fn insert(&self, key: K, value: V) {
        let access = self.access_counter.fetch_add(1, Ordering::Relaxed);
        let mut entries = self.entries.write().unwrap();

        // Evict if at capacity
        if entries.len() >= self.capacity && !entries.contains_key(&key) {
            self.evict_lru(&mut entries);
        }

        entries.insert(key, (value, access));
    }

    /// Evict the least recently used entry
    fn evict_lru(&self, entries: &mut HashMap<K, (V, u64)>) {
        if let Some(lru_key) = entries
            .iter()
            .min_by_key(|(_, (_, access))| *access)
            .map(|(k, _)| k.clone())
        {
            entries.remove(&lru_key);
            self.evictions.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Clear the cache
    pub fn clear(&self) {
        let mut entries = self.entries.write().unwrap();
        entries.clear();
    }

    /// Get cache statistics
    pub fn stats(&self) -> CacheStats {
        let entries = self.entries.read().unwrap();
        CacheStats {
            hits: self.hits.load(Ordering::Relaxed),
            misses: self.misses.load(Ordering::Relaxed),
            size: entries.len(),
            capacity: self.capacity,
            evictions: self.evictions.load(Ordering::Relaxed),
        }
    }

    /// Get the number of items in the cache
    pub fn len(&self) -> usize {
        self.entries.read().unwrap().len()
    }

    /// Check if the cache is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Cache for parsed HGVS variants
pub struct ParseCache {
    cache: LruCache<String, HgvsVariant>,
}

impl ParseCache {
    /// Create a new parse cache with the given capacity
    pub fn new(capacity: usize) -> Self {
        Self {
            cache: LruCache::new(capacity),
        }
    }

    /// Get a parsed variant from cache or parse it
    pub fn get_or_parse(&self, input: &str) -> Result<HgvsVariant, FerroError> {
        if let Some(variant) = self.cache.get(&input.to_string()) {
            return Ok(variant);
        }

        let variant = parse_hgvs(input)?;
        self.cache.insert(input.to_string(), variant.clone());
        Ok(variant)
    }

    /// Get a variant from the cache only (no parsing)
    pub fn get(&self, input: &str) -> Option<HgvsVariant> {
        self.cache.get(&input.to_string())
    }

    /// Pre-populate the cache with a parsed variant
    pub fn insert(&self, input: &str, variant: HgvsVariant) {
        self.cache.insert(input.to_string(), variant);
    }

    /// Clear the cache
    pub fn clear(&self) {
        self.cache.clear();
    }

    /// Get cache statistics
    pub fn stats(&self) -> CacheStats {
        self.cache.stats()
    }

    /// Get the number of cached variants
    pub fn len(&self) -> usize {
        self.cache.len()
    }

    /// Check if the cache is empty
    pub fn is_empty(&self) -> bool {
        self.cache.is_empty()
    }
}

impl Default for ParseCache {
    fn default() -> Self {
        Self::new(1000)
    }
}

/// Cache for transcript lookups
pub struct TranscriptCache<T: Clone> {
    cache: LruCache<String, T>,
}

impl<T: Clone> TranscriptCache<T> {
    /// Create a new transcript cache with the given capacity
    pub fn new(capacity: usize) -> Self {
        Self {
            cache: LruCache::new(capacity),
        }
    }

    /// Get a transcript from the cache
    pub fn get(&self, accession: &str) -> Option<T> {
        self.cache.get(&accession.to_string())
    }

    /// Insert a transcript into the cache
    pub fn insert(&self, accession: &str, transcript: T) {
        self.cache.insert(accession.to_string(), transcript);
    }

    /// Get or compute a transcript
    pub fn get_or_insert_with<F>(&self, accession: &str, f: F) -> Option<T>
    where
        F: FnOnce() -> Option<T>,
    {
        if let Some(transcript) = self.get(accession) {
            return Some(transcript);
        }

        if let Some(transcript) = f() {
            self.insert(accession, transcript.clone());
            Some(transcript)
        } else {
            None
        }
    }

    /// Clear the cache
    pub fn clear(&self) {
        self.cache.clear();
    }

    /// Get cache statistics
    pub fn stats(&self) -> CacheStats {
        self.cache.stats()
    }
}

impl<T: Clone> Default for TranscriptCache<T> {
    fn default() -> Self {
        Self::new(1000)
    }
}

/// Combined cache for all operations
pub struct OperationCache {
    /// Parse cache
    pub parse: ParseCache,
    /// Normalized variant cache (keyed by original string)
    pub normalized: LruCache<String, HgvsVariant>,
}

impl OperationCache {
    /// Create a new operation cache with the given capacities
    pub fn new(parse_capacity: usize, normalize_capacity: usize) -> Self {
        Self {
            parse: ParseCache::new(parse_capacity),
            normalized: LruCache::new(normalize_capacity),
        }
    }

    /// Clear all caches
    pub fn clear_all(&self) {
        self.parse.clear();
        self.normalized.clear();
    }

    /// Get combined statistics
    pub fn all_stats(&self) -> HashMap<&'static str, CacheStats> {
        let mut stats = HashMap::new();
        stats.insert("parse", self.parse.stats());
        stats.insert("normalized", self.normalized.stats());
        stats
    }
}

impl Default for OperationCache {
    fn default() -> Self {
        Self::new(1000, 1000)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lru_cache_basic() {
        let cache: LruCache<String, i32> = LruCache::new(3);

        cache.insert("a".to_string(), 1);
        cache.insert("b".to_string(), 2);
        cache.insert("c".to_string(), 3);

        assert_eq!(cache.get(&"a".to_string()), Some(1));
        assert_eq!(cache.get(&"b".to_string()), Some(2));
        assert_eq!(cache.get(&"c".to_string()), Some(3));
        assert_eq!(cache.get(&"d".to_string()), None);
    }

    #[test]
    fn test_lru_cache_eviction() {
        let cache: LruCache<String, i32> = LruCache::new(2);

        cache.insert("a".to_string(), 1);
        cache.insert("b".to_string(), 2);
        cache.insert("c".to_string(), 3); // Should evict "a"

        assert_eq!(cache.len(), 2);
        assert!(cache.get(&"a".to_string()).is_none()); // Evicted
        assert!(cache.get(&"b".to_string()).is_some());
        assert!(cache.get(&"c".to_string()).is_some());
    }

    #[test]
    fn test_cache_stats() {
        let cache: LruCache<String, i32> = LruCache::new(10);

        cache.insert("a".to_string(), 1);
        cache.get(&"a".to_string()); // Hit
        cache.get(&"b".to_string()); // Miss

        let stats = cache.stats();
        assert_eq!(stats.hits, 1);
        assert_eq!(stats.misses, 1);
        assert_eq!(stats.size, 1);
        assert!((stats.hit_rate() - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_parse_cache() {
        let cache = ParseCache::new(100);

        // First parse should miss and parse
        let v1 = cache.get_or_parse("NM_000088.3:c.459A>G").unwrap();

        // Second should hit cache
        let v2 = cache.get_or_parse("NM_000088.3:c.459A>G").unwrap();

        assert_eq!(format!("{}", v1), format!("{}", v2));

        let stats = cache.stats();
        assert_eq!(stats.hits, 1);
        assert_eq!(stats.misses, 1);
    }

    #[test]
    fn test_transcript_cache() {
        let cache: TranscriptCache<String> = TranscriptCache::new(100);

        cache.insert("NM_000088.3", "test transcript".to_string());
        assert_eq!(
            cache.get("NM_000088.3"),
            Some("test transcript".to_string())
        );
        assert_eq!(cache.get("NM_000099.1"), None);
    }

    #[test]
    fn test_operation_cache() {
        let cache = OperationCache::default();

        let _ = cache.parse.get_or_parse("NM_000088.3:c.459A>G");
        cache.normalized.insert(
            "test".to_string(),
            parse_hgvs("NC_000001.11:g.12345A>G").unwrap(),
        );

        let all_stats = cache.all_stats();
        assert!(all_stats.contains_key("parse"));
        assert!(all_stats.contains_key("normalized"));
    }

    #[test]
    fn test_cache_stats_miss_rate() {
        let stats = CacheStats {
            hits: 75,
            misses: 25,
            size: 100,
            capacity: 1000,
            evictions: 0,
        };
        assert!((stats.hit_rate() - 75.0).abs() < 0.01);
        assert!((stats.miss_rate() - 25.0).abs() < 0.01);
    }

    #[test]
    fn test_cache_stats_zero_total() {
        let stats = CacheStats::default();
        assert!((stats.hit_rate() - 0.0).abs() < 0.01);
        assert!((stats.miss_rate() - 100.0).abs() < 0.01);
    }

    #[test]
    fn test_lru_cache_clear() {
        let cache: LruCache<String, i32> = LruCache::new(10);
        cache.insert("a".to_string(), 1);
        cache.insert("b".to_string(), 2);
        cache.insert("c".to_string(), 3);

        assert_eq!(cache.len(), 3);
        assert!(!cache.is_empty());

        cache.clear();

        assert_eq!(cache.len(), 0);
        assert!(cache.is_empty());
        assert!(cache.get(&"a".to_string()).is_none());
    }

    #[test]
    fn test_lru_cache_is_empty() {
        let cache: LruCache<String, i32> = LruCache::new(10);
        assert!(cache.is_empty());

        cache.insert("a".to_string(), 1);
        assert!(!cache.is_empty());
    }

    #[test]
    fn test_lru_cache_update_existing() {
        let cache: LruCache<String, i32> = LruCache::new(3);

        cache.insert("a".to_string(), 1);
        cache.insert("b".to_string(), 2);
        cache.insert("c".to_string(), 3);

        // Update existing key should not trigger eviction
        cache.insert("a".to_string(), 10);

        assert_eq!(cache.len(), 3);
        assert_eq!(cache.get(&"a".to_string()), Some(10));
        assert_eq!(cache.get(&"b".to_string()), Some(2));
        assert_eq!(cache.get(&"c".to_string()), Some(3));
    }

    #[test]
    fn test_lru_cache_stats_eviction_count() {
        let cache: LruCache<String, i32> = LruCache::new(2);

        cache.insert("a".to_string(), 1);
        cache.insert("b".to_string(), 2);
        cache.insert("c".to_string(), 3); // Evicts "a"
        cache.insert("d".to_string(), 4); // Evicts "b"

        let stats = cache.stats();
        assert_eq!(stats.evictions, 2);
        assert_eq!(stats.size, 2);
        assert_eq!(stats.capacity, 2);
    }

    #[test]
    fn test_parse_cache_get_only() {
        let cache = ParseCache::new(100);

        // Get should return None when nothing is cached
        assert!(cache.get("NM_000088.3:c.459A>G").is_none());

        // Parse and cache it
        let _ = cache.get_or_parse("NM_000088.3:c.459A>G").unwrap();

        // Now get should return the cached value
        assert!(cache.get("NM_000088.3:c.459A>G").is_some());
    }

    #[test]
    fn test_parse_cache_insert() {
        let cache = ParseCache::new(100);

        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        cache.insert("NC_000001.11:g.12345A>G", variant.clone());

        let retrieved = cache.get("NC_000001.11:g.12345A>G");
        assert!(retrieved.is_some());
        assert_eq!(format!("{}", retrieved.unwrap()), format!("{}", variant));
    }

    #[test]
    fn test_parse_cache_clear() {
        let cache = ParseCache::new(100);

        let _ = cache.get_or_parse("NM_000088.3:c.459A>G").unwrap();
        assert!(!cache.is_empty());

        cache.clear();
        assert!(cache.is_empty());
        assert!(cache.get("NM_000088.3:c.459A>G").is_none());
    }

    #[test]
    fn test_parse_cache_len_and_is_empty() {
        let cache = ParseCache::new(100);

        assert!(cache.is_empty());
        assert_eq!(cache.len(), 0);

        let _ = cache.get_or_parse("NM_000088.3:c.459A>G").unwrap();
        assert!(!cache.is_empty());
        assert_eq!(cache.len(), 1);

        let _ = cache.get_or_parse("NC_000001.11:g.12345A>G").unwrap();
        assert_eq!(cache.len(), 2);
    }

    #[test]
    fn test_parse_cache_default() {
        let cache = ParseCache::default();
        // Default should have capacity of 1000
        let stats = cache.stats();
        assert_eq!(stats.capacity, 1000);
    }

    #[test]
    fn test_transcript_cache_get_or_insert_with() {
        let cache: TranscriptCache<String> = TranscriptCache::new(100);

        // Should compute and cache
        let result =
            cache.get_or_insert_with("NM_000088.3", || Some("computed_transcript".to_string()));
        assert_eq!(result, Some("computed_transcript".to_string()));

        // Should return cached value (closure not called)
        let result = cache.get_or_insert_with("NM_000088.3", || {
            panic!("Closure should not be called for cached value");
        });
        assert_eq!(result, Some("computed_transcript".to_string()));
    }

    #[test]
    fn test_transcript_cache_get_or_insert_with_none() {
        let cache: TranscriptCache<String> = TranscriptCache::new(100);

        // Closure returns None
        let result = cache.get_or_insert_with("NM_000088.3", || None);
        assert!(result.is_none());

        // Should still be empty (None is not cached)
        let stats = cache.stats();
        assert_eq!(stats.size, 0);
    }

    #[test]
    fn test_transcript_cache_clear() {
        let cache: TranscriptCache<String> = TranscriptCache::new(100);

        cache.insert("NM_000088.3", "test1".to_string());
        cache.insert("NM_000099.1", "test2".to_string());

        let stats_before = cache.stats();
        assert_eq!(stats_before.size, 2);

        cache.clear();

        let stats_after = cache.stats();
        assert_eq!(stats_after.size, 0);
    }

    #[test]
    fn test_transcript_cache_default() {
        let cache: TranscriptCache<String> = TranscriptCache::default();
        let stats = cache.stats();
        assert_eq!(stats.capacity, 1000);
    }

    #[test]
    fn test_operation_cache_clear_all() {
        let cache = OperationCache::new(100, 100);

        // Add some items
        let _ = cache.parse.get_or_parse("NM_000088.3:c.459A>G");
        cache.normalized.insert(
            "test".to_string(),
            parse_hgvs("NC_000001.11:g.12345A>G").unwrap(),
        );

        assert!(!cache.parse.is_empty());
        assert!(!cache.normalized.is_empty());

        cache.clear_all();

        assert!(cache.parse.is_empty());
        assert!(cache.normalized.is_empty());
    }

    #[test]
    fn test_operation_cache_custom_capacities() {
        let cache = OperationCache::new(50, 75);

        let parse_stats = cache.parse.stats();
        let normalized_stats = cache.normalized.stats();

        assert_eq!(parse_stats.capacity, 50);
        assert_eq!(normalized_stats.capacity, 75);
    }

    #[test]
    fn test_parse_cache_invalid_input() {
        let cache = ParseCache::new(100);

        // Invalid HGVS should return an error
        let result = cache.get_or_parse("invalid-hgvs-string");
        assert!(result.is_err());

        // Invalid inputs should not be cached
        assert!(cache.is_empty());
    }
}
