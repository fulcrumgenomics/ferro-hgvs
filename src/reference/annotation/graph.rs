//! Two-pass feature graph: ingest then resolve. See spec §6 Stage 3.

use std::collections::HashMap;

use smol_str::SmolStr;

use super::feature::{Feature, FeatureType};
use super::record::AnnotationRecord;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FeatureId(pub u32);

#[derive(Debug, Default)]
pub struct FeatureGraph {
    pub(crate) features: Vec<Feature>,
    pub(crate) by_id: HashMap<SmolStr, FeatureId>,
    pub(crate) children: HashMap<FeatureId, Vec<FeatureId>>,
    pub(crate) roots: Vec<FeatureId>,
    pub(crate) orphan_count: usize,
    pub(crate) orphan_diagnostics: Vec<(u64, String, Option<String>)>, // (line, parent_id, child_id)
}

impl FeatureGraph {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn ingest<R: AnnotationRecord>(&mut self, rec: R) {
        let fid = FeatureId(
            u32::try_from(self.features.len()).expect("annotation file exceeds u32::MAX features"),
        );
        let id = rec.id().map(SmolStr::new);
        let parent_ids: Vec<SmolStr> = rec.parents().iter().map(SmolStr::new).collect();
        let feature_type = rec.feature_type().clone();
        let seqid = SmolStr::new(rec.seqid());
        let range = (rec.start(), rec.end());
        let strand = rec.strand();
        let phase = rec.phase();
        let source_line = rec.source_line();
        let attrs = rec.into_attrs();
        let feature = Feature {
            id: id.clone(),
            parent_ids,
            feature_type,
            seqid,
            range,
            strand,
            phase,
            attrs,
            source_line,
        };
        self.features.push(feature);
        if let Some(id) = id {
            self.by_id.insert(id, fid);
        }
    }

    pub fn resolve(&mut self) -> usize {
        // Derived state is rebuilt from scratch; calling resolve() twice on
        // the same graph must not duplicate edges or diagnostics.
        self.children.clear();
        self.roots.clear();
        self.orphan_diagnostics.clear();
        self.orphan_count = 0;

        let mut orphan_count = 0;
        for i in 0..self.features.len() {
            // Safe by invariant: `ingest` panics if a feature push would exceed
            // u32::MAX, so here `i < features.len() <= u32::MAX`.
            let fid = FeatureId(i as u32);
            let parent_ids: Vec<SmolStr> = self.features[i].parent_ids.clone();
            if parent_ids.is_empty() {
                self.roots.push(fid);
                continue;
            }
            let mut any_resolved = false;
            for pid in &parent_ids {
                if let Some(&parent_fid) = self.by_id.get(pid) {
                    self.children.entry(parent_fid).or_default().push(fid);
                    any_resolved = true;
                }
            }
            if !any_resolved {
                // A node whose parents are all missing is an orphan only when it
                // has no ID of its own.  If it has an ID, other features may
                // reference it as a parent, so we keep it as a de-facto root.
                if self.features[i].id.is_none() {
                    orphan_count += 1;
                    for pid in &parent_ids {
                        self.orphan_diagnostics.push((
                            self.features[i].source_line,
                            pid.to_string(),
                            None,
                        ));
                    }
                } else {
                    self.roots.push(fid);
                }
            }
        }
        self.orphan_count = orphan_count;
        orphan_count
    }

    /// Look up a feature by its string ID.
    ///
    /// Named `feature_by_id` (not `by_id`) to avoid shadowing the `by_id` field.
    /// Used in tests and Phase 2+ builder code.
    #[allow(dead_code)]
    pub fn feature_by_id(&self, id: &str) -> Option<FeatureId> {
        self.by_id.get(id).copied()
    }

    pub fn children_of(&self, fid: FeatureId) -> &[FeatureId] {
        self.children.get(&fid).map(|v| v.as_slice()).unwrap_or(&[])
    }

    pub fn feature(&self, fid: FeatureId) -> &Feature {
        &self.features[fid.0 as usize]
    }

    pub fn descendants_of_type(&self, fid: FeatureId, ty: &FeatureType) -> Vec<FeatureId> {
        let mut out = Vec::new();
        let mut stack = vec![fid];
        while let Some(node) = stack.pop() {
            for &c in self.children_of(node) {
                if &self.features[c.0 as usize].feature_type == ty {
                    out.push(c);
                }
                stack.push(c);
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::annotation::record::Gff3Record;

    fn add_line(graph: &mut FeatureGraph, line: &str, n: u64) {
        let rec = Gff3Record::parse(line, n).unwrap().unwrap();
        graph.ingest(rec);
    }

    #[test]
    fn graph_resolves_parent_child_basic() {
        let mut g = FeatureGraph::new();
        add_line(
            &mut g,
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=g1",
            1,
        );
        add_line(&mut g, "chr1\t.\texon\t100\t200\t.\t+\t.\tParent=tx1", 2);
        add_line(&mut g, "chr1\t.\texon\t300\t500\t.\t+\t.\tParent=tx1", 3);
        let dropped = g.resolve();
        assert_eq!(dropped, 0);
        let tx_id = g.feature_by_id("tx1").unwrap();
        let children = g.children_of(tx_id);
        assert_eq!(children.len(), 2);
    }

    #[test]
    fn graph_handles_forward_reference() {
        let mut g = FeatureGraph::new();
        add_line(&mut g, "chr1\t.\texon\t100\t200\t.\t+\t.\tParent=tx1", 1);
        add_line(&mut g, "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1", 2);
        let dropped = g.resolve();
        assert_eq!(dropped, 0);
        let tx_id = g.feature_by_id("tx1").unwrap();
        assert_eq!(g.children_of(tx_id).len(), 1);
    }

    #[test]
    fn graph_drops_orphans() {
        let mut g = FeatureGraph::new();
        add_line(
            &mut g,
            "chr1\t.\texon\t100\t200\t.\t+\t.\tParent=missing",
            1,
        );
        let dropped = g.resolve();
        assert_eq!(dropped, 1);
    }

    #[test]
    fn resolve_is_idempotent() {
        let mut g = FeatureGraph::new();
        add_line(&mut g, "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1", 1);
        add_line(&mut g, "chr1\t.\texon\t100\t200\t.\t+\t.\tParent=tx1", 2);
        add_line(&mut g, "chr1\t.\texon\t300\t500\t.\t+\t.\tParent=tx1", 3);
        add_line(
            &mut g,
            "chr1\t.\texon\t100\t200\t.\t+\t.\tParent=missing",
            4,
        );

        let first = g.resolve();
        let tx = g.feature_by_id("tx1").unwrap();
        let children_first = g.children_of(tx).len();
        let roots_first = g.roots.len();
        let orphans_first = g.orphan_diagnostics.len();

        let second = g.resolve();
        assert_eq!(first, second, "orphan count must be stable across calls");
        assert_eq!(g.children_of(tx).len(), children_first);
        assert_eq!(g.roots.len(), roots_first);
        assert_eq!(g.orphan_diagnostics.len(), orphans_first);
    }

    #[test]
    fn graph_supports_multi_parent_exon() {
        let mut g = FeatureGraph::new();
        add_line(&mut g, "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1", 1);
        add_line(&mut g, "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx2", 2);
        add_line(
            &mut g,
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1,tx2",
            3,
        );
        let dropped = g.resolve();
        assert_eq!(dropped, 0);
        let tx1 = g.feature_by_id("tx1").unwrap();
        let tx2 = g.feature_by_id("tx2").unwrap();
        assert_eq!(g.children_of(tx1).len(), 1);
        assert_eq!(g.children_of(tx2).len(), 1);
    }
}
