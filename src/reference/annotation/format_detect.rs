//! Format detection for GFF3 / GTF inputs. See spec §6 Stage 1.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::FerroError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotationFormat {
    Gff3,
    Gtf,
}

pub fn detect_format(path: &Path) -> Result<AnnotationFormat, FerroError> {
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        match ext {
            "gff" | "gff3" => return Ok(AnnotationFormat::Gff3),
            "gtf" => return Ok(AnnotationFormat::Gtf),
            _ => {}
        }
    }
    let file = File::open(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", path.display(), e),
    })?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("read: {}", e),
        })?;
        if line.starts_with('#') || line.trim().is_empty() {
            if line.contains("##gff-version") {
                return Ok(AnnotationFormat::Gff3);
            }
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 9 {
            return Ok(if fields[8].contains('=') && !fields[8].contains('"') {
                AnnotationFormat::Gff3
            } else {
                AnnotationFormat::Gtf
            });
        }
    }
    Ok(AnnotationFormat::Gff3)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_temp(ext: &str, content: &str) -> NamedTempFile {
        let mut tf = tempfile::Builder::new()
            .suffix(&format!(".{}", ext))
            .tempfile()
            .unwrap();
        tf.write_all(content.as_bytes()).unwrap();
        tf
    }

    #[test]
    fn detects_gff3_by_extension() {
        let f = write_temp(
            "gff3",
            "##gff-version 3\nchr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n",
        );
        let fmt = detect_format(f.path()).unwrap();
        assert_eq!(fmt, AnnotationFormat::Gff3);
    }

    #[test]
    fn detects_gtf_by_extension() {
        let f = write_temp("gtf", "chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\";\n");
        let fmt = detect_format(f.path()).unwrap();
        assert_eq!(fmt, AnnotationFormat::Gtf);
    }

    #[test]
    fn detects_gff3_by_content_when_extension_ambiguous() {
        let f = write_temp(
            "txt",
            "##gff-version 3\nchr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n",
        );
        let fmt = detect_format(f.path()).unwrap();
        assert_eq!(fmt, AnnotationFormat::Gff3);
    }

    #[test]
    fn detects_gtf_by_content_attribute_shape() {
        let f = write_temp("txt", "chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\";\n");
        let fmt = detect_format(f.path()).unwrap();
        assert_eq!(fmt, AnnotationFormat::Gtf);
    }
}
