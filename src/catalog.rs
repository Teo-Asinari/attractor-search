// Save/load discovered attractors.

use crate::ode::{Coeffs, State};
use serde::{Deserialize, Serialize};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::path::Path;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Entry {
    pub id: u64,
    pub coeffs: Vec<f64>,
    pub spectrum: [f64; 3],
    pub ky_dim: f64,
    pub trajectory: Vec<[f64; 3]>,
    #[serde(default)]
    pub method: String,
}

impl Entry {
    pub fn new(
        coeffs: &Coeffs,
        spectrum: [f64; 3],
        ky_dim: f64,
        traj: &[State],
        method: &str,
    ) -> Self {
        let id = coeff_hash(coeffs);
        Entry {
            id,
            coeffs: coeffs.to_vec(),
            spectrum,
            ky_dim,
            trajectory: traj.to_vec(),
            method: method.to_string(),
        }
    }
}

/// Hash coefficients for quick ID.
pub fn coeff_hash(c: &Coeffs) -> u64 {
    let mut h = DefaultHasher::new();
    for &v in c.iter() {
        v.to_bits().hash(&mut h);
    }
    h.finish()
}

/// Save entry to JSON file in dir.
pub fn save(
    dir: &Path,
    entry: &Entry,
) -> std::io::Result<()> {
    std::fs::create_dir_all(dir)?;
    let path = dir.join(
        format!("{:016x}.json", entry.id),
    );
    let json = serde_json::to_string_pretty(entry)
        .map_err(|e| {
            std::io::Error::new(
                std::io::ErrorKind::Other,
                e,
            )
        })?;
    std::fs::write(path, json)
}

/// Load all entries from dir.
pub fn load_all(
    dir: &Path,
) -> std::io::Result<Vec<Entry>> {
    let mut entries = Vec::new();
    if !dir.exists() {
        return Ok(entries);
    }
    for item in std::fs::read_dir(dir)? {
        let item = item?;
        let p = item.path();
        if p.extension().map_or(false, |e| e == "json")
        {
            let data = std::fs::read_to_string(&p)?;
            if let Ok(e) =
                serde_json::from_str::<Entry>(&data)
            {
                entries.push(e);
            }
        }
    }
    Ok(entries)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ode::NCOEFFS;
    use std::path::PathBuf;

    #[test]
    fn roundtrip() {
        let c = [0.5; NCOEFFS];
        let traj = vec![[1.0, 2.0, 3.0]; 5];
        let entry = Entry::new(
            &c,
            [0.9, 0.0, -14.0],
            2.06,
            &traj,
            "test",
        );
        let dir = PathBuf::from("/tmp/attractor_test");
        let _ = std::fs::remove_dir_all(&dir);
        save(&dir, &entry).unwrap();
        let loaded = load_all(&dir).unwrap();
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded[0].id, entry.id);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
