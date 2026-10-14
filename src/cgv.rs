//! The cgv program (cgv/ in this repository: exact GV invariants in C), built by build.rs and
//! bundled with the crate. [`executable`] writes it to a cache directory on first use and returns
//! its path; callers run it like the standalone program (see cgv/README.md).

use std::fs;
use std::io;
use std::path::PathBuf;

#[cfg(not(windows))]
static CGV: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/cgv"));
#[cfg(windows)]
static CGV: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/cgv.exe"));

/// Path to the bundled cgv program, written to the cache directory on first use.
pub fn executable() -> io::Result<PathBuf> {
    // FNV-1a of the program: a new build never reuses an old file
    let hash = CGV.iter().fold(0xcbf29ce484222325u64, |h, &b| {
        (h ^ b as u64).wrapping_mul(0x100000001b3)
    });
    let dir = cache_dir().join("cygv");
    let name = format!(
        "cgv-{}-{hash:016x}{}",
        env!("CARGO_PKG_VERSION"),
        std::env::consts::EXE_SUFFIX
    );
    let path = dir.join(&name);
    if !path.exists() {
        fs::create_dir_all(&dir)?;
        let tmp = dir.join(format!("{name}.{}.tmp", std::process::id()));
        fs::write(&tmp, CGV)?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            fs::set_permissions(&tmp, fs::Permissions::from_mode(0o755))?;
        }
        fs::rename(&tmp, &path)?; // atomic, so concurrent first uses are fine
    }
    Ok(path)
}

fn cache_dir() -> PathBuf {
    let var = |k: &str| {
        std::env::var_os(k)
            .filter(|v| !v.is_empty())
            .map(PathBuf::from)
    };
    if cfg!(windows) {
        var("LOCALAPPDATA").unwrap_or_else(std::env::temp_dir)
    } else {
        var("XDG_CACHE_HOME")
            .or_else(|| var("HOME").map(|h| h.join(".cache")))
            .unwrap_or_else(std::env::temp_dir)
    }
}
