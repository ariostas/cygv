//! Builds the cgv program (cgv/, a standalone C program) when the `cgv` feature is on.
//!
//! Same as the `cgv` target of cgv/Makefile: gv.c is compiled once per number of prime lanes
//! (2, 3, 4) and linked with main.c. The only difference is portable instead of native CPU
//! flags, so the result runs on any machine of the target architecture (~2% slower).
//!
//! On Linux it also builds the GPU variant (the `cgv_gpu` and `cgv_hip` targets of the Makefile)
//! if a GPU compiler is found: CUDA's nvcc (`NVCC`, `CUDA_HOME`/`CUDA_PATH`, or `PATH`) or HIP's
//! hipcc (`HIPCC`, `ROCM_PATH`, or `PATH`). GPU architectures: `CGV_CUDA_ARCH` /
//! `CGV_HIP_ARCH` (comma-separated, or `native`; default `native` for the GPUs present, and for
//! CUDA without a GPU every architecture from sm_75 that nvcc supports). With both compilers present,
//! the GPU in the machine decides. `CGV_GPU=0` skips the GPU variant; `CGV_GPU=cuda` / `hip` picks
//! the toolkit, and `CGV_GPU=1` / `cuda` / `hip` make a failure to build it an error, not a warning.

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    for f in [
        "cgv/main.c",
        "cgv/gv.c",
        "cgv/cgv.h",
        "cgv/gpu.cu",
        "cgv/gpu.h",
        "cgv/gpu_compat.h",
    ] {
        println!("cargo:rerun-if-changed={f}");
    }
    for v in [
        "CGV_GPU",
        "CGV_CUDA_ARCH",
        "CGV_HIP_ARCH",
        "NVCC",
        "HIPCC",
        "CUDA_HOME",
        "CUDA_PATH",
        "ROCM_PATH",
        "PATH",
    ] {
        println!("cargo:rerun-if-env-changed={v}");
    }
    println!("cargo::rustc-check-cfg=cfg(cgv_gpu, cgv_hip)");
    if env::var_os("CARGO_FEATURE_CGV").is_none() {
        return;
    }
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let windows = env::var("CARGO_CFG_TARGET_OS").unwrap() == "windows";
    let linux = env::var("CARGO_CFG_TARGET_OS").unwrap() == "linux";
    let x86_64 = env::var("CARGO_CFG_TARGET_ARCH").unwrap() == "x86_64";
    let compiler = cc::Build::new().get_compiler();
    if compiler.is_like_msvc() {
        panic!("the cgv feature needs gcc or clang (cgv uses 128-bit integers, which MSVC lacks); build with MinGW, or disable the feature");
    }
    let flags: &[&str] = if x86_64 {
        &["-O3", "-march=x86-64-v2"]
    } else {
        &["-O3"]
    };
    let cc = |args: &[String]| {
        let mut cmd = compiler.to_command();
        cmd.args(flags).args(args);
        run(cmd)
    };
    let s = |p: &Path| p.display().to_string();

    // CPU program
    let mut objects = vec![];
    for nl in [2, 3, 4] {
        let obj = out.join(format!("gv_nl{nl}.o"));
        cc(&[
            format!("-DNL={nl}"),
            format!("-DCGV_ENTRY=cgv_entry_nl{nl}"),
            "-c".into(),
            "cgv/gv.c".into(),
            "-o".into(),
            s(&obj),
        ])
        .expect("building cgv failed");
        objects.push(s(&obj));
    }
    let main_obj = s(&out.join("main.o"));
    cc(&[
        "-c".into(),
        "cgv/main.c".into(),
        "-o".into(),
        main_obj.clone(),
    ])
    .expect("building cgv failed");
    let exe = s(&out.join(if windows { "cgv.exe" } else { "cgv" }));
    let mut link = vec!["-o".to_string(), exe, main_obj.clone()];
    link.extend(objects);
    link.extend(["-lpthread".to_string(), "-lm".to_string()]);
    cc(&link).expect("building cgv failed");

    // GPU program (optional)
    let want = env::var("CGV_GPU").unwrap_or_default();
    if want == "0" || !linux {
        return;
    }
    // follow the hardware: CUDA for a visible NVIDIA GPU, HIP for an AMD GPU (/dev/kfd); without a GPU,
    // CUDA if nvcc is there (all architectures), HIP only if CGV_HIP_ARCH names them
    let nvcc = find_tool("NVCC", "nvcc", &["CUDA_HOME", "CUDA_PATH"]);
    let hipcc = find_tool("HIPCC", "hipcc", &["ROCM_PATH"]);
    let nvidia = nvidia_gpu_visible();
    let amd = Path::new("/dev/kfd").exists();
    let hip_arch_set = env::var_os("CGV_HIP_ARCH").is_some();
    let choice = match (want.as_str(), &nvcc, &hipcc) {
        ("cuda", Some(n), _) => Some((n.clone(), false)),
        ("hip", _, Some(h)) => Some((h.clone(), true)),
        ("cuda" | "hip", _, _) => None,
        (_, Some(n), _) if nvidia => Some((n.clone(), false)),
        (_, _, Some(h)) if amd => Some((h.clone(), true)),
        (_, Some(n), _) => Some((n.clone(), false)),
        (_, _, Some(h)) if hip_arch_set => Some((h.clone(), true)),
        _ => None,
    };
    let result = match choice {
        Some((gpucc, hip)) => build_gpu(&out, &gpucc, hip, &cc, &main_obj),
        None => Err("no GPU compiler (nvcc or hipcc) found".to_string()),
    };
    match result {
        Ok(hip) => {
            println!("cargo:rustc-cfg=cgv_gpu");
            if hip {
                println!("cargo:rustc-cfg=cgv_hip");
            }
        }
        Err(e) if !want.is_empty() => {
            panic!("CGV_GPU={want}, but the GPU variant of cgv was not built: {e}")
        }
        Err(e) if e.starts_with("no GPU compiler") => {}
        Err(e) => {
            println!("cargo:warning=cgv: GPU variant not built ({e}); the CPU program is built")
        }
    }
}

/// Whether nvidia-smi lists a GPU.
fn nvidia_gpu_visible() -> bool {
    Command::new("nvidia-smi")
        .arg("-L")
        .output()
        .is_ok_and(|o| o.status.success() && !o.stdout.is_empty())
}

/// Runs a command; Err with the command on failure.
fn run(mut cmd: Command) -> Result<(), String> {
    match cmd.status() {
        Ok(st) if st.success() => Ok(()),
        Ok(_) => Err(format!("{cmd:?} failed")),
        Err(e) => Err(format!("cannot run {cmd:?}: {e}")),
    }
}

/// A tool from an explicit variable (e.g. NVCC), a toolkit root (e.g. CUDA_HOME/bin), or PATH.
fn find_tool(var: &str, name: &str, roots: &[&str]) -> Option<PathBuf> {
    if let Some(p) = env::var_os(var) {
        return Some(PathBuf::from(p));
    }
    let dirs = roots
        .iter()
        .filter_map(|r| env::var_os(r).map(|d| PathBuf::from(d).join("bin")))
        .chain(
            env::var_os("PATH")
                .map(|p| env::split_paths(&p).collect::<Vec<_>>())
                .unwrap_or_default(),
        );
    dirs.map(|d| d.join(name)).find(|p| p.is_file())
}

/// Builds cgv_gpu like the Makefile's cgv_gpu (CUDA) or cgv_hip (HIP) target. Ok(hip) on success.
fn build_gpu(
    out: &Path,
    gpucc: &Path,
    hip: bool,
    cc: &dyn Fn(&[String]) -> Result<(), String>,
    main_obj: &str,
) -> Result<bool, String> {
    let arch = arch_flags(gpucc, hip)?;
    let mut objects = vec![main_obj.to_string()];
    for nl in [2, 3, 4] {
        let renames = [
            format!("-DNL={nl}"),
            format!("-Dgpu_extract=gpu_extract_nl{nl}"),
            format!("-Dgpu_free_gb=gpu_free_gb_nl{nl}"),
        ];
        let obj = out.join(format!("gvg_nl{nl}.o")).display().to_string();
        let mut args = vec![
            "-DUSE_GPU".to_string(),
            format!("-DCGV_ENTRY=cgv_entry_nl{nl}"),
        ];
        args.extend(renames.iter().cloned());
        args.extend(["-c".into(), "cgv/gv.c".into(), "-o".into(), obj.clone()]);
        cc(&args)?;
        objects.push(obj);
        let obj = out.join(format!("gpu_nl{nl}.o")).display().to_string();
        let mut cmd = Command::new(gpucc);
        cmd.arg("-O3");
        if hip {
            cmd.args(["-x", "hip"]);
        }
        cmd.args(&arch)
            .args(&renames)
            .args(["-c", "cgv/gpu.cu", "-o", &obj]);
        run(cmd)?;
        objects.push(obj);
    }
    let exe = out.join("cgv_gpu").display().to_string();
    let mut cmd = Command::new(gpucc);
    cmd.args(&arch)
        .arg("-o")
        .arg(&exe)
        .args(&objects)
        .args(["-lpthread", "-lm"]);
    run(cmd)?;
    Ok(hip)
}

/// GPU architecture flags: CGV_CUDA_ARCH / CGV_HIP_ARCH, else `native` (CUDA without a visible GPU:
/// every architecture from sm_75 that this nvcc supports, plus PTX for the newest).
fn arch_flags(gpucc: &Path, hip: bool) -> Result<Vec<String>, String> {
    let var = if hip { "CGV_HIP_ARCH" } else { "CGV_CUDA_ARCH" };
    let mut archs: Vec<String> = env::var(var)
        .unwrap_or_default()
        .split(',')
        .map(|a| a.trim().to_string())
        .filter(|a| !a.is_empty())
        .collect();
    if hip {
        if archs.is_empty() {
            archs.push("native".into());
        }
        return Ok(archs
            .iter()
            .map(|a| format!("--offload-arch={a}"))
            .collect());
    }
    if archs.is_empty() {
        if nvidia_gpu_visible() {
            archs.push("native".into());
        } else {
            let list = Command::new(gpucc)
                .arg("--list-gpu-arch")
                .output()
                .map_err(|e| format!("cannot run nvcc: {e}"))?;
            archs = String::from_utf8_lossy(&list.stdout)
                .lines()
                .filter_map(|l| {
                    l.trim()
                        .strip_prefix("compute_")
                        .and_then(|n| n.parse::<u32>().ok())
                })
                .filter(|&n| n >= 75)
                .map(|n| format!("sm_{n}"))
                .collect();
            if archs.is_empty() {
                return Err("nvcc lists no GPU architecture from sm_75 on".into());
            }
        }
    }
    if archs == ["native"] {
        return Ok(vec!["-arch=native".into()]);
    }
    let mut flags: Vec<String> = archs
        .iter()
        .map(|a| {
            let n = a.trim_start_matches("sm_");
            format!("-gencode=arch=compute_{n},code=sm_{n}")
        })
        .collect();
    // PTX for the newest architecture, so newer GPUs can still run it (JIT)
    if let Some(n) = archs
        .iter()
        .filter_map(|a| a.trim_start_matches("sm_").parse::<u32>().ok())
        .max()
    {
        flags.push(format!("-gencode=arch=compute_{n},code=compute_{n}"));
    }
    Ok(flags)
}
