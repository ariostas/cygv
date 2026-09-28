//! Builds the cgv program (cgv/, a standalone C program) when the `cgv` feature is on.
//!
//! Same as the `cgv` target of cgv/Makefile: gv.c is compiled once per number of prime lanes
//! (2, 3, 4) and linked with main.c. The only difference is portable instead of native CPU
//! flags, so the result runs on any machine of the target architecture (~2% slower).

use std::env;
use std::path::PathBuf;

fn main() {
    for f in ["cgv/main.c", "cgv/gv.c", "cgv/cgv.h"] {
        println!("cargo:rerun-if-changed={f}");
    }
    if env::var_os("CARGO_FEATURE_CGV").is_none() {
        return;
    }
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let windows = env::var("CARGO_CFG_TARGET_OS").unwrap() == "windows";
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
    let run = |args: &[String]| {
        let mut cmd = compiler.to_command();
        cmd.args(flags).args(args);
        let status = cmd
            .status()
            .unwrap_or_else(|e| panic!("cannot run the C compiler: {e}"));
        assert!(status.success(), "building cgv failed: {cmd:?}");
    };
    let mut objects = vec![];
    for nl in [2, 3, 4] {
        let obj = out.join(format!("gv_nl{nl}.o"));
        run(&[
            format!("-DNL={nl}"),
            format!("-DCGV_ENTRY=cgv_entry_nl{nl}"),
            "-c".into(),
            "cgv/gv.c".into(),
            "-o".into(),
            obj.display().to_string(),
        ]);
        objects.push(obj.display().to_string());
    }
    let main_obj = out.join("main.o").display().to_string();
    run(&[
        "-c".into(),
        "cgv/main.c".into(),
        "-o".into(),
        main_obj.clone(),
    ]);
    let exe = out
        .join(if windows { "cgv.exe" } else { "cgv" })
        .display()
        .to_string();
    let mut link = vec!["-o".to_string(), exe, main_obj];
    link.extend(objects);
    link.extend(["-lpthread".to_string(), "-lm".to_string()]);
    run(&link);
}
