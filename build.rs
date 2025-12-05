fn main() {
    if std::env::var("CARGO_FEATURE_GPU").is_err() {
        println!("cargo:warning=GPU feature disabled → skipping CUDA build");
        return;
    }

    println!("cargo:warning=GPU feature enabled → checking nvcc...");

    let nvcc = match which::which("nvcc") {
        Ok(path) => path,
        Err(_) => {
            println!("cargo:warning=CUDA toolkit not found → GPU backend disabled");
            return;
        }
    };

    println!("cargo:rerun-if-changed=src/annotate/cuda/ani_kernel.cu");

    let out_dir = std::env::var("OUT_DIR").unwrap();
    let output_ptx = format!("{}/ani_kernel.ptx", out_dir);

    println!("cargo:warning=Compiling CUDA kernel with nvcc: {:?}", nvcc);

    let status = std::process::Command::new(nvcc)
        .args(&[
            "-ptx",
            "src/annotate/cuda/ani_kernel.cu",
            "-o",
            &output_ptx,
            "-arch=sm_75",
            "-O3",
        ])
        .status()
        .expect("Failed to run nvcc");

    if !status.success() {
        panic!("nvcc failed to compile ani_kernel.cu");
    }

    println!("cargo:warning=CUDA kernel generated: {}", output_ptx);
}
