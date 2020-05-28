extern crate cc;

// is_x86_feature_detected("adx")
// is_x86_feature_detected("bmi2")

fn main() {
    println!("cargo:rerun-if-changed=./src/mul_4.S");
    
    let target_arch = std::env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    if target_arch == "x86_64" {
        cc::Build::new()
            .flag("-c")
            .file("./src/mul_4.S")
            .compile("libff-derive-crypto.a");
    }
}