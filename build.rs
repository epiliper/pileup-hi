fn main() {
    cc::Build::new()
        .file("src/overlap_wrapper.c")
        .include("htslib") // Just use "htslib" - the linker will find it
        .opt_level(3)
        .flag_if_supported("-march=native")
        .compile("overlap_wrapper");
}
