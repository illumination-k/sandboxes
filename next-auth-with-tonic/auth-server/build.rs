fn main() {
    tonic_build::configure()
        .compile(&["auth.proto"], &["proto"]).unwrap();
}