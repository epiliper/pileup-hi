build:
	cargo build --release && cp target/release/viggo $$HOME/.cargo/bin

test:
	cargo test -- --show-output
