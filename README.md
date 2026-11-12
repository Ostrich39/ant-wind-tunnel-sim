# ant-wind-tunnel-sim

A desktop program, written in Rust, to design and analyze airfoils.

- Goals and stages: [GOALS.md](GOALS.md)
- Sources: [REFERENCES.md](REFERENCES.md)

## Status

Stage 1 is complete: an inviscid (no friction) panel method for NACA 4-digit
airfoils, with validation tests against exact theory.

## Use

Install Rust from https://rustup.rs, then:

```
cargo test                                      # validation tests
cargo run --release --example polar -- 2412     # results table for NACA 2412
```

## License

MIT. See [LICENSE](LICENSE).
