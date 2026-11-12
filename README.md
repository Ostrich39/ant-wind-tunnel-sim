# ant-wind-tunnel-sim

A desktop program, written in Rust, to design and analyze airfoils.

- Goals and stages: [GOALS.md](GOALS.md)
- Sources: [REFERENCES.md](REFERENCES.md)

## Status

![Program window](docs/screenshot.png)

- Stage 1 is complete: an inviscid (no friction) panel method for NACA 4-digit
  airfoils, with validation tests against exact theory.
- Stage 2 is complete: a desktop window with an airfoil view, a pressure
  plot, and live controls.

The model has no friction. It gives no drag value and it does not predict stall.

## Use

Install Rust from https://rustup.rs, then:

```
cargo run --release -p ant-wind-tunnel          # start the desktop program
cargo test                                      # validation tests
cargo run --release --example polar -- 2412     # results table for NACA 2412
```

## License

MIT. See [LICENSE](LICENSE).
