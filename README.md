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

- Stage 3a is complete: a one-way boundary layer calculation (not used by
  the program after stage 3b).
- Stage 3b is complete for attached flow: the flow and the boundary layer
  are solved together (Drela's method). The program shows lift with
  friction, drag, the start of turbulent flow, and flow separation.

Accuracy for the NACA 0012 at a Reynolds number of 6 million, with a trip
strip, against wind tunnel data: lift within 8%, drag within 8%, for angles
to 14°. The method is not reliable near stall and after stall. See
[GOALS.md](GOALS.md) for the full validation table and the limits.

## Use

Install Rust from https://rustup.rs, then:

```
cargo run --release -p ant-wind-tunnel          # start the desktop program
cargo test                                      # validation tests
cargo run --release --example polar -- 2412     # results table for NACA 2412
```

## License

MIT. See [LICENSE](LICENSE).
