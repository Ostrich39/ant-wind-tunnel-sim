# Project goals

A desktop program, written in Rust, to design and analyze airfoils.
The project is for personal use and learning. It is not for sale.

## Function goals

1. **Analyze 2D airfoils.** Calculate lift, drag, and moment coefficients for
   a given shape, angle of attack, and Reynolds number.
2. **Design airfoils.** Let the user edit a shape and see the results change.
3. **Show the flow.** Show pressure and velocity around the airfoil.
4. **Validate results.** Compare each solver stage with theory or with
   published wind tunnel data. Automated tests do this comparison.

## Performance goals

| Goal | Target |
|---|---|
| High-quality render | Sharp lines and smooth color fields at any window size, drawn on the GPU |
| Small file size | Desktop program below 20 MB |
| Smooth operation | 60 frames per second while the user moves or edits a shape |
| Quick results | New results in less than 0.1 second after each change |

## Decisions

- **Platform:** desktop app. A browser version is not in the scope at this time.
- **Window library:** egui and eframe with the wgpu renderer (GPU).
- **Update mode:** live. The program calculates a new result when an input changes.
- **Language:** Rust for all core calculations.
- **Order of work:** airfoils first (panel method). Arbitrary shapes later
  (live flow solver, labeled as approximate).
- **License:** MIT (see `LICENSE`). No code is copied from GPL software such as XFOIL or XFLR5. Methods are written from the published papers.

## Stages

| Stage | Content | Status |
|---|---|---|
| 1 | Inviscid panel method, NACA 4-digit shapes, validation tests | Done |
| 2 | Desktop window: airfoil view, pressure plot, live controls for angle and NACA shape | Done |
| 3a | One-way boundary layer: friction drag, transition, Reynolds number, separation warning | Done |
| 3b | Two-way (coupled) boundary layer from Drela 1989: lift correction, separation bubbles, stall | Not started |
| 4 | Shape editor and airfoil file import (.dat) | Not started |
| 5 | Flow field display (streamlines, pressure colors) | Not started |
| 6 | Live flow solver for arbitrary shapes (approximate) | Not started |

## Limits of the current model (stages 1 to 3a)

- The boundary layer calculation is one-way. It uses the pressure from the
  friction-free solution and does not change it.
- Thus the lift has no friction correction. Lift values are too high
  compared with wind tunnel data, and the program cannot predict stall.
- The drag value is valid only for attached flow. If the turbulent boundary
  layer separates before 98% of the chord, the program gives no drag value
  and shows where the flow separates.
- The laminar separation bubble is not modeled. Laminar separation starts
  turbulent flow at the same position.
- The drag is compared with wind tunnel data for one airfoil (NACA 0012) at
  one Reynolds number (6 million). See the validation table.
- Stage 3b corrects the first four limits.

## Assumptions in stage 3a

| Assumption | Value | Basis |
|---|---|---|
| Shape factor at the start of turbulent flow | 1.4 | Assumption. No source is recorded |
| Shape factor at turbulent separation | 2.4 | Published range is 1.8 to 2.4 (REFERENCES.md 1.11) |
| Trip strip position for the Ladson comparison | 5% of chord | Assumption. The position is not in the data file |
| Drag reference station | 98% of chord | Assumption. The friction-free flow has a stagnation point at the trailing edge that a real flow does not have |

## Validation status of stage 3a

| Quantity | Compared with | Result |
|---|---|---|
| Laminar flat plate: thickness, friction, shape factor | Blasius exact solution | Difference 1.0%, 1.2%, 0.8% |
| Stagnation flow thickness | Thwaites closed form; Hiemenz exact solution | Equal to the closed form; 6% below the exact solution (known property of Thwaites' method) |
| Decelerating (Howarth) flow | Thwaites closed form; published separation position | Equal to the closed form; separation at x = 0.123 (published 0.123; exact 0.1199) |
| Turbulent flat plate drag | Prandtl-Schlichting formula | Difference 7% at Re 1e6, 6% at Re 1e7, 12% at Re 1e8 |
| NACA 0012 drag with trip strip, Re 6 million, angles −4° to 8.3° | Ladson wind tunnel data (26 points) | Difference −5.4% to +4.1%. The scatter between the three data sets is approximately 4% |
| NACA 0012 drag, free transition, Re 6 million, lift coefficient to ±0.46 | Abbott and von Doenhoff data (digitized) | Difference −1.1% to +5.3% |
| NACA 0012 drag, free transition, lift coefficient 0.55 to 1.08 | Abbott and von Doenhoff data (digitized) | Model is 8% to 17% low |
| NACA 0012 lift at the same angle, 2° to 8.3° | Ladson wind tunnel data | Model is 9% to 17% high (no friction correction) |
| Drag at angles above approximately 9° | Ladson wind tunnel data | No model value. The data continue to 17° |
| Other airfoils, other Reynolds numbers | Wind tunnel data | Not compared |
| Airfoil drag, internal checks | Symmetry, trend with Reynolds number, panel count | Pass. Change of 1% to 2% when the panel count doubles |

## Status against the performance goals (stage 2)

| Goal | Target | Measured | Status |
|---|---|---|---|
| Small file size | Below 20 MB | 10.1 MB (Linux release build) | Met on Linux. Not measured on Windows or macOS |
| Quick results | Below 0.1 s | 2 ms to 3 ms at 160 panels. An automated test applies the limit at 400 panels | Met |
| Smooth operation | 60 frames per second | Not measured. The test computer has no GPU | Open |
| High-quality render | Sharp at any window size | Lines are vector shapes with anti-aliasing. Checked by eye on one screen capture | Partly checked |

## Validation status of stage 1

| Quantity | Compared with | Result |
|---|---|---|
| Lift | Exact Joukowski solution | Error 0.01% (symmetric, 200 panels); 0.6% to 1.3% (cambered, 400 panels) |
| Surface pressure, front 90% of chord | Exact Joukowski solution | Largest error 0.04 in pressure coefficient (400 panels) |
| Surface pressure, last 10% of chord | Exact Joukowski solution | Not validated. The error is large at the cusp trailing edge of the Joukowski shape |
| NACA 4-digit shapes | No exact solution is available | Shape checked against the defining equation. Flow results checked only for internal consistency |
| Wind tunnel data | Not compared | Stage 3 |

## Rules for changes

- Each change must pass `cargo fmt --check`, `cargo clippy` with no
  warnings, and `cargo test`. GitHub runs these checks automatically.
- Each number in a test must have a source in `REFERENCES.md` or come from
  an exact solution.
- Each new source goes into `REFERENCES.md` in the same change.

## Open items for stage 3a

- The trip strip position in the Ladson tests is not confirmed. The model
  uses 5% of chord.
- Compare a cambered airfoil and a lower Reynolds number with wind tunnel data.

## Open decisions

- Accuracy target against wind tunnel data (for example, lift within 5%).

## References

- Drela, M. (1989). XFOIL: An Analysis and Design System for Low Reynolds
  Number Airfoils. https://doi.org/10.1007/978-3-642-84010-4_1
- Drela, M., Giles, M. B. (1987). Viscous-Inviscid Analysis of Transonic and
  Low Reynolds Number Airfoils. AIAA Journal 25(10), 1347-1355.
- Kuethe, A. M., Chow, C.-Y. Foundations of Aerodynamics, 5th ed., section 5.10
  (linear-vortex panel method used in stage 1).
