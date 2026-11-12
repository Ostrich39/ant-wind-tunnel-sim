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
| Small file size | Desktop program below 35 MB (changed from 20 MB on 12 November 2026) |
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
| 3b-1 | Source panels for the displacement effect; fast repeated solution of the panel equations | Done |
| 3b-2 | Two-way calculation by direct iteration (experimental) | Removed. It was not accurate enough. Stage 3b-3 replaces it |
| 3b-3 | Wake panels and simultaneous solution of the flow and the boundary layer (Drela's method): lift with friction, drag, mild separation | Done. The program uses it |
| 3b-4a | Solution at high angles: start march with inverse mode, moving stagnation point, update limits of the paper | Done |
| 3b-4b | Prediction of abrupt stall; separation bubbles; transition between stations | Not started |
| 4 | Shape editor and airfoil file import (.dat) | Not started |
| 5 | Flow field display (streamlines, pressure colors) | Not started |
| 6 | Live flow solver for arbitrary shapes (approximate) | Not started |

## Limits of the current model (stage 3b-3, used by the program)

- **Stall.** The calculation does not predict an abrupt stall. It gives a
  slow growth of the separation from the trailing edge, and its lift does
  not decrease. For the NACA 0012 at a Reynolds number of 6 million, the
  wind tunnel lift decreases from 1.6 to 1.0 between 17° and 18°. The model
  lift at 18° is 1.74 (74% too high), and its drag is 87% too low. The
  window shows a warning when the flow separates before 95% of the chord.
- **Near maximum lift** (15° to 17° for that airfoil): the lift is 4% to 6%
  too high and the drag is 3% to 17% too low.
- **No solution.** For some conditions in deep stall there is no solution
  (for example a thin airfoil at a high angle). The window then shows a
  message and no drag value.
- **Separation bubbles.** A laminar separation bubble is not modeled.
  Turbulent flow starts at the laminar separation point or before it. This
  is most important at Reynolds numbers below approximately 500 000.
- **Transition position.** Turbulent flow starts at a station (a panel
  midpoint), not at a position between two stations. In one solution the
  start of turbulent flow can only move toward the leading edge.
- **Free transition.** The drag with free transition depends on the
  transition factor N, which describes the turbulence of the wind tunnel.
  With the standard value N = 9, the drag is 5% to 17% below the Abbott and
  von Doenhoff data. See the validation table.
- **Compared airfoils.** Wind tunnel data are compared for one airfoil
  (NACA 0012) at one Reynolds number (6 million).
- **Speed range.** Low speed only (no compressibility).
- **Trailing edge.** Sharp trailing edge only.

## Differences from XFOIL

The method is the same in principle (Drela 1989). These parts are different:

| Part | This program | XFOIL |
|---|---|---|
| Boundary layer stations | Panel midpoints | Panel end points |
| Transition | At a station | At a position between two stations |
| Jacobian matrix | Finite differences, full matrix solution | Analytic, block solution |
| Stagnation point | Moves to a different panel when the surface speed near it changes sign. Its position inside a panel is fixed | Moves in each iteration |
| Amplification equation | With the modification of Fidkowski (mfoil) | Original |
| Panel method | Linear vortex panels, velocity form | Linear vortex panels, stream function form |

No XFOIL code and no mfoil code is in this program. The equations are from
the papers (REFERENCES.md, items 2.1, 2.2, 2.4).

## Assumptions in stage 3b-3

| Assumption | Value | Basis |
|---|---|---|
| Wake length | 1 chord, 20 panels | Assumption (mfoil uses 1 chord) |
| Trip strip is effective only if the momentum-thickness Reynolds number is 200 or more | 200 | Assumption. No source is recorded |
| Trip strip position for the Ladson comparison | 5% of chord | Assumption. The position is not in the data file |
| Edge speed for the start march in the last 10% of the chord | Straight continuation of the speed between 80% and 90% of the chord | Assumption. The friction-free flow has a stagnation point at the trailing edge. This changes the start values only, not the solution |
| Laminar separation starts turbulent flow | Shape factor 3.8 (limit of the paper for the direct march) | No separation bubble model. The rule applies near convergence, or after 10 iterations |
| Panels at the trailing edge | Approximately 1% of chord (120 panels) | Smaller panels made the calculation unstable |

## Validation status of stage 3b-3

NACA 0012, Reynolds number 6 million, 120 panels.

| Quantity | Compared with | Result |
|---|---|---|
| Lift with trip strip, angles −4° to 14.3° (38 points) | Ladson wind tunnel data | Largest difference 8.2% (at 2°). At 4° and above the model is 0.5% to 5% high. Stage 3a was 9% to 17% high |
| Drag with trip strip, angles −4° to 14.3° (38 points) | Ladson wind tunnel data | −8.7% to +5.1%. The scatter between the three data sets is approximately 4% |
| Drag with free transition, N = 9 (18 points) | Abbott and von Doenhoff data (digitized) | −17.7% to +0.8%. Most points are 10% to 15% low |
| Effect of N on drag at 0°, free transition | Abbott data value 0.0059 | N = 9: 0.0049. N = 7: 0.0053. N = 5: 0.0058. N = 3: 0.0064 |
| Lift and drag with trip strip, 15.2° to 17.2° (9 points, to the maximum lift) | Ladson wind tunnel data | Lift 3.7% to 5.5% high. Drag 3% to 17% low |
| Lift with trip strip at 18° (after the stall) | Ladson wind tunnel data: 1.00 | Model: 1.74. The model does not predict the stall |
| Conditions with a solution from the friction-free start: 7 airfoil and Reynolds number sets, 13 angles from −4° to 20° | Number of solutions | 77 of 91. Before stage 3b-4a: 60 of 91. The 14 conditions without a solution are in deep stall |
| Iterations for one solution | Same 91 conditions | Median 5, 90% in 20 or less, largest 48 |
| Panel count 120 and 200 | Same case | Lift changes 0.1%, drag changes 1.3% |
| Laminar flat plate | Blasius exact solution | Thickness, friction, and shape factor within 1% |
| Stagnation flow | Hiemenz exact solution | Thickness within 0.3% (stage 3a: 6%) |
| Decelerating (Howarth) flow | Exact separation position 0.1199 | Shape factor limit reached at 0.1236 |
| Turbulent flat plate drag | Prandtl-Schlichting formula | 5% to 7% low at Reynolds numbers of 10 and 100 million |
| Other airfoils, other Reynolds numbers | Wind tunnel data | Not compared |

Speed: one solution takes 30 ms to 60 ms for 120 panels in the normal
angle range (target: below 100 ms). Near stall it can take 150 ms to
350 ms. A condition without a solution takes approximately 400 ms in the
window.

The calculation runs on a second thread. The window does not wait for
it: the last result stays on the screen until the new result is ready,
and only the newest inputs are calculated. Thus a slow calculation
delays the numbers, but it does not stop the window.

## Check of the equations against the paper (12 November 2026)

The closure relations and the residual equations in the code were compared
with the pages of the paper (REFERENCES.md, item 2.4: equations 5 to 13, 17,
and 45 to 79, and Table 2).

| Item | Result |
|---|---|
| Five details that were first completed from memory (turbulent H*, tanh argument in the skin friction, outer dissipation term, factor 2δ in the shear lag equation, constant G_B in the slip velocity) | All five agree with the paper |
| All other closure relations and all constants of Table 2 | Agree with the paper |
| Exponent limit in the turbulent skin friction (equation 58) | Was different (hard limit). Corrected to the smooth limit of the paper. No change in the test results |
| Extra amplification term (equation 76) | Was not in the code. Added. No change in the test results |
| Wake dissipation (equation 66) | The code had a factor of 2 from memory. The paper has no such factor. Corrected to the paper. Drag changed by 0.5% to 1.6%, lift by less than 0.5% |

Parts of the method that are different from the paper on purpose:

| Part | Paper | This program | Cause |
|---|---|---|---|
| Transition (equations 16 to 19) | At a position between two stations | At a station | The first version of the paper's method in this program made convergence worse. This is the work of stage 3b-4 |
| First wake station (equation 14) | Sum of the two boundary layers at the trailing edge point | Wake equations from the last airfoil stations to the first wake station | The stations of this program are panel midpoints, thus the last airfoil stations and the first wake station are at different positions |
| First station after the stagnation point | Extrapolation of the first two nodes | Edge speed proportional to the distance, slope from the second station | Simpler, and stable when the stagnation point moves |
| Start values (section V.A) | March with direct mode, then inverse mode, then extrapolation | The same, and a smooth edge speed near the trailing edge; laminar separation starts turbulent flow | The stations of this program are near a trailing edge stagnation point that the node speeds of the paper do not have |
| Update limits (section V.B) | One factor for all unknowns | The same limits | No difference |
| Number of wake panels | N/10 + 10 (22 for 120 panels) | 20 | No reason. Not important for the results |

## Stage 3a (one-way calculation)

The solver still contains the one-way calculation (`boundary_layer.rs`,
`viscous.rs`) and its tests. The program does not use it.

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
| NACA 0012 drag with trip strip, Re 6 million, angles −4° to 8.3° | Ladson wind tunnel data (23 points) | Difference −5.4% to +4.1%. The scatter between the three data sets is approximately 4% |
| NACA 0012 drag, free transition, Re 6 million, lift coefficient to ±0.46 | Abbott and von Doenhoff data (digitized) | Difference −1.1% to +5.3% |
| NACA 0012 drag, free transition, lift coefficient 0.55 to 1.08 | Abbott and von Doenhoff data (digitized) | Model is 8% to 17% low |
| NACA 0012 lift at the same angle, 2° to 8.3° | Ladson wind tunnel data | Model is 9% to 17% high (no friction correction) |
| Drag at angles above approximately 9° | Ladson wind tunnel data | No model value. The data continue to 17° |
| Other airfoils, other Reynolds numbers | Wind tunnel data | Not compared |
| Airfoil drag, internal checks | Symmetry, trend with Reynolds number, panel count | Pass. Change of 1% to 2% when the panel count doubles |

## Status against the performance goals (stage 2)

| Goal | Target | Measured | Status |
|---|---|---|---|
| Small file size | Below 35 MB | 10.2 MB (Linux release build) | Met on Linux. Not measured on Windows or macOS |
| Quick results | Below 0.1 s | 2 ms to 3 ms at 160 panels. An automated test applies the limit at 400 panels | Met |
| Smooth operation | 60 frames per second | Not measured. The test computer has no GPU. The calculation does not stop the window (second thread, automated test) | Open |
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

## Result of stage 3b-2 (experimental two-way calculation, removed)

The calculation puts a source on each panel to model the displacement of
the outer flow by the boundary layer, and repeats the flow and boundary
layer solutions until the lift is stable.

Measured for the NACA 0012 with a trip strip, Re 6 million, against the
Ladson data:

| Quantity | One-way (stage 3a) | Two-way (stage 3b-2) |
|---|---|---|
| Lift at the same angle, 2° to 11° | 9% to 17% high | 2% to 13% high |
| Drag difference | −5.4% to +4.1%, angles to 8.3° | −6.2% to +4.8%, angles to 11.2° |
| Largest angle with a drag value | Approximately 9° | Approximately 11° |
| Time for one result (200 panels) | 3 ms | Approximately 10 ms |

Why the program does not use it:

- The lift depends on the cut-off station for the sources, which is an
  arbitrary value. Lift error at 4° and 8° for different cut-off values:
  90% of chord: +10% and +9%. 95% of chord: +8% and +7%. 98% of chord: +6%
  and +4%. At 99.5% of chord the iteration does not converge.
- To select the cut-off value that agrees best with the data is curve
  fitting. It is not a prediction. The code uses 95%, which was selected
  before the comparison.
- The cause is known: the model has no wake, and the direct iteration is
  not stable near the trailing edge and near separation. Stage 3b-3 removes
  the cut-off.

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
