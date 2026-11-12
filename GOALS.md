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
- **Language:** Rust for all core calculations.
- **Order of work:** airfoils first (panel method). Arbitrary shapes later
  (live flow solver, labeled as approximate).
- **License:** MIT (see `LICENSE`). No code is copied from GPL software such as XFOIL or XFLR5. Methods are written from the published papers.

## Stages

| Stage | Content | Status |
|---|---|---|
| 1 | Inviscid panel method, NACA 4-digit shapes, validation tests | Done |
| 2 | Desktop window: draw the airfoil, pressure plot, angle-of-attack control | Not started |
| 3 | Boundary layer (friction drag, transition, Reynolds number) from Drela 1989 | Not started |
| 4 | Shape editor and airfoil file import (.dat) | Not started |
| 5 | Flow field display (streamlines, pressure colors) | Not started |
| 6 | Live flow solver for arbitrary shapes (approximate) | Not started |

## Limits of stage 1

- The solver ignores friction. Thus the drag value is not physical; theory
  says it is zero, and the printed value shows the numerical error.
- Lift values are too high compared with wind tunnel data, and the solver
  cannot predict stall. Stage 3 corrects this.

## Open decisions

- Accuracy target against wind tunnel data (for example, lift within 5%).
- Live update while the user drags a shape, or update after each change.

## References

- Drela, M. (1989). XFOIL: An Analysis and Design System for Low Reynolds
  Number Airfoils. https://doi.org/10.1007/978-3-642-84010-4_1
- Drela, M., Giles, M. B. (1987). Viscous-Inviscid Analysis of Transonic and
  Low Reynolds Number Airfoils. AIAA Journal 25(10), 1347-1355.
- Kuethe, A. M., Chow, C.-Y. Foundations of Aerodynamics, 5th ed., section 5.10
  (linear-vortex panel method used in stage 1).
