# References

This file lists each source used or planned for this project.
Add a new source here in the same change that uses it.

Status words:
- **Used**: the code or tests use this source now.
- **Planned**: a later stage will use this source.
- **Background**: for comparison or reading only.

Check words:
- **Checked**: the citation was found in a search of the web.
- **Not checked**: the citation is from memory. Confirm it before you cite it in other work.

## 1. Theory used in the code (stage 1)

| # | Source | Used for | Where | Check |
|---|---|---|---|---|
| 1.1 | Kuethe, A. M., Chow, C.-Y. *Foundations of Aerodynamics*, 5th ed., Wiley, 1998, section 5.10. | Linear-vortex panel method | `crates/solver/src/panel.rs` | Not checked |
| 1.2 | Jacobs, E. N., Ward, K. E., Pinkerton, R. M. *The Characteristics of 78 Related Airfoil Sections from Tests in the Variable-Density Wind Tunnel*. NACA Report 460, 1933. | NACA 4-digit thickness and camber equations | `crates/solver/src/geometry.rs` | Not checked |
| 1.3 | Abbott, I. H., von Doenhoff, A. E. *Theory of Wing Sections*. Dover, 1959. | NACA 4-digit equations (second source); wind tunnel data for later validation | `crates/solver/src/geometry.rs` | Not checked |
| 1.4 | Joukowski transformation and Kutta-Joukowski theorem. Standard theory, for example Anderson, J. D. *Fundamentals of Aerodynamics*, McGraw-Hill. | Exact lift solution for the validation tests; lift from circulation | `crates/solver/tests/validation.rs`, `panel.rs` | Not checked |
| 1.5 | Thick-airfoil lift slope estimate, 2π(1 + 0.77 t). | Approximate check of NACA 0012 lift | `crates/solver/tests/validation.rs` | Source not identified |

## 2. Theory planned for later stages

| # | Source | Planned use | Check |
|---|---|---|---|
| 2.1 | Drela, M. (1989). XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils. In: Mueller, T. J. (ed.), *Low Reynolds Number Aerodynamics*, Lecture Notes in Engineering 54, Springer. https://doi.org/10.1007/978-3-642-84010-4_1 | Stage 3: boundary layer, transition, coupling with the panel method | Checked |
| 2.2 | Drela, M., Giles, M. B. (1987). Viscous-Inviscid Analysis of Transonic and Low Reynolds Number Airfoils. *AIAA Journal* 25(10), 1347-1355. | Stage 3: boundary layer equations | Checked |
| 2.3 | Drela, M. (1989). An Integral Boundary Layer Formulation for Blunt Trailing Edges. AIAA paper. | Stage 3: blunt trailing edges | Checked, but the paper number is not certain: one source gives AIAA-89-2166, a different source gives AIAA-89-2200 |
| 2.4 | Fidkowski, K. mfoil: coupled integral boundary layer and vortex-panel solver (Matlab). https://websites.umich.edu/~kfid/codes/mfoil/mfoil.pdf | Stage 3: modern description of the XFOIL-type method | Checked |
| 2.5 | Griebel, M., Dornseifer, T., Neunhoeffer, T. (1998). *Numerical Simulation in Fluid Dynamics: A Practical Introduction*. SIAM. | Stage 6: Navier-Stokes solver on a staggered grid | Checked |
| 2.6 | Succi, S. (2001). *The Lattice Boltzmann Equation for Fluid Dynamics and Beyond*. Oxford University Press. | Stage 6: lattice Boltzmann method (alternative) | Not checked |

## 3. Software

| # | Software | License | Link | Status |
|---|---|---|---|---|
| 3.1 | XFOIL (Drela, Youngren) | GPL | https://web.mit.edu/drela/Public/web/xfoil | Background. Comparison results for validation |
| 3.2 | XFLR5 | GPL v3 | https://sourceforge.net/projects/xflr5/ | Background |
| 3.3 | QBlade | Open source | Not recorded | Background |
| 3.4 | OpenFOAM | GPL v3 | https://openfoam.org | Background |
| 3.5 | ANSYS Fluent | Commercial | Not recorded | Background |
| 3.6 | RWIND (Dlubal) | Commercial | https://www.dlubal.com | Background |
| 3.7 | flow2d-rs (Rust, Navier-Stokes) | Not recorded | https://github.com/1n0r1/flow2d-rs | Background for stage 6 |
| 3.8 | stroemung (Rust, CFD) | Not recorded | https://github.com/wickedchicken/stroemung | Background for stage 6 |
| 3.9 | bloe (Rust, lattice Boltzmann) | MIT | https://github.com/ndbaker1/bloe | Background for stage 6 |
| 3.10 | electra (Rust and WebAssembly, lattice Boltzmann) | Apache-2.0 or MIT | https://github.com/s2011r2593/electra | Background for stage 6 |
| 3.11 | Flexcompute Foil (Rust reimplementation of XFOIL; source release announced) | GPL (announced) | https://hs.flexcompute.com/blog/flexfoil | Background for stage 3 |

No code from this software is copied into this project.

## 4. Validation studies

| # | Source | Content | Check |
|---|---|---|---|
| 4.1 | Řezníček, T. (2023). *Open source and commercial CFD codes comparison on 2D airfoil analysis*. Bachelor's thesis, Brno University of Technology. https://theses.cz/id/nw0p3h/TomasReznicek_BP_2023.pdf | OpenFOAM and ANSYS Fluent compared with wind tunnel data | Checked |

## 5. Items to complete

- Confirm each "Not checked" citation (edition, year, section number).
- Find the source of item 1.5, or replace that test with wind tunnel data.
- Resolve the paper number of item 2.3.
- Select the wind tunnel data set for NACA 0012 validation in stage 3.
