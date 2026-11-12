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
| 1.1 | Kuethe, A. M., Chow, C.-Y. *Foundations of Aerodynamics: Bases of Aerodynamic Design*, 5th ed., Wiley, section 5.10, page 156. | Linear-vortex panel method | `crates/solver/src/panel.rs` | Checked (title, edition, section, page). Publication year not checked |
| 1.2 | Jacobs, E. N., Ward, K. E., Pinkerton, R. M. *The Characteristics of 78 Related Airfoil Sections from Tests in the Variable-Density Wind Tunnel*. NACA Report 460, 1933. https://ntrs.nasa.gov/citations/19930091108 | NACA 4-digit thickness and camber equations | `crates/solver/src/geometry.rs` | Checked |
| 1.3 | Abbott, I. H., von Doenhoff, A. E. *Theory of Wing Sections: Including a Summary of Airfoil Data*. Dover, 1959. ISBN 9780486605869. | NACA 4-digit equations (second source); wind tunnel data for later validation | `crates/solver/src/geometry.rs` | Checked |
| 1.4 | Joukowski transformation and Kutta-Joukowski theorem. Standard theory, for example Anderson, J. D. *Fundamentals of Aerodynamics*, McGraw-Hill. | Exact lift solution for the validation tests; lift from circulation | `crates/solver/tests/validation.rs`, `panel.rs` | Not checked |

### Boundary layer (stage 3a)

| # | Source | Used for | Where | Check |
|---|---|---|---|---|
| 1.5 | Thwaites, B. (1949). Approximate calculation of the laminar boundary layer. *Aeronautical Quarterly* 1, 245. | Laminar boundary layer method | `crates/solver/src/boundary_layer.rs` | Checked (citation). Paper not read |
| 1.6 | Thwaites correlations H(λ) and l(λ), as given in: Nordin, A., Lim, C. W. Thwaites' method in laminar boundary layer. *Jurnal Teknologi*, UTM. https://journals.utm.my/jurnalteknologi/article/download/995/965/1005 | Correlation constants | `boundary_layer.rs` | Checked, with one difference: that paper prints 1.420 where this code uses 1.402. The value 1.402 makes the slope continuous at λ = 0 (a test checks this) and agrees with item 1.8 from memory |
| 1.7 | Sankar, L. *Boundary Layer Analyses*, AE 3903/4903 Airfoil Design, Georgia Institute of Technology (course slides). https://www.slideshare.net/slideshow/boundarylayeranalysis/4360255 | Michel criterion, Head's method inverse relation, Ludwieg-Tillmann law, separation at λ = −0.09 | `boundary_layer.rs` | Checked |
| 1.8 | Moran, J. *An Introduction to Theoretical and Computational Aerodynamics*. Dover, 2003. ISBN 9780486428796. Chapter 7. | Thwaites, Head, transition, drag method (the same set of methods) | `boundary_layer.rs`, `viscous.rs` | Checked (book and chapter list). Equations from memory: Head's entrainment function 0.0306 (H1 − 3)^−0.6169 and the Squire-Young exponent (H + 5)/2 |
| 1.9 | Head, M. R. (1958). Entrainment in the turbulent boundary layer. ARC R&M 3152. | Turbulent boundary layer method | `boundary_layer.rs` | Not checked |
| 1.10 | Ludwieg, H., Tillmann, W. (1949). Skin friction law for turbulent boundary layers. | Turbulent skin friction | `boundary_layer.rs` | Checked (authors and year only) |
| 1.11 | Cebeci, T., Mosinskis, G. J., Smith, A. M. O. Calculation of separation points in incompressible turbulent flows. *Journal of Aircraft*. https://typeset.io/pdf/calculation-of-separation-points-in-incompressible-turbulent-4llvp2elsu.pdf | Shape factor at turbulent separation, 1.8 to 2.4 | `boundary_layer.rs` | Checked (statement). Author list and year not checked |
| 1.12 | Michel, R. (1951). ONERA report. | Transition criterion | `boundary_layer.rs` | Not checked. The formula is from item 1.7 |
| 1.13 | Squire, H. B., Young, A. D. (1938). The calculation of the profile drag of aerofoils. ARC R&M 1838. | Drag formula | `viscous.rs` | Checked (authors and year). Report number not checked |
| 1.14 | Stern, F. *Momentum Integral Methods*, course notes, University of Iowa. https://stern.lab.uiowa.edu/sites/stern.lab.uiowa.edu/files/2025-11/Chapter-7.4-Momentum-Integral-Methods.pdf | Test values: Howarth flow closed form, separation at 0.123, exact 0.1199 | `crates/solver/tests/boundary_layer.rs` | Checked |
| 1.15 | Blasius flat plate solution (θ = 0.664 x/√Re, cf = 0.664/√Re, δ* = 1.72 x/√Re) and Prandtl-Schlichting formula CF = 0.455/(log10 Re)^2.58. Course notes, University of Iowa. https://user.engineering.uiowa.edu/~fluids/posting/Schedule/Turbulent_BL_12-02-2015.pdf | Test values | `crates/solver/tests/boundary_layer.rs` | Checked |
| 1.16 | Hiemenz stagnation flow, θ = 0.292 √(ν/k). | Test value | `crates/solver/tests/boundary_layer.rs` | Not checked (from memory) |

### Displacement sources (stage 3b)

| # | Source | Used for | Where | Check |
|---|---|---|---|---|
| 1.17 | Velocity of a constant-strength source panel. Standard result, for example Katz, J., Plotkin, A. *Low-Speed Aerodynamics*, Cambridge University Press; Hess, J. L., Smith, A. M. O. (1967). | Source panels | `crates/solver/src/panel.rs` | Not checked (from memory). Tests check the equations: flow through a far circle and flow out of the surface |
| 1.18 | Displacement effect as a surface source with strength d(ue δ*)/ds (Lighthill 1958, transpiration model), as used in XFOIL (item 2.1). | Two-way calculation | `crates/solver/src/coupled.rs` | Not checked (from memory) |

## 2. Theory for the two-way calculation and for later stages

| # | Source | Planned use | Check |
|---|---|---|---|
| 2.1 | Drela, M. (1989). XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils. In: Mueller, T. J. (ed.), *Low Reynolds Number Aerodynamics*, Lecture Notes in Engineering 54, Springer. https://doi.org/10.1007/978-3-642-84010-4_1 | **Used** in stage 3b-3 (method). The equations are from item 2.4 | Checked (citation). Paper not read |
| 2.2 | Drela, M., Giles, M. B. (1987). Viscous-Inviscid Analysis of Transonic and Low Reynolds Number Airfoils. *AIAA Journal* 25(10), 1347-1355. | **Used** in stage 3b-3 (boundary layer formulation). The equations are from item 2.4 | Checked |
| 2.3 | Drela, M. (1989). An Integral Boundary Layer Formulation for Blunt Trailing Edges. AIAA paper. | Stage 3b: blunt trailing edges | Checked, but the paper number is not certain: one source gives AIAA-89-2166, a different source gives AIAA-89-2200 |
| 2.4 | Fidkowski, K. (2022). A Coupled Inviscid-Viscous Airfoil Analysis Solver, Revisited (mfoil). *AIAA Journal* (accepted November 2021). https://public.websites.umich.edu/~kfid/codes/mfoil/mfoil.pdf | **Used** in stage 3b-3: boundary layer equations, closure relations, amplification rate, constants (`crates/solver/src/closure.rs`, `ibl.rs`, `interaction.rs`) | Checked against the pages of the paper on 12 November 2026 (equations 5 to 13, 17, 45 to 79, Table 2). Result: GOALS.md, section "Check of the equations against the paper". Journal volume and pages not checked |
| 2.5 | Griebel, M., Dornseifer, T., Neunhoeffer, T. (1998). *Numerical Simulation in Fluid Dynamics: A Practical Introduction*. SIAM. | Stage 6: Navier-Stokes solver on a staggered grid | Checked |
| 2.6 | Succi, S. (2001). *The Lattice Boltzmann Equation for Fluid Dynamics and Beyond*. Oxford University Press. https://doi.org/10.1093/oso/9780198503989.001.0001 | Stage 6: lattice Boltzmann method (alternative) | Checked |

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

### Libraries that the program uses

| # | Library | License | Link | Used for |
|---|---|---|---|---|
| 3.12 | egui and eframe 0.37 | MIT or Apache-2.0 | https://github.com/emilk/egui | Window, controls, and drawing (`crates/app`) |
| 3.13 | wgpu (through eframe) | MIT or Apache-2.0 | https://github.com/gfx-rs/wgpu | GPU rendering |

## 4. Validation studies

| # | Source | Content | Check |
|---|---|---|---|
| 4.1 | Řezníček, T. (2023). *Open source and commercial CFD codes comparison on 2D airfoil analysis*. Bachelor's thesis, Brno University of Technology. https://theses.cz/id/nw0p3h/TomasReznicek_BP_2023.pdf | OpenFOAM and ANSYS Fluent compared with wind tunnel data | Checked |

| 4.2 | NASA Turbulence Modeling Resource. 2D NACA 0012 Airfoil Validation Case. https://tmbwg.github.io/turbmodels/naca0012_val.html | Wind tunnel data for the NACA 0012 at Re = 6 million. Ladson data with a trip strip: `NACA0012_validation/CLCD_Ladson_expdata.dat`. Abbott and von Doenhoff data without a trip strip: `NACA0012_validation/0012.abbottdata.cd.dat` | Checked. The two data files are in `data/naca0012/` and a test reads them |
| 4.3 | McCroskey, W. J. (1987). A Critical Assessment of Wind Tunnel Results for the NACA 0012 Airfoil. NASA TM 100019. https://ntrs.nasa.gov/citations/19880002254 | Review of NACA 0012 wind tunnel data | Checked (citation). Report not read |

| 4.4 | Ladson, C. L. (1988). NASA TM 4074. (Citation as given in the data file. Title not checked.) | Source of the tripped NACA 0012 data in item 4.2 | Checked (report number and year, from the data file) |

## 5. Removed items

- Thick-airfoil lift slope estimate, 2π(1 + 0.77 t). No source was identified,
  thus the test that used it was removed.

## 6. Items to complete

- Confirm item 1.4 (Anderson: edition and section) and the publication year of item 1.1.
- Resolve the paper number of item 2.3.
- Confirm the trip strip position in item 4.4.
- Confirm items 1.9, 1.12, 1.16 and the equations from memory in item 1.8.
- Find a source for the shape factor 1.4 at the start of turbulent flow.
