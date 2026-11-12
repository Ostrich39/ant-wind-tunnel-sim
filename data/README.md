# Data

## naca0012

Wind tunnel data for the NACA 0012 airfoil at a Reynolds number of 6 million.
The files are from the NASA Turbulence Modeling Resource:
https://tmbwg.github.io/turbmodels/naca0012_val.html

| File | Content | Source given in the file |
|---|---|---|
| `CLCD_Ladson_expdata.dat` | Angle, lift, drag. Transition tripped. Three grit sizes. Mach 0.15 | Ladson, NASA TM 4074, 1988 |
| `0012.abbottdata.cd.dat` | Lift, drag. No trip strip. Digitized from a plot (approximate) | Abbott and von Doenhoff, *Theory of Wing Sections* |

Do not change these files. The test `crates/solver/tests/wind_tunnel.rs` reads them.
