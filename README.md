# Numerical Methods in C++

C++ implementations of core numerical methods applied to discrete tracking data and 3D trajectory simulation:
numerical differentiation (finite differences), interpolation (Newton divided differences), numerical integration (quadrature), and ODE integration (RK4).

📄 Full write-up: `docs/Report.pdf`

---

## Quickstart

```bash
cmake -S . -B build
cmake --build build

./build/q2
./build/q3
./build/q4 25.0
```

## What each executable does

* q2 — Kinematics from tracking data
Reads position samples and estimates velocity/acceleration using second-order centred finite differences, producing (e.g.) speed vs time outputs.

* q3 — Power interpolation + energy estimation
Interpolates a power–speed relationship using Newton divided differences and estimates total energy via composite quadrature (trapezoid / Simpson-type).

* q4 — Trajectory simulation
Simulates a 3D trajectory by integrating a first-order ODE system with RK4. Optional drag / Magnus terms are enabled via constants/flags in the code.

## Project structure
```bash
.
├── CMakeLists.txt
├── include/          # headers (constants + interpolation utilities)
├── src/              # q2/q3/q4 implementations
├── data/             # input dataset(s) (and optional stored outputs)
└── docs/             # report (PDF)
```

## Data and outputs
By default, the code expects the input/output files in the repository root:

* input: tracking_data.dat

* output: speed_data.dat

All outputs are plain-text files suitable for plotting/diagnostics.

## If your files live in data/ (symlink option)
If you keep inputs/outputs under data/, create symlinks so the executables can find them:

```bash
# from the repository root
ln -s data/tracking_data.dat tracking_data.dat

# generate speed_data.dat in the repo root
./build/q2

# (optional) store it under data/
mv speed_data.dat data/speed_data.dat

# make q3 still find it (since q3 reads speed_data.dat from repo root)
ln -s data/speed_data.dat speed_data.dat

# now run q3
./build/q3
```
## Build notes
Clean rebuild:
```bash
rm -rf build
cmake -S . -B build
cmake --build build
```
## Verification (expected results)
- `q2`: generates `speed_data.dat` (≈999 rows if using centred differences)
- Peak speed (from report): ~6.83 m/s
- `q3`: energy estimates ~111529 J (trapezoid) and ~111532 J (Simpson-type)


## Example output
<img width="1280" height="720" alt="speed_plot" src="https://github.com/user-attachments/assets/25a85d8e-501d-4d9b-85aa-97b3ffb2469c" />
