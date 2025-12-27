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

