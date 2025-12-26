# Numerical Methods in C++

C++ implementations of numerical differentiation, interpolation, numerical integration, and trajectory simulation.

## Build
cmake -S . -B build
cmake --build build

## Run
./build/q2
./build/q3
./build/q4 25.0

## Report
See `docs/Report.pdf`.

## Data
This code expects `tracking_data.dat` and writes `speed_data.dat` in the repository root.
If your files are in `data/`, create symlinks:

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
