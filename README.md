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
