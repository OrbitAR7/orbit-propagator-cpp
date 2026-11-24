# orbit-propagator-cpp

[![CI](https://github.com/OrbitAR7/orbit-propagator-cpp/actions/workflows/CI.yml/badge.svg)](https://github.com/OrbitAR7/orbit-propagator-cpp/actions/workflows/CI.yml)

A header-only C++17 orbit propagator: two-body gravity plus the J2 zonal harmonic, integrated with fixed-step RK4. It includes Keplerian ↔ Cartesian conversion and CSV export. It is checked against analytic results (see [Validation](#validation)).

## Quick example

```cpp
#include "orbit_propagator.hpp"
using namespace orbitprop;

int main() {
    KeplerianElements kep;
    kep.a = 6781000.0;                    // m (~408 km altitude)
    kep.e = 0.0001;
    kep.i = 51.6 * constants::DEG2RAD;
    kep.raan = 0.0; kep.argp = 0.0; kep.nu = 0.0;

    OrbitalState initial = keplerianToCartesian(kep);

    ForceModelConfig config;              // two-body + J2 by default
    OrbitPropagator propagator(config);
    propagator.default_timestep = 30.0;   // s

    auto trajectory = propagator.propagate(initial, 86400.0, 60.0);  // 1 day, output every 60 s
    exportToCSV(trajectory, "trajectory.csv");
}
```

## Build

```bash
cmake -S . -B build && cmake --build build
./build/iss_propagation      # ISS-like orbit with and without J2
./build/sunsync_orbit        # sun-synchronous inclination design
./build/validation           # checks below
ctest --test-dir build --output-on-failure
```

Or copy `include/orbit_propagator.hpp` into your project. It has no dependencies beyond the standard library.

## Validation

`examples/validation.cpp` runs two checks on a 700 km, 98.19° sun-synchronous orbit with a 10 s step:

| Check | Result |
|---|---|
| Two-body only: position after exactly 10 orbital periods vs. initial position | 0.21 m closure error |
| J2: RAAN drift over one day vs. analytic secular rate −1.5 n J2 (Re/p)² cos i | 0.9874° numeric vs. 0.9859° analytic (sun-synchronous target 0.9856°/day) |

The small J2 difference is expected: the numerical result includes short-period J2 oscillations that the secular formula averages out. For an ISS-like orbit, `iss_propagation` shows about 1,100 km position difference after one day between two-body and two-body + J2, and a nodal drift of about −5.0°/day.

### Tests

`tests/test_propagator.cpp` runs with `ctest` and in CI on Ubuntu and macOS. It checks:

| Check | Measured | Pass if |
|---|---|---|
| Keplerian → Cartesian: radius, energy −μ/2a, angular momentum, inclination | relative errors ~1e-16 | < 1e-13 |
| Two-body closure after 10 periods (10 s step) | 0.21 m | < 1 m |
| Two-body energy and angular momentum drift over 10 periods | 2e-10 | < 1e-9 |
| RK4 order of convergence (800 vs 1600 steps per orbit) | 4.18 | 4 ± 0.3 |
| J2: energy (with the J2 potential) and h_z drift over 1 day | 3e-10 | < 1e-9 |
| J2: RAAN drift over 1 day vs secular rate | 0.15% | < 1% |

The J2 energy check uses the full potential. Without the J2 term the energy changes by 1e-3 over the day, so the check would catch a sign or factor error in the J2 acceleration. ctest also runs the three examples to make sure they finish.

## Examples

- **`iss_propagation.cpp`**: propagates an ISS-like orbit (408 km, 51.6°) with and without J2 and compares the numerical RAAN drift with theory.
- **`sunsync_orbit.cpp`**: solves for the sun-synchronous inclination at 700 km and checks the precession rate.
- **`validation.cpp`**: the checks above.

Each writes CSV output for plotting.

## API summary

| Item | Description |
|---|---|
| `KeplerianElements` | `a` [m], `e`, `i`, `raan`, `argp`, `nu` [rad] |
| `OrbitalState` | time `t` [s], ECI position `r` [m] and velocity `v` [m/s] |
| `ForceModelConfig` | `two_body`, `j2_perturbation`, `mu`, `re`, `j2` |
| `OrbitPropagator` | `propagate(initial, duration, output_step)`, `propagateTo(initial, t)`, `default_timestep` |
| Utilities | `keplerianToCartesian`, `orbitalPeriod`, `exportToCSV` |

Constants: μ = 3.986004418 × 10¹⁴ m³/s², Rₑ = 6 378 137 m, J2 = 1.08262668 × 10⁻³.

## Limitations and next steps

- Fixed-step RK4 with J2 only: there is no drag, solar radiation pressure or third-body gravity, and no higher-order harmonics.
- No Cartesian → Keplerian conversion yet, and no adaptive step-size control.
- Planned: an exponential-atmosphere drag model, an adaptive Runge–Kutta (e.g., Dormand–Prince 5(4)), and osculating-element output.

## References

- Vallado, *Fundamentals of Astrodynamics and Applications*, 4th ed., 2013
- Montenbruck & Gill, *Satellite Orbits: Models, Methods and Applications*, Springer, 2000

## License

MIT. See `LICENSE`.
