# Orbit Propagator C++

A simple C++ library for propagating spacecraft orbits with J2 perturbations. It's header-only, so you can just drop it into your project and start using it.

## What it does

- Propagates orbits using **two-body dynamics** (basic Keplerian motion)
- Includes **J2 perturbation** to account for Earth's equatorial bulge
- Uses **RK4 integration** for decent accuracy without being too complicated
- Convert between **Keplerian elements** and Cartesian coordinates
- Export trajectories to **CSV files** for plotting in Excel, Python, MATLAB, etc.
- No external dependencies - just C++17 standard library

## Quick Example

```cpp
#include "orbit_propagator.hpp"

using namespace orbitprop;

int main() {
    // Set up an ISS-like orbit
    KeplerianElements kep;
    kep.a = 6781000.0;                    // about 408 km altitude
    kep.e = 0.0001;                       // nearly circular
    kep.i = 51.6 * constants::DEG2RAD;   // ISS inclination
    kep.raan = 0.0;                       
    kep.argp = 0.0;                       
    kep.nu = 0.0;                         

    OrbitalState initial = keplerianToCartesian(kep);

    // Turn on J2 perturbation
    ForceModelConfig config;
    config.two_body = true;
    config.j2_perturbation = true;

    OrbitPropagator propagator(config);
    propagator.default_timestep = 30.0;  // integrate every 30 seconds

    // Run for 1 day, save every minute
    auto trajectory = propagator.propagate(initial, 86400.0, 60.0);

    // Save to CSV
    exportToCSV(trajectory, "trajectory.csv");

    return 0;
}
```

## Building
How to Build

You'll need:
- A C++17 compiler (GCC, Clang, or MSVC)
- CMake 3.14 or newer

```bash
mkdir build && cd build
cmake ..
make

# Run the examples
./iss_propagation
./sunsync_orbit
```

Or just copy `include/orbit_propagator.hpp` into your project and include it directly - no build system needed.
## Examples

### 1. ISS-like Orbit Propagation (`examples/iss_propagation.cpp`)

Propagates an ISS-like orbit (408 km, 51.6° inclination) with and without J2, demonstrating:
- Position divergence due to J2 over time
- Theoretical vs. numerical RAAN drift rate

### 2. Sun-Synchronous Orbit Design (`examples/sunsync_orbit.cpp`)

Calculates the required inclination for a sun-synchronous orbit and verifies the RAAN precession rate matches the expected ~0.9856°/day.

## API Reference
ISS Orbit (`examples/iss_propagation.cpp`)

Simulates the ISS orbit for 1 day, comparing simple two-body physics vs. the more realistic J2 model. You'll see how much Earth's bulge affects the orbit - over 1000 km difference after just one day!

### Sun-Synchronous Orbit (`examples/sunsync_orbit.cpp`)

Designs a sun-synchronous orbit at 700 km altitude. These are the orbits that Earth observation satellites use to keep consistent lighting. The code calculates what inclination you need and verifies the drift rate.

Both examples output CSV files that you can plot in Excel, Python, MATLAB, or whatever you prefer
- `Vector3 r` — Position in ECI [m]
- `Vector3 v` — Velocity in ECI [m/s]

#### `KeplerianElements`
Classical orbital elements:
- `a` — Semi-major axis [m]
- `e` — Eccentricity [-]
- `i` — Inclination [rad]
- `raan` — Right ascension of ascending node [rad]
- `argp` — Argument of periapsis [rad]
- `nu` — True anomaly [rad]

#### `ForceModelConfig`
Configuration for force models:
- `bool two_body` — Enable central body gravity
- `bool j2_perturbation` — Enable J2 oblateness
- `double mu` — Gravitational parameter (default: Earth)
- `double re` — Equatorial radius (default: Earth)
- `double j2` — J2 coefficient (default: Earth)

#### `OrbitPropagator`
Main propagator class:
- `propagate(initial, duration, output_step)` — Returns trajectory vector
- `propagateTo(initial, target_time)` — Returns state at specific time
- `default_timestep` — Integration step size [s]

### Utility Functions

```cpp
// Convert Keplerian elements to Cartesian state
OrbitalState keplerianToCartesian(const KeplerianElements& kep, double mu);

// Compute orbital period from semi-major axis
double orbitalPeriod(double a, double mu);

// Export trajectory to CSV file
void exportToCSV(const std::vector<OrbitalState>& trajectory, const std::string& filename);
```

### Constants (`orbitprop::constants`)

| Constant | Value | Description |
|----------|-------|-------------|
| `MU_EARTH` | 3.986004418×10¹⁴ m³/s² | Earth gravitational parameter |
| `R_EARTH` | 6,378,137 m | Earth equatorial radius |
| `J2` | 1.08262668×10⁻³ | J2 zonal harmonic |
| `OMEGA_EARTH` | 7.2921159×10⁻⁵ rad/s | Earth rotation rate |

## J2 Perturbation

The J2 term accounts for Earth's equatorial bulge, causing:

1. **RAAN precession** — Orbital plane rotates about Earth's axis
   - Westward for prograde orbits (i < 90°)
   - Enables sun-synchronous orbits at specific inclinations

2. **Argument of perigee drift** — Apse line rotates within orbital plane

3. **Mean motion change** — Slight modification to orbital period
About J2 Perturbation

J2 is the biggest effect after simple two-body gravity. It accounts for Earth's equatorial bulge (Earth is slightly fatter at the equator). This causes:

1. **RAAN precession** - your orbital plane slowly rotates around Earth's axis
   - Goes westward for normal orbits (inclination < 90°)
   - This is actually useful! Sun-synchronous satellites use this to track the sun

2. **Argument of perigee drift** - the point where your orbit is closest to Earth rotates


This is good enough for most student projects, simulations, and understanding orbital mechanics. The RK4 integrator with 30-60 second steps works well for Low Earth Orbit satellites.

If you need higher accuracy (like for real mission planning), you'd want:
- Smaller timesteps
- Better integrators
- More perturbations (atmospheric drag, solar pressure, Moon/Sun gravity)

But for learning and most applications, this is solid.

## Possible Improvements

Things I might add later:
- Atmospheric drag model
- Solar radiation pressure
- Moon and Sun gravitational effects
- Converting back from Cartesian to Keplerian elements
- Adaptive timestep (automatically adjusts for accuracy)
## References

1. Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications* (4th ed.)
If you want to learn more about orbital mechanics:
- Vallado's "Fundamentals of Astrodynamics and Applications" - the classic textbook
- Montenbruck & Gill's "Satellite Orbits" - very practical
- NASA's orbital mechanics documentation

## License

MIT License - use it however you want, just don't blame me if something breaks.

