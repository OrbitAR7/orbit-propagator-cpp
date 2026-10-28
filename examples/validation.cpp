// Validation checks for the propagator:
//  1. Two-body only: after an integer number of periods the orbit should return to its initial position.
//  2. J2: numerical RAAN drift over one day vs. the analytic secular rate
//     dΩ/dt = -1.5 n J2 (Re/p)^2 cos i  (for a 700 km sun-synchronous orbit, ≈ +0.9856 deg/day).
#include "orbit_propagator.hpp"
#include <cmath>
#include <cstdio>

using namespace orbitprop;

int main() {
    const double PI = constants::PI, mu = constants::MU_EARTH, Re = constants::R_EARTH, J2 = constants::J2;
    KeplerianElements k{Re + 700e3, 0.001, 98.19 * PI / 180.0, 0.3, 0.5, 0.0};
    OrbitalState s0 = keplerianToCartesian(k);
    s0.t = 0.0;
    const double T = orbitalPeriod(k.a);

    ForceModelConfig two_body;
    two_body.j2_perturbation = false;
    OrbitPropagator p2(two_body);
    p2.default_timestep = 10.0;
    OrbitalState s = p2.propagateTo(s0, 10.0 * T);
    std::printf("Two-body, 10 periods, dt = 10 s: position closure error = %.3f m\n", (s.r - s0.r).norm());

    OrbitPropagator pj{ForceModelConfig{}};
    pj.default_timestep = 10.0;
    const double tf = 86400.0;
    OrbitalState sj = pj.propagateTo(s0, tf);
    auto raan = [](const OrbitalState& x) { Vector3 h = x.r.cross(x.v); return std::atan2(h.x, -h.y); };
    double numeric = raan(sj) - raan(s0);
    while (numeric > PI) numeric -= 2 * PI;
    while (numeric < -PI) numeric += 2 * PI;
    const double n = std::sqrt(mu / std::pow(k.a, 3)), p = k.a * (1 - k.e * k.e);
    const double analytic = -1.5 * n * J2 * std::pow(Re / p, 2) * std::cos(k.i) * tf;
    std::printf("J2 RAAN drift over 1 day: numeric %.5f deg, analytic secular %.5f deg\n",
                numeric * 180 / PI, analytic * 180 / PI);
    return 0;
}
