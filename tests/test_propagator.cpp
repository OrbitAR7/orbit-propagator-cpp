// Checks the propagator against results known in closed form.
// Run with ctest, or directly: ./build/test_propagator
#include "orbit_propagator.hpp"
#include <cmath>
#include <cstdio>

using namespace orbitprop;

static int failures = 0;
static void check(bool ok, const char* what) {
    std::printf("%s  %s\n", ok ? "pass" : "FAIL", what);
    if (!ok) ++failures;
}

static const double PI = constants::PI, MU = constants::MU_EARTH, RE = constants::R_EARTH, J2 = constants::J2;

// Specific energy including the J2 term of the potential. J2 is conservative,
// so this should stay constant when the J2 force is on.
static double energy(const OrbitalState& s, bool with_j2) {
    const double r = s.r.norm();
    double U = -MU / r;
    if (with_j2) {
        const double sz = s.r.z / r;
        U += MU / r * J2 * (RE / r) * (RE / r) * 0.5 * (3.0 * sz * sz - 1.0);
    }
    return 0.5 * s.v.dot(s.v) + U;
}

static double raan(const OrbitalState& s) {
    const Vector3 h = s.r.cross(s.v);
    return std::atan2(h.x, -h.y);
}

static double wrap(double a) {
    while (a > PI) a -= 2 * PI;
    while (a < -PI) a += 2 * PI;
    return a;
}

int main() {
    // 700 km sun-synchronous orbit, the same as examples/validation.cpp
    const KeplerianElements k{RE + 700e3, 0.001, 98.19 * PI / 180.0, 0.3, 0.5, 0.0};
    OrbitalState s0 = keplerianToCartesian(k);
    s0.t = 0.0;
    const double T = orbitalPeriod(k.a);

    // 1. Keplerian -> Cartesian: radius from the conic equation, energy -mu/2a
    {
        const KeplerianElements ke{RE + 2000e3, 0.2, 0.9, 1.1, 2.3, 0.7};
        const OrbitalState s = keplerianToCartesian(ke);
        const double p = ke.a * (1 - ke.e * ke.e);
        const double r_expected = p / (1 + ke.e * std::cos(ke.nu));
        const double dr = std::fabs(s.r.norm() - r_expected);
        const double dE = std::fabs(energy(s, false) + MU / (2 * ke.a)) / (MU / (2 * ke.a));
        const double dh = std::fabs(s.r.cross(s.v).norm() - std::sqrt(MU * p)) / std::sqrt(MU * p);
        const double di = std::fabs(std::acos(s.r.cross(s.v).z / s.r.cross(s.v).norm()) - ke.i);
        std::printf("conversion: radius %.2e m, energy %.2e, |h| %.2e, inclination %.2e rad\n", dr, dE, dh, di);
        check(dr < 1e-6 && dE < 1e-13 && dh < 1e-13 && di < 1e-12, "Keplerian to Cartesian matches the conic");
    }

    // 2. Two-body: back to the start after exactly 10 periods
    ForceModelConfig two_body;
    two_body.j2_perturbation = false;
    OrbitPropagator p2(two_body);
    p2.default_timestep = 10.0;
    const OrbitalState s10 = p2.propagateTo(s0, 10.0 * T);
    const double closure = (s10.r - s0.r).norm();
    std::printf("two-body closure after 10 periods, dt = 10 s: %.3f m\n", closure);
    check(closure < 1.0, "two-body orbit closes to within 1 m after 10 periods");

    // 3. Two-body: energy and angular momentum are conserved
    {
        const double dE = std::fabs(energy(s10, false) / energy(s0, false) - 1.0);
        const Vector3 h0 = s0.r.cross(s0.v), h1 = s10.r.cross(s10.v);
        const double dh = (h1 - h0).norm() / h0.norm();
        std::printf("two-body drift over 10 periods: energy %.2e, angular momentum %.2e\n", dE, dh);
        check(dE < 1e-9 && dh < 1e-9, "two-body energy and angular momentum conserved");
    }

    // 4. RK4 is fourth order: halving the step cuts the closure error by about 16.
    //    The steps divide the period exactly, so there is no shorter final step.
    //    Coarser steps overshoot (4.7 at 100 -> 200 steps per orbit), and finer
    //    ones run into round-off (2.9 at 1600 -> 3200), so use 800 -> 1600.
    {
        OrbitPropagator coarse(two_body), fine(two_body);
        coarse.default_timestep = T / 800.0;
        fine.default_timestep = T / 1600.0;
        const double e_coarse = (coarse.propagateTo(s0, 10.0 * T).r - s0.r).norm();
        const double e_fine = (fine.propagateTo(s0, 10.0 * T).r - s0.r).norm();
        const double order = std::log2(e_coarse / e_fine);
        std::printf("closure error: 800 steps/orbit %.2e m, 1600 steps/orbit %.2e m, observed order %.2f\n",
                    e_coarse, e_fine, order);
        check(std::fabs(order - 4.0) < 0.3, "RK4 converges at fourth order");
    }

    // 5. J2: energy with the J2 potential and the polar angular momentum h_z are
    //    conserved (the J2 field is symmetric about the z axis)
    OrbitPropagator pj{ForceModelConfig{}};
    pj.default_timestep = 10.0;
    const double day = 86400.0;
    const OrbitalState sd = pj.propagateTo(s0, day);
    {
        const double dE = std::fabs(energy(sd, true) / energy(s0, true) - 1.0);
        const double hz0 = s0.r.cross(s0.v).z, hz1 = sd.r.cross(sd.v).z;
        const double dhz = std::fabs(hz1 / hz0 - 1.0);
        const double dE_kepler = std::fabs(energy(sd, false) / energy(s0, false) - 1.0);
        std::printf("J2 drift over 1 day: energy %.2e (vs %.2e without the J2 term), h_z %.2e\n",
                    dE, dE_kepler, dhz);
        check(dE < 1e-9 && dhz < 1e-9, "J2 energy and polar angular momentum conserved");
    }

    // 6. J2 nodal drift over one day vs the secular rate -1.5 n J2 (Re/p)^2 cos i
    {
        const double numeric = wrap(raan(sd) - raan(s0));
        const double n = std::sqrt(MU / std::pow(k.a, 3)), p = k.a * (1 - k.e * k.e);
        const double analytic = -1.5 * n * J2 * std::pow(RE / p, 2) * std::cos(k.i) * day;
        std::printf("RAAN drift over 1 day: numeric %.5f deg, analytic %.5f deg\n",
                    numeric * 180 / PI, analytic * 180 / PI);
        check(std::fabs(numeric / analytic - 1.0) < 0.01, "J2 nodal drift within 1% of the secular rate");
    }

    std::printf("%s\n", failures ? "SOME CHECKS FAILED" : "all checks passed");
    return failures ? 1 : 0;
}
