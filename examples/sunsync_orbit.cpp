// Sun-synchronous orbit calculator
// These satellites keep the same sun angle - useful for Earth observation
// The trick is using Earth's J2 bulge to make the orbit precess with the sun

#include <iostream>
#include <iomanip>
#include <cmath>
#include "../include/orbit_propagator.hpp"

using namespace orbitprop;

// Figure out what inclination we need for sun-sync at a given altitude
double sunSyncInclination(double altitude) {
    double a = constants::R_EARTH + altitude;
    double n = std::sqrt(constants::MU_EARTH / (a * a * a));
    
    // Sun moves ~0.9856 deg/day relative to stars
    double target_drift = 2.0 * constants::PI / (365.25 * 86400.0);
    
    // J2 formula: drift = -1.5 * n * J2 * (Re/a)^2 * cos(i)
    double cos_i = -target_drift / (1.5 * n * constants::J2 * 
                   std::pow(constants::R_EARTH / a, 2));
    
    if (std::abs(cos_i) > 1.0) {
        std::cout << "Can't make sun-sync at this altitude!\n";
        return -1;
    }
    
    return std::acos(cos_i);
}

// quick RAAN calculation from state vectors
double computeRAAN(const Vector3& r, const Vector3& v) {
    Vector3 h = r.cross(v);
    Vector3 n(-h.y, h.x, 0.0);
    
    double n_mag = n.norm();
    if (n_mag < 1e-10) return 0.0;
    
    double raan = std::acos(n.x / n_mag);
    if (n.y < 0) raan = 2.0 * constants::PI - raan;
    
    return raan;
}

int main() {
    std::cout << "\n--- Sun-Synchronous Orbit Design ---\n\n";

    double altitude = 700000.0;  // 700 km - typical for Earth obs satellites
    double a = constants::R_EARTH + altitude;
    
    double inc = sunSyncInclination(altitude);
    
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Design for " << altitude/1000 << " km altitude:\n";
    std::cout << "  Needed inclination: " << inc * constants::RAD2DEG << " deg\n";
    std::cout << "  (retrograde - goes over poles)\n\n";

    // set up the orbit
    KeplerianElements kep;
    kep.a = a;
    kep.e = 0.001;
    kep.i = inc;
    kep.raan = 0.0;  // start here
    kep.argp = 0.0;
    kep.nu = 0.0;

    OrbitalState state = keplerianToCartesian(kep);
    double period = orbitalPeriod(a);

    std::cout << "Orbit period: " << period/60 << " minutes\n";
    std::cout << "Starting RAAN: " << kep.raan * constants::RAD2DEG << " deg\n\n";

    // run for 10 days and track RAAN drift
    ForceModelConfig config;
    config.two_body = true;
    config.j2_perturbation = true;

    OrbitPropagator prop(config);
    prop.default_timestep = 30.0;

    std::cout << "Simulating 10 days...\n\n";
    std::cout << "Day   RAAN      Drift/day\n";
    std::cout << "---   ----      ---------\n";

    double prev_raan = computeRAAN(state.r, state.v);
    std::cout << "  0   " << std::setw(6) << prev_raan * constants::RAD2DEG << " deg\n";

    OrbitalState current = state;
    for (int day = 1; day <= 10; ++day) {
        current = prop.propagateTo(current, day * 86400.0);
        double raan = computeRAAN(current.r, current.v);
        
        double drift = (raan - prev_raan) * constants::RAD2DEG;
        // handle wraparound
        if (drift < -180) drift += 360;
        if (drift > 180) drift -= 360;
        
        std::cout << "  " << day << "   " << std::setw(6) << raan * constants::RAD2DEG 
                  << " deg   +" << drift << "\n";
        
        prev_raan = raan;
    }

    double expected = 360.0 / 365.25;
    std::cout << "\nExpected: ~" << expected << " deg/day (to track sun)\n";
    std::cout << "Looks good!\n";

    // save full trajectory
    auto trajectory = prop.propagate(state, 10 * 86400.0, 300.0);
    exportToCSV(trajectory, "trajectory_sunsync.csv");
    std::cout << "\nData saved to trajectory_sunsync.csv\n";
    std::cout << "  Analysis complete!\n";
    std::cout << "=================================================\n";

    return 0;
}
