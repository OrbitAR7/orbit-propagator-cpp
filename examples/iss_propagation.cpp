// ISS orbit propagation example
// Compares simple 2-body physics vs J2 perturbation effects
// Run this to see how Earth's oblateness affects the orbit over time

#include <iostream>
#include <iomanip>
#include "../include/orbit_propagator.hpp"

using namespace orbitprop;

int main() {
    std::cout << "\n--- ISS Orbit Propagation ---\n\n";

    // Setting up ISS orbital parameters
    // Got these from NASA's website - pretty cool stuff
    KeplerianElements iss_kep;
    iss_kep.a = 6781000.0;                    // about 408 km up
    iss_kep.e = 0.0001;                       // almost circular
    iss_kep.i = 51.6 * constants::DEG2RAD;    // classic ISS inclination
    iss_kep.raan = 0.0;                       
    iss_kep.argp = 0.0;                       
    iss_kep.nu = 0.0;                         // starting position

    OrbitalState initial_state = keplerianToCartesian(iss_kep);
    double period = orbitalPeriod(iss_kep.a);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Initial setup:\n";
    std::cout << "  Altitude: " << (iss_kep.a - constants::R_EARTH) / 1000.0 << " km\n";
    std::cout << "  Inclination: " << iss_kep.i * constants::RAD2DEG << " degrees\n";
    std::cout << "  Period: " << period / 60.0 << " minutes\n\n";

    std::cout << "Starting position (x,y,z): [" << initial_state.r.x/1000 << ", " 
              << initial_state.r.y/1000 << ", " 
              << initial_state.r.z/1000 << "] km\n\n";

    // let's propagate for one full day
    double duration = 86400.0;  // 1 day in seconds
    double output_step = 60.0;  // save every minute

    // First run - simple 2-body physics (idealized)
    std::cout << "Running simulation for 1 day...\n";

    ForceModelConfig simple_config;
    simple_config.two_body = true;
    simple_config.j2_perturbation = false;

    OrbitPropagator simple_prop(simple_config);
    simple_prop.default_timestep = 30.0;
    auto trajectory_2body = simple_prop.propagate(initial_state, duration, output_step);

    // Second run - adding J2 (Earth's bulge effect)
    ForceModelConfig j2_config;
    j2_config.two_body = true;
    j2_config.j2_perturbation = true;

    OrbitPropagator j2_prop(j2_config);
    j2_prop.default_timestep = 30.0;
    auto trajectory_j2 = j2_prop.propagate(initial_state, duration, output_step);

    // Compare the results
    const auto& final_simple = trajectory_2body.back();
    const auto& final_j2 = trajectory_j2.back();

    std::cout << "\n--- Results after 1 day ---\n";
    
    std::cout << "\nSimple model final position: [" << final_simple.r.x/1000 << ", " 
              << final_simple.r.y/1000 << ", " 
              << final_simple.r.z/1000 << "] km\n";

    std::cout << "With J2 final position: [" << final_j2.r.x/1000 << ", " 
              << final_j2.r.y/1000 << ", " 
              << final_j2.r.z/1000 << "] km\n";

    Vector3 difference = final_j2.r - final_simple.r;
    std::cout << "\nDifference: " << difference.norm()/1000 << " km";
    std::cout << " (that's pretty significant!)\n";

    // quick calculation of the drift rate
    double n = std::sqrt(constants::MU_EARTH / (iss_kep.a * iss_kep.a * iss_kep.a));
    double raan_dot = -1.5 * n * constants::J2 * std::pow(constants::R_EARTH / iss_kep.a, 2) 
                      * std::cos(iss_kep.i);
    double drift_per_day = raan_dot * 86400.0 * constants::RAD2DEG;
    
    std::cout << "RAAN drift: " << drift_per_day << " deg/day\n";
    std::cout << "(orbital plane rotates westward)\n";

    // save the data
    exportToCSV(trajectory_2body, "trajectory_2body.csv");
    exportToCSV(trajectory_j2, "trajectory_j2.csv");

    std::cout << "\nSaved " << trajectory_2body.size() << " data points to CSV files\n";
    std::cout << "Open the CSV files in Excel, Python, or MATLAB to visualize the orbits.\n";

    return 0;
}
