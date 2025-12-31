#ifndef ORBIT_PROPAGATOR_HPP
#define ORBIT_PROPAGATOR_HPP

#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace orbitprop {

// WGS84 / EGM96 Constants
namespace constants {
    constexpr double MU_EARTH = 3.986004418e14;      // Earth gravitational parameter [m^3/s^2]
    constexpr double R_EARTH = 6378137.0;            // Earth equatorial radius [m]
    constexpr double J2 = 1.08262668e-3;             // J2 zonal harmonic coefficient
    constexpr double OMEGA_EARTH = 7.2921159e-5;     // Earth rotation rate [rad/s]
    constexpr double PI = 3.14159265358979323846;
    constexpr double DEG2RAD = PI / 180.0;
    constexpr double RAD2DEG = 180.0 / PI;
}

/**
 * @brief Simple 3D vector class for position/velocity operations
 */
class Vector3 {
public:
    double x, y, z;

    Vector3() : x(0.0), y(0.0), z(0.0) {}
    Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // Vector operations
    Vector3 operator+(const Vector3& v) const { return Vector3(x + v.x, y + v.y, z + v.z); }
    Vector3 operator-(const Vector3& v) const { return Vector3(x - v.x, y - v.y, z - v.z); }
    Vector3 operator*(double s) const { return Vector3(x * s, y * s, z * s); }
    Vector3 operator/(double s) const { return Vector3(x / s, y / s, z / s); }
    
    Vector3& operator+=(const Vector3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vector3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }

    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    double dot(const Vector3& v) const { return x*v.x + y*v.y + z*v.z; }
    Vector3 cross(const Vector3& v) const {
        return Vector3(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
    }
    Vector3 normalized() const { 
        double n = norm();
        if (n < 1e-15) return Vector3(0, 0, 0);
        return *this / n; 
    }
};

inline Vector3 operator*(double s, const Vector3& v) { return v * s; }

/**
 * @brief Classical orbital elements
 */
struct KeplerianElements {
    double a;       // Semi-major axis [m]
    double e;       // Eccentricity [-]
    double i;       // Inclination [rad]
    double raan;    // Right ascension of ascending node [rad]
    double argp;    // Argument of periapsis [rad]
    double nu;      // True anomaly [rad]
};

/**
 * @brief Spacecraft state in Cartesian coordinates (ECI frame)
 */
struct OrbitalState {
    double t;       // Time since epoch [s]
    Vector3 r;      // Position [m]
    Vector3 v;      // Velocity [m/s]

    OrbitalState() : t(0.0) {}
    OrbitalState(double t_, const Vector3& r_, const Vector3& v_) : t(t_), r(r_), v(v_) {}

    // Convert to state array for integration
    std::array<double, 6> toArray() const {
        return {r.x, r.y, r.z, v.x, v.y, v.z};
    }

    // Create from state array
    static OrbitalState fromArray(double t, const std::array<double, 6>& arr) {
        return OrbitalState(t, Vector3(arr[0], arr[1], arr[2]), Vector3(arr[3], arr[4], arr[5]));
    }
};

/**
 * @brief Configuration for force models
 */
struct ForceModelConfig {
    bool two_body = true;       // Central body gravity
    bool j2_perturbation = true; // J2 oblateness effect
    double mu = constants::MU_EARTH;
    double re = constants::R_EARTH;
    double j2 = constants::J2;
};

/**
 * @brief Computes gravitational acceleration with optional J2 perturbation
 */
class ForceModel {
public:
    ForceModelConfig config;

    ForceModel() = default;
    explicit ForceModel(const ForceModelConfig& cfg) : config(cfg) {}

    /**
     * @brief Compute total acceleration at given position
     * @param r Position vector in ECI [m]
     * @return Acceleration vector [m/s^2]
     */
    Vector3 computeAcceleration(const Vector3& r) const {
        Vector3 acc(0, 0, 0);

        if (config.two_body) {
            acc += twoBodyAcceleration(r);
        }

        if (config.j2_perturbation) {
            acc += j2Acceleration(r);
        }

        return acc;
    }

private:
    Vector3 twoBodyAcceleration(const Vector3& r) const {
        double r_mag = r.norm();
        if (r_mag < 1e-10) {
            throw std::runtime_error("Position magnitude too small");
        }
        return -config.mu / (r_mag * r_mag * r_mag) * r;
    }

    Vector3 j2Acceleration(const Vector3& r) const {
        double r_mag = r.norm();
        if (r_mag < 1e-10) return Vector3(0, 0, 0);

        double r2 = r_mag * r_mag;
        double r5 = r2 * r2 * r_mag;
        double re2 = config.re * config.re;
        double z2 = r.z * r.z;

        double factor = 1.5 * config.j2 * config.mu * re2 / r5;
        double z2_r2 = z2 / r2;

        double ax = factor * r.x * (5.0 * z2_r2 - 1.0);
        double ay = factor * r.y * (5.0 * z2_r2 - 1.0);
        double az = factor * r.z * (5.0 * z2_r2 - 3.0);

        return Vector3(ax, ay, az);
    }
};

/**
 * @brief Runge-Kutta 4th order integrator for orbit propagation
 */
class RK4Integrator {
public:
    /**
     * @brief Perform single RK4 step
     * @param state Current orbital state
     * @param dt Time step [s]
     * @param force_model Force model for acceleration computation
     * @return New orbital state after dt
     */
    static OrbitalState step(const OrbitalState& state, double dt, const ForceModel& force_model) {
        auto y = state.toArray();
        
        auto k1 = computeDerivative(y, force_model);
        auto y_temp = addScaled(y, k1, 0.5 * dt);
        
        auto k2 = computeDerivative(y_temp, force_model);
        y_temp = addScaled(y, k2, 0.5 * dt);
        
        auto k3 = computeDerivative(y_temp, force_model);
        y_temp = addScaled(y, k3, dt);
        
        auto k4 = computeDerivative(y_temp, force_model);

        // Combine: y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        std::array<double, 6> y_new;
        for (int i = 0; i < 6; ++i) {
            y_new[i] = y[i] + dt / 6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        }

        return OrbitalState::fromArray(state.t + dt, y_new);
    }

private:
    static std::array<double, 6> computeDerivative(const std::array<double, 6>& y, 
                                                    const ForceModel& force_model) {
        Vector3 r(y[0], y[1], y[2]);
        Vector3 v(y[3], y[4], y[5]);
        Vector3 a = force_model.computeAcceleration(r);
        return {v.x, v.y, v.z, a.x, a.y, a.z};
    }

    static std::array<double, 6> addScaled(const std::array<double, 6>& y,
                                           const std::array<double, 6>& dy,
                                           double scale) {
        std::array<double, 6> result;
        for (int i = 0; i < 6; ++i) {
            result[i] = y[i] + scale * dy[i];
        }
        return result;
    }
};

/**
 * @brief Main orbit propagator class
 */
class OrbitPropagator {
public:
    ForceModel force_model;
    double default_timestep = 60.0;  // Default 60 second timestep

    OrbitPropagator() = default;
    explicit OrbitPropagator(const ForceModelConfig& config) : force_model(config) {}

    /**
     * @brief Propagate orbit from initial state
     * @param initial_state Initial orbital state
     * @param duration Total propagation time [s]
     * @param output_step Time between output states [s] (0 = every integration step)
     * @return Vector of orbital states
     */
    std::vector<OrbitalState> propagate(const OrbitalState& initial_state,
                                        double duration,
                                        double output_step = 0.0) const {
        std::vector<OrbitalState> trajectory;
        trajectory.push_back(initial_state);

        double dt = default_timestep;
        double t_end = initial_state.t + duration;
        double next_output = initial_state.t + (output_step > 0 ? output_step : dt);

        OrbitalState current = initial_state;

        while (current.t < t_end - 1e-10) {
            // Adjust last step to hit end time exactly
            double step = std::min(dt, t_end - current.t);
            current = RK4Integrator::step(current, step, force_model);

            // Store at output intervals
            if (output_step <= 0 || current.t >= next_output - 1e-10) {
                trajectory.push_back(current);
                next_output += output_step;
            }
        }

        return trajectory;
    }

    /**
     * @brief Propagate to a specific time
     * @param initial_state Initial orbital state
     * @param target_time Target time [s]
     * @return Orbital state at target time
     */
    OrbitalState propagateTo(const OrbitalState& initial_state, double target_time) const {
        double duration = target_time - initial_state.t;
        if (duration < 0) {
            throw std::runtime_error("Target time is before initial state time");
        }
        
        OrbitalState current = initial_state;
        double dt = default_timestep;

        while (current.t < target_time - 1e-10) {
            double step = std::min(dt, target_time - current.t);
            current = RK4Integrator::step(current, step, force_model);
        }

        return current;
    }
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Convert Keplerian elements to Cartesian state
 */
inline OrbitalState keplerianToCartesian(const KeplerianElements& kep, 
                                          double mu = constants::MU_EARTH) {
    double a = kep.a, e = kep.e, i = kep.i;
    double raan = kep.raan, argp = kep.argp, nu = kep.nu;

    // Position and velocity in perifocal frame
    double p = a * (1.0 - e * e);
    double r_mag = p / (1.0 + e * std::cos(nu));

    double r_pqw_x = r_mag * std::cos(nu);
    double r_pqw_y = r_mag * std::sin(nu);

    double v_pqw_x = -std::sqrt(mu / p) * std::sin(nu);
    double v_pqw_y = std::sqrt(mu / p) * (e + std::cos(nu));

    // Rotation matrix elements
    double cos_raan = std::cos(raan), sin_raan = std::sin(raan);
    double cos_argp = std::cos(argp), sin_argp = std::sin(argp);
    double cos_i = std::cos(i), sin_i = std::sin(i);

    // Transform to ECI
    Vector3 r, v;
    
    r.x = (cos_raan*cos_argp - sin_raan*sin_argp*cos_i) * r_pqw_x +
          (-cos_raan*sin_argp - sin_raan*cos_argp*cos_i) * r_pqw_y;
    r.y = (sin_raan*cos_argp + cos_raan*sin_argp*cos_i) * r_pqw_x +
          (-sin_raan*sin_argp + cos_raan*cos_argp*cos_i) * r_pqw_y;
    r.z = (sin_argp*sin_i) * r_pqw_x + (cos_argp*sin_i) * r_pqw_y;

    v.x = (cos_raan*cos_argp - sin_raan*sin_argp*cos_i) * v_pqw_x +
          (-cos_raan*sin_argp - sin_raan*cos_argp*cos_i) * v_pqw_y;
    v.y = (sin_raan*cos_argp + cos_raan*sin_argp*cos_i) * v_pqw_x +
          (-sin_raan*sin_argp + cos_raan*cos_argp*cos_i) * v_pqw_y;
    v.z = (sin_argp*sin_i) * v_pqw_x + (cos_argp*sin_i) * v_pqw_y;

    return OrbitalState(0.0, r, v);
}

/**
 * @brief Compute orbital period from semi-major axis
 */
inline double orbitalPeriod(double a, double mu = constants::MU_EARTH) {
    return 2.0 * constants::PI * std::sqrt(a * a * a / mu);
}

/**
 * @brief Export trajectory to CSV file
 */
inline void exportToCSV(const std::vector<OrbitalState>& trajectory,
                        const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    file << std::setprecision(12);
    file << "t_sec,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,r_mag_m,v_mag_m_s\n";

    for (const auto& state : trajectory) {
        file << state.t << ","
             << state.r.x << "," << state.r.y << "," << state.r.z << ","
             << state.v.x << "," << state.v.y << "," << state.v.z << ","
             << state.r.norm() << "," << state.v.norm() << "\n";
    }

    file.close();
}

} // namespace orbitprop

#endif // ORBIT_PROPAGATOR_HPP
