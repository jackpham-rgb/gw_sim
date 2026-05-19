// binary_system.hpp -- A two-body compact binary in quasi-circular orbit,
// shrinking under the back-reaction of gravitational-wave emission
// (Peters 1964 / quadrupole formula).
//
// State is parametrized by:
//   r   = orbital separation (m)
//   phi = orbital phase (rad)
// and evolved with RK4. The 3D positions of the two bodies are recovered
// from (r, phi) when needed.

#pragma once

#include "rk4.hpp"

namespace gw {

// ---- Physical constants (SI) ----
constexpr double G_NEWTON   = 6.67430e-11;       // m^3 kg^-1 s^-2
constexpr double C_LIGHT    = 2.99792458e8;      // m/s
constexpr double M_SUN      = 1.98892e30;        // kg
constexpr double PARSEC     = 3.0857e16;         // m
constexpr double MEGAPARSEC = 1.0e6 * PARSEC;    // m

class BinarySystem {
public:
    // m1_solar, m2_solar -- component masses in solar masses
    // f_gw_initial       -- starting GW frequency in Hz (sets initial separation)
    BinarySystem(double m1_solar, double m2_solar, double f_gw_initial);

    // Advance the system by one RK4 step of size dt seconds.
    void step(double dt);

    // ---- Accessors ----
    double time()              const { return t_; }
    double separation()        const { return r_; }    // m
    double phase()             const { return phi_; }  // rad
    double m1()                const { return m1_; }   // kg
    double m2()                const { return m2_; }   // kg
    double total_mass()        const { return m1_ + m2_; }
    double orbital_frequency() const;                  // rad/s, omega
    double gw_frequency()      const;                  // Hz, f_gw = 2 f_orb
    double chirp_mass()        const;                  // kg

    // ISCO of a Schwarzschild BH of mass M = m1+m2: r_ISCO = 6 G M / c^2.
    // Used as the cutoff where Newtonian + quadrupole stops being trustworthy.
    double schwarzschild_isco() const;

    // True once separation has fallen below ISCO -- end of the inspiral phase.
    bool merged() const { return r_ <= schwarzschild_isco(); }

    // Quadrupole strain at an observer at luminosity distance D, with the
    // orbital angular momentum tilted by inclination iota from the line of
    // sight (iota = 0 means face-on -> maximum amplitude).
    void strain(double distance,
                double inclination,
                double& h_plus,
                double& h_cross) const;

private:
    double m1_, m2_;   // kg
    double r_;         // separation (m)
    double phi_;       // orbital phase (rad)
    double t_;         // simulation time (s)

    // Right-hand side of (dr/dt, dphi/dt) for the RK4 integrator.
    State<2> derivatives(double t, const State<2>& y) const;
};

} // namespace gw
