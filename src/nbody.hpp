// nbody.hpp -- 3D N-body gravitational simulator.
//
// Newtonian gravity for all body pairs, with optional 2.5 post-Newtonian
// radiation reaction force on a designated "inner binary." The 2.5PN
// term uses the Iyer-Will / Damour-Deruelle form for the relative
// acceleration in harmonic coordinates, which reproduces Peters' (1964)
// orbit-averaged energy-loss rate <dE/dt> = -(32/5) G^4 (m1 m2)^2
// (m1+m2) / (c^5 r^5) for circular orbits.
//
// Integration is classical RK4. Good enough for visualization; for
// production work you'd want a symplectic integrator (e.g., leapfrog or
// Wisdom-Holman) so secular angular momentum and energy errors don't
// accumulate.

#pragma once

#include "vec3.hpp"

#include <cstddef>
#include <vector>

namespace gw {

// Physical constants (SI). Duplicated from binary_system.hpp deliberately
// so this module can stand alone.
constexpr double G_NEWTON_3D    = 6.67430e-11;       // m^3 kg^-1 s^-2
constexpr double C_LIGHT_3D     = 2.99792458e8;      // m/s
constexpr double M_SUN_3D       = 1.98892e30;        // kg
constexpr double PARSEC_3D      = 3.0857e16;         // m
constexpr double MEGAPARSEC_3D  = 1.0e6 * PARSEC_3D; // m

struct Body {
    double mass;   // kg
    Vec3   pos;    // m
    Vec3   vel;    // m/s
};

class NBodySystem {
public:
    NBodySystem();

    // Returns the index of the new body.
    std::size_t add_body(const Body& b);

    // Designate which two bodies form the "inner binary." The 2.5PN
    // radiation reaction force will be applied between these two when
    // enabled. Pass -1 to disable.
    void set_inner_binary(int i, int j) {
        binary_i_ = i;
        binary_j_ = j;
    }
    void enable_radiation_reaction(bool on) { include_pn25_ = on; }

    // Advance the system by dt seconds using one classical RK4 step.
    void step_rk4(double dt);

    // ---- Accessors ----
    double                    time()   const { return t_; }
    const std::vector<Body>&  bodies() const { return bodies_; }
    std::size_t               size()   const { return bodies_.size(); }
    bool                      uses_pn25() const { return include_pn25_; }

    // Instantaneous orbital observables for the inner binary. If no inner
    // binary is set, these return 0. They use the bound-Keplerian
    // formulas from r_rel and v_rel.
    double inner_separation()       const;  // m
    double inner_semi_major_axis()  const;  // m  (assumes bound)
    double inner_eccentricity()     const;
    double inner_orbital_omega()    const;  // rad/s
    double inner_gw_frequency()     const;  // Hz, = orbital_omega / pi

    // Total energy (kinetic + Newtonian potential) in the simulation.
    // Useful as a diagnostic: with PN off it should be conserved by RK4
    // to a few parts in 1e8 over a few orbits.
    double energy() const;

private:
    std::vector<Body> bodies_;
    double            t_           = 0.0;
    int               binary_i_    = -1;
    int               binary_j_    = -1;
    bool              include_pn25_ = false;

    // Newtonian acceleration of body i given a state snapshot.
    Vec3 newtonian_accel(std::size_t i,
                         const std::vector<Body>& state) const;

    // 2.5PN acceleration on body i due to the inner-binary partner.
    // Returns zero if i is not part of the inner binary or PN is off.
    Vec3 pn25_accel(std::size_t i,
                    const std::vector<Body>& state) const;

    // Full acceleration of body i (Newtonian + optional 2.5PN).
    Vec3 total_accel(std::size_t i,
                     const std::vector<Body>& state) const;
};

} // namespace gw
