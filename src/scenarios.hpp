// scenarios.hpp -- preset N-body initial conditions.
//
// Each scenario configures an NBodySystem with bodies, sets up the
// inner-binary indices for radiation reaction (when applicable), and
// returns a recommended integration time step and viewport scale for
// the visualizer.

#pragma once

#include "nbody.hpp"

#include <cmath>
#include <string>

namespace gw {

struct Scenario {
    std::string  name;
    std::string  blurb;            // one-liner shown in the HUD
    NBodySystem  system;
    double       suggested_dt;     // s
    double       view_scale;       // m -- distance from origin to fit on screen
    bool         pn25_on;
};

// ---- Helper: place a circular Keplerian binary at the origin -----------
//
// Returns (r1, v1, r2, v2) for two masses on a circular orbit of
// separation r in the xy-plane, centered at COM with COM velocity zero.
// The orbital plane normal is +z; orbit is counter-clockwise viewed from
// +z.
inline void make_circular_binary(double m1, double m2, double sep,
                                 Vec3& r1, Vec3& v1,
                                 Vec3& r2, Vec3& v2)
{
    const double M = m1 + m2;
    const double v_orb = std::sqrt(G_NEWTON_3D * M / sep);
    // Body 1 at +(m2/M)*sep along +x; Body 2 at -(m1/M)*sep along +x
    r1 = Vec3{ +(m2 / M) * sep, 0.0, 0.0 };
    r2 = Vec3{ -(m1 / M) * sep, 0.0, 0.0 };
    // Velocities perpendicular (along +/- y), circular orbit speeds:
    //   v1 = (m2/M) v_orb,  v2 = -(m1/M) v_orb
    v1 = Vec3{ 0.0, +(m2 / M) * v_orb, 0.0 };
    v2 = Vec3{ 0.0, -(m1 / M) * v_orb, 0.0 };
}

// ============================================================
// Scenario 1: pure binary inspiral (GW150914-like)
// ============================================================
inline Scenario make_binary_inspiral() {
    Scenario s;
    s.name  = "Binary Inspiral (GW150914-like)";
    s.blurb = "Two black holes spiral together, losing energy to gravitational waves.";
    const double m1 = 36.0 * M_SUN_3D;
    const double m2 = 29.0 * M_SUN_3D;
    const double f_gw = 35.0;
    const double M = m1 + m2;
    const double sep = std::cbrt(G_NEWTON_3D * M / std::pow(M_PI * f_gw, 2.0));

    Vec3 r1, v1, r2, v2;
    make_circular_binary(m1, m2, sep, r1, v1, r2, v2);
    auto i = s.system.add_body({ m1, r1, v1 });
    auto j = s.system.add_body({ m2, r2, v2 });
    s.system.set_inner_binary(static_cast<int>(i), static_cast<int>(j));
    s.system.enable_radiation_reaction(true);

    s.suggested_dt = 1.0 / 16384.0;        // ~60 us
    s.view_scale   = sep * 1.4;
    s.pn25_on      = true;
    return s;
}

// ============================================================
// Scenario 2: Hierarchical triple (Kozai-Lidov candidate)
// ============================================================
//
// Inner binary: two ~30 M_sun BHs in a tight orbit.
// Outer body: a third ~50 M_sun perturber on a wide, *inclined* orbit
// at ~6x the inner separation. The inclination drives Kozai-Lidov
// oscillations in the inner binary's eccentricity.
//
// Real astrophysical KL cycles take megayears -- way too slow to
// visualize alongside a sub-second GW inspiral. To keep both
// timescales visible, we shrink the outer separation (so the outer
// period is only ~10x the inner period) and run with PN OFF, i.e.
// purely Newtonian. You will see the inner binary's orbit precess and
// flatten as the third body pumps angular momentum in and out.
inline Scenario make_hierarchical_triple() {
    Scenario s;
    s.name  = "Hierarchical Triple (Kozai-Lidov)";
    s.blurb = "Inclined third body perturbs the inner binary's orbit.";
    const double m1 = 30.0 * M_SUN_3D;
    const double m2 = 30.0 * M_SUN_3D;
    const double m3 = 50.0 * M_SUN_3D;
    const double f_gw_inner = 30.0;
    const double M_in = m1 + m2;
    const double a_in = std::cbrt(G_NEWTON_3D * M_in /
                                  std::pow(M_PI * f_gw_inner, 2.0));
    const double a_out = 6.0 * a_in;     // outer semi-major axis
    const double M_total = M_in + m3;

    // Build inner binary first (in xy plane)
    Vec3 r1, v1, r2, v2;
    make_circular_binary(m1, m2, a_in, r1, v1, r2, v2);
    // Now displace the inner binary's COM and add COM velocity so the
    // two (inner-COM, outer body) system itself orbits its mutual COM.
    const double v_outer = std::sqrt(G_NEWTON_3D * M_total / a_out);
    // Inner COM at +x_offset, outer body at -x_offset along the outer
    // separation axis. Choose offsets so total COM is at origin.
    const double x_in_com  = +(m3      / M_total) * a_out;
    const double x_outer   = -(M_in    / M_total) * a_out;
    const double v_in_com  = +(m3      / M_total) * v_outer;
    const double v_out     = -(M_in    / M_total) * v_outer;

    // Inclination of outer orbit relative to inner-binary plane.
    // 60 deg gives strong Kozai-Lidov effects.
    const double inc = 60.0 * M_PI / 180.0;
    // Offset the inner binary's COM along +x, give it +y velocity.
    // Place the outer body in a plane tilted by 'inc' around the x-axis:
    //   outer position is along -x,
    //   outer velocity is in (y,z) plane, tilted by 'inc' from the
    //   inner-binary plane (xy).
    Vec3 com_offset    { x_in_com, 0.0, 0.0 };
    Vec3 com_velocity  { 0.0, v_in_com, 0.0 };
    Vec3 outer_pos     { x_outer, 0.0, 0.0 };
    Vec3 outer_vel     { 0.0, v_out * std::cos(inc), v_out * std::sin(inc) };

    // Apply the inner-binary COM offset/velocity:
    r1 += com_offset;  v1 += com_velocity;
    r2 += com_offset;  v2 += com_velocity;

    auto i = s.system.add_body({ m1, r1, v1 });
    auto j = s.system.add_body({ m2, r2, v2 });
    /*auto k =*/ s.system.add_body({ m3, outer_pos, outer_vel });

    s.system.set_inner_binary(static_cast<int>(i), static_cast<int>(j));
    s.system.enable_radiation_reaction(false);   // pure Newtonian here

    s.suggested_dt = 1.0 / 16384.0;
    s.view_scale   = a_out * 1.4;
    s.pn25_on      = false;
    return s;
}

// ============================================================
// Scenario 3: Chenciner-Montgomery figure-8 three-body orbit
// ============================================================
//
// An exact periodic solution to the Newtonian three-body problem
// discovered numerically by Moore (1993) and proved to exist by
// Chenciner & Montgomery (Annals of Math., 2000). Three equal
// masses chase each other around a figure-eight curve.
//
// The initial conditions below (in units G = m = 1) are the standard
// ones from the literature; we rescale into SI so the orbit covers
// roughly 1 AU and one period takes a few seconds of wall-clock time
// at the visualizer's sim-rate.
inline Scenario make_figure8() {
    Scenario s;
    s.name  = "Figure-8 Three-Body (Chenciner-Montgomery)";
    s.blurb = "Three equal masses chase each other on a figure-eight orbit.";

    // Dimensionless ICs in G = m = 1 units:
    //   r1 = -r2 = ( 0.97000436, -0.24308753, 0)
    //   r3 = (0, 0, 0)
    //   v3 = (-0.93240737, -0.86473146, 0)
    //   v1 = v2 = -v3 / 2
    const double L = 1.0e10;                  // length scale, m
    const double m = 1.0e30;                  // mass scale, kg (~0.5 M_sun)
    // Time scale from G m / L^3 = 1/T^2 -> T = sqrt(L^3 / (G m))
    const double T = std::sqrt(std::pow(L, 3.0) / (G_NEWTON_3D * m));
    const double V = L / T;

    Vec3 r1{  0.97000436 * L, -0.24308753 * L, 0.0 };
    Vec3 r2{ -0.97000436 * L,  0.24308753 * L, 0.0 };
    Vec3 r3{  0.0,             0.0,            0.0 };
    Vec3 v3{ -0.93240737 * V, -0.86473146 * V, 0.0 };
    Vec3 v1 = -v3 * 0.5;
    Vec3 v2 = -v3 * 0.5;

    s.system.add_body({ m, r1, v1 });
    s.system.add_body({ m, r2, v2 });
    s.system.add_body({ m, r3, v3 });
    s.system.set_inner_binary(-1, -1);
    s.system.enable_radiation_reaction(false);

    s.suggested_dt = T / 1024.0;     // ~1024 steps per period
    s.view_scale   = 1.5 * L;
    s.pn25_on      = false;
    return s;
}

} // namespace gw
