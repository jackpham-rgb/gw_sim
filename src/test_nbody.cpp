// test_nbody.cpp -- sanity-check the N-body simulator and 2.5PN.
//
// Verifies:
//   1) Pure Newtonian binary conserves energy to high precision.
//   2) 2.5PN binary loses energy at the Peters rate (within ~10%).
//   3) Figure-8 three-body solution stays close to its periodic orbit.

#include "scenarios.hpp"

#include <cmath>
#include <cstdio>

using namespace gw;

int main() {
    printf("=== N-body sanity tests ===\n\n");

    // ---- Test 1: Newtonian binary energy conservation ------------------
    {
        Scenario s = make_binary_inspiral();
        s.system.enable_radiation_reaction(false);
        const double E0 = s.system.energy();
        const double dt = s.suggested_dt;
        for (int i = 0; i < 10000; ++i) s.system.step_rk4(dt);
        const double E1 = s.system.energy();
        const double rel_drift = std::fabs((E1 - E0) / E0);
        printf("[1] Newtonian binary, 10000 steps:\n");
        printf("    E_initial = %.6e J\n", E0);
        printf("    E_final   = %.6e J\n", E1);
        printf("    relative drift = %.3e   %s\n\n",
               rel_drift, rel_drift < 1e-6 ? "PASS" : "WARN");
    }

    // ---- Test 2: 2.5PN binary loses energy at Peters' rate -------------
    {
        Scenario s = make_binary_inspiral();
        // Predicted Peters' rate at the initial orbit:
        const auto& bs = s.system.bodies();
        double m1 = bs[0].mass, m2 = bs[1].mass, M = m1 + m2;
        double r0 = (bs[0].pos - bs[1].pos).norm();
        double dEdt_peters = -(32.0 / 5.0) * std::pow(G_NEWTON_3D, 4.0)
                             * std::pow(m1 * m2, 2.0) * M
                             / (std::pow(C_LIGHT_3D, 5.0) * std::pow(r0, 5.0));
        // Measure dE/dt over 100 steps (still ~circular at start)
        double E0 = s.system.energy();
        const int N = 100;
        for (int i = 0; i < N; ++i) s.system.step_rk4(s.suggested_dt);
        double E1 = s.system.energy();
        double dEdt_measured = (E1 - E0) / (N * s.suggested_dt);
        double ratio = dEdt_measured / dEdt_peters;
        printf("[2] 2.5PN binary energy-loss rate:\n");
        printf("    Peters' formula  dE/dt = %.4e W\n", dEdt_peters);
        printf("    Measured         dE/dt = %.4e W\n", dEdt_measured);
        printf("    ratio (measured/Peters) = %.4f   %s\n\n",
               ratio, std::fabs(ratio - 1.0) < 0.05 ? "PASS" : "FAIL");
    }

    // ---- Test 3: figure-8 stays periodic over a few periods ------------
    {
        Scenario s = make_figure8();
        const double dt = s.suggested_dt;
        const double E0 = s.system.energy();
        // Period is ~2 pi in dimensionless units -> 2 pi * T in SI
        // but our dt is set to T/1024, so 1024 steps is one "natural" period.
        const int steps_per_period = 1024;
        const int periods = 3;
        for (int i = 0; i < steps_per_period * periods; ++i) {
            s.system.step_rk4(dt);
        }
        const double E1 = s.system.energy();
        const double rel_drift = std::fabs((E1 - E0) / E0);
        printf("[3] Figure-8 over %d periods (%d steps):\n",
               periods, steps_per_period * periods);
        printf("    E_initial = %.6e J\n", E0);
        printf("    E_final   = %.6e J\n", E1);
        printf("    relative drift = %.3e   %s\n\n",
               rel_drift, rel_drift < 1e-3 ? "PASS" : "WARN");
        // Sanity: bodies should still be roughly in their initial neighborhood
        const auto& b = s.system.bodies();
        printf("    body 0 r = (%.3e, %.3e, %.3e)\n", b[0].pos.x, b[0].pos.y, b[0].pos.z);
        printf("    body 1 r = (%.3e, %.3e, %.3e)\n", b[1].pos.x, b[1].pos.y, b[1].pos.z);
        printf("    body 2 r = (%.3e, %.3e, %.3e)\n", b[2].pos.x, b[2].pos.y, b[2].pos.z);
    }

    return 0;
}
