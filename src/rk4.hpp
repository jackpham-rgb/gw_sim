// rk4.hpp -- Generic 4th-order Runge-Kutta integrator for systems of ODEs.
//
// Solves dy/dt = f(t, y) where y is a fixed-size state vector.
// Used by BinarySystem to advance the inspiral equations of motion.

#pragma once

#include <array>
#include <cstddef>
#include <functional>

namespace gw {

template <std::size_t N>
using State = std::array<double, N>;

template <std::size_t N>
using Derivative = std::function<State<N>(double, const State<N>&)>;

// One classical RK4 step. Returns y(t + dt) given y(t) and f.
template <std::size_t N>
State<N> rk4_step(const Derivative<N>& f,
                  double t,
                  const State<N>& y,
                  double dt)
{
    State<N> k1 = f(t, y);

    State<N> tmp;
    for (std::size_t i = 0; i < N; ++i) tmp[i] = y[i] + 0.5 * dt * k1[i];
    State<N> k2 = f(t + 0.5 * dt, tmp);

    for (std::size_t i = 0; i < N; ++i) tmp[i] = y[i] + 0.5 * dt * k2[i];
    State<N> k3 = f(t + 0.5 * dt, tmp);

    for (std::size_t i = 0; i < N; ++i) tmp[i] = y[i] + dt * k3[i];
    State<N> k4 = f(t + dt, tmp);

    State<N> y_next;
    for (std::size_t i = 0; i < N; ++i) {
        y_next[i] = y[i] + (dt / 6.0) *
            (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    return y_next;
}

} // namespace gw
