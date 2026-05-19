// nbody.cpp -- implementation of the N-body simulator with 2.5PN.
//
// References:
//   Peters, P. C., Phys. Rev. 136, B1224 (1964).
//   Damour, T. and Deruelle, N., Phys. Lett. A 87, 81 (1981).
//   Iyer, B. R. and Will, C. M., Phys. Rev. D 52, 6882 (1995).
//   Blanchet, L., "Gravitational Radiation from Post-Newtonian Sources
//     and Inspiralling Compact Binaries", Living Rev. Relativity 17, 2
//     (2014).

#include "nbody.hpp"

#include <cmath>

namespace gw {

NBodySystem::NBodySystem() = default;

std::size_t NBodySystem::add_body(const Body& b) {
    bodies_.push_back(b);
    return bodies_.size() - 1;
}

Vec3 NBodySystem::newtonian_accel(std::size_t i,
                                  const std::vector<Body>& state) const
{
    Vec3 a{};
    const Vec3& ri = state[i].pos;
    for (std::size_t j = 0; j < state.size(); ++j) {
        if (j == i) continue;
        Vec3 dr  = state[j].pos - ri;
        double r2 = dr.norm2();
        if (r2 == 0.0) continue;
        double r  = std::sqrt(r2);
        a += dr * (G_NEWTON_3D * state[j].mass / (r2 * r));
    }
    return a;
}

Vec3 NBodySystem::pn25_accel(std::size_t i,
                             const std::vector<Body>& state) const
{
    if (!include_pn25_) return {};
    if (binary_i_ < 0 || binary_j_ < 0) return {};

    // Are we one of the inner-binary bodies?
    int partner = -1;
    if (static_cast<int>(i) == binary_i_) partner = binary_j_;
    else if (static_cast<int>(i) == binary_j_) partner = binary_i_;
    else return {};

    const Body& self = state[i];
    const Body& other = state[partner];

    // Relative coordinates of the binary, "self - other" so that the
    // acceleration acts to drain *self's* relative motion.
    const Vec3 r_vec = self.pos - other.pos;
    const Vec3 v_vec = self.vel - other.vel;
    const double r   = r_vec.norm();
    if (r == 0.0) return {};
    const Vec3 n_hat = r_vec / r;
    const double r_dot = v_vec.dot(n_hat);          // dr/dt (scalar)
    const double v2    = v_vec.norm2();

    const double M  = self.mass + other.mass;       // total mass
    const double mu = self.mass * other.mass / M;   // reduced mass
    const double nu = mu / M;                       // symmetric mass ratio

    // 2.5PN relative acceleration in harmonic coordinates
    // (Blanchet, Living Rev. Relativity 17 (2014), eq. 203):
    //
    //   a_rel = +(8/5) G^2 M^2 nu / (c^5 r^3) *
    //           {[3 v^2 + 17 G M / (3 r)] r_dot n_hat
    //            - [v^2 + 3 G M / r] v}
    //
    // For circular orbits (r_dot = 0, v^2 = G M / r) this reduces to
    //   a_rel = -(32/5) G^3 M^3 nu / (c^5 r^4) v,
    // which gives <dE/dt> = -(32/5) G^4 (m1 m2)^2 (m1+m2) / (c^5 r^5),
    // matching Peters (1964).
    const double GM      = G_NEWTON_3D * M;
    const double prefactor = (8.0 / 5.0) * GM * GM * nu /
                             (std::pow(C_LIGHT_3D, 5.0) * std::pow(r, 3.0));

    const Vec3 term_n = n_hat * (r_dot * (3.0 * v2 + 17.0 * GM / (3.0 * r)));
    const Vec3 term_v = v_vec * (v2 + 3.0 * GM / r);

    Vec3 a_rel = (term_n - term_v) * prefactor;

    // a_rel is the acceleration of the relative coordinate. It splits
    // between the two bodies as a_self = (m_other / M) * a_rel,
    // a_other = -(m_self / M) * a_rel  (Newton's 3rd law on the binary's
    // internal force). We return the contribution to *self*.
    return a_rel * (other.mass / M);
}

Vec3 NBodySystem::total_accel(std::size_t i,
                              const std::vector<Body>& state) const
{
    return newtonian_accel(i, state) + pn25_accel(i, state);
}

void NBodySystem::step_rk4(double dt) {
    const std::size_t N = bodies_.size();

    // RK4 on the 6N-dimensional state (positions + velocities).
    // We store derivatives as N pairs (v, a).
    auto deriv = [this, N](const std::vector<Body>& st) {
        std::vector<Body> d(N);
        for (std::size_t i = 0; i < N; ++i) {
            d[i].mass = st[i].mass;        // unused, but keep masses aligned
            d[i].pos  = st[i].vel;          // dr/dt = v
            d[i].vel  = total_accel(i, st); // dv/dt = a
        }
        return d;
    };

    // k1
    auto k1 = deriv(bodies_);

    // k2
    auto y2 = bodies_;
    for (std::size_t i = 0; i < N; ++i) {
        y2[i].pos += k1[i].pos * (0.5 * dt);
        y2[i].vel += k1[i].vel * (0.5 * dt);
    }
    auto k2 = deriv(y2);

    // k3
    auto y3 = bodies_;
    for (std::size_t i = 0; i < N; ++i) {
        y3[i].pos += k2[i].pos * (0.5 * dt);
        y3[i].vel += k2[i].vel * (0.5 * dt);
    }
    auto k3 = deriv(y3);

    // k4
    auto y4 = bodies_;
    for (std::size_t i = 0; i < N; ++i) {
        y4[i].pos += k3[i].pos * dt;
        y4[i].vel += k3[i].vel * dt;
    }
    auto k4 = deriv(y4);

    // Combine
    for (std::size_t i = 0; i < N; ++i) {
        bodies_[i].pos += (k1[i].pos + k2[i].pos * 2.0 + k3[i].pos * 2.0 + k4[i].pos)
                          * (dt / 6.0);
        bodies_[i].vel += (k1[i].vel + k2[i].vel * 2.0 + k3[i].vel * 2.0 + k4[i].vel)
                          * (dt / 6.0);
    }
    t_ += dt;
}

double NBodySystem::inner_separation() const {
    if (binary_i_ < 0 || binary_j_ < 0) return 0.0;
    return (bodies_[binary_i_].pos - bodies_[binary_j_].pos).norm();
}

double NBodySystem::inner_semi_major_axis() const {
    if (binary_i_ < 0 || binary_j_ < 0) return 0.0;
    const Body& a = bodies_[binary_i_];
    const Body& b = bodies_[binary_j_];
    Vec3 r = a.pos - b.pos;
    Vec3 v = a.vel - b.vel;
    double r_norm = r.norm();
    double v2     = v.norm2();
    double M      = a.mass + b.mass;
    // Specific orbital energy: eps = v^2/2 - GM/r
    // For bound orbit: eps = -GM / (2 a) -> a = -GM / (2 eps)
    double eps = 0.5 * v2 - G_NEWTON_3D * M / r_norm;
    if (eps >= 0.0) return r_norm;     // unbound; return current separation
    return -G_NEWTON_3D * M / (2.0 * eps);
}

double NBodySystem::inner_eccentricity() const {
    if (binary_i_ < 0 || binary_j_ < 0) return 0.0;
    const Body& a = bodies_[binary_i_];
    const Body& b = bodies_[binary_j_];
    Vec3 r = a.pos - b.pos;
    Vec3 v = a.vel - b.vel;
    double M = a.mass + b.mass;
    double GM = G_NEWTON_3D * M;
    Vec3 h = r.cross(v);             // specific angular momentum
    Vec3 e_vec = v.cross(h) / GM - r / r.norm();
    return e_vec.norm();
}

double NBodySystem::inner_orbital_omega() const {
    if (binary_i_ < 0 || binary_j_ < 0) return 0.0;
    double a = inner_semi_major_axis();
    if (a <= 0.0) return 0.0;
    double M = bodies_[binary_i_].mass + bodies_[binary_j_].mass;
    return std::sqrt(G_NEWTON_3D * M / (a * a * a));
}

double NBodySystem::inner_gw_frequency() const {
    return inner_orbital_omega() / M_PI;
}

double NBodySystem::energy() const {
    double T = 0.0, U = 0.0;
    for (std::size_t i = 0; i < bodies_.size(); ++i) {
        T += 0.5 * bodies_[i].mass * bodies_[i].vel.norm2();
        for (std::size_t j = i + 1; j < bodies_.size(); ++j) {
            double r = (bodies_[i].pos - bodies_[j].pos).norm();
            if (r > 0.0) {
                U -= G_NEWTON_3D * bodies_[i].mass * bodies_[j].mass / r;
            }
        }
    }
    return T + U;
}

} // namespace gw
