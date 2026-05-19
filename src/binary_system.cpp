// binary_system.cpp -- implementation of the inspiral equations.
//
// References:
//   Peters, P. C., "Gravitational Radiation and the Motion of Two Point
//     Masses", Phys. Rev. 136, B1224 (1964).
//   Maggiore, M., "Gravitational Waves Vol. 1: Theory and Experiments",
//     Oxford University Press (2008), Ch. 4.

#include "binary_system.hpp"

#include <cmath>

namespace gw {

BinarySystem::BinarySystem(double m1_solar, double m2_solar, double f_gw_initial)
    : m1_(m1_solar * M_SUN),
      m2_(m2_solar * M_SUN),
      phi_(0.0),
      t_(0.0)
{
    // For a circular orbit f_gw = 2 f_orb -> omega = pi * f_gw.
    // Kepler's third law gives the corresponding separation:
    //   omega^2 = G M / r^3   =>   r = (G M / omega^2)^(1/3)
    const double omega = M_PI * f_gw_initial;
    const double M     = m1_ + m2_;
    r_ = std::cbrt(G_NEWTON * M / (omega * omega));
}

double BinarySystem::orbital_frequency() const {
    const double M = m1_ + m2_;
    return std::sqrt(G_NEWTON * M / (r_ * r_ * r_));
}

double BinarySystem::gw_frequency() const {
    // f_gw = 2 * f_orbital = omega / pi
    return orbital_frequency() / M_PI;
}

double BinarySystem::chirp_mass() const {
    // M_c = (m1*m2)^(3/5) / (m1+m2)^(1/5)
    const double M = m1_ + m2_;
    return std::pow(m1_ * m2_, 3.0 / 5.0) / std::pow(M, 1.0 / 5.0);
}

double BinarySystem::schwarzschild_isco() const {
    const double M = m1_ + m2_;
    return 6.0 * G_NEWTON * M / (C_LIGHT * C_LIGHT);
}

State<2> BinarySystem::derivatives(double /*t*/, const State<2>& y) const {
    const double r = y[0];
    // y[1] is phi -- not needed on the right-hand side, the system is
    // autonomous and isotropic.

    const double M = m1_ + m2_;

    // Peters (1964) quadrupole orbit-shrinkage rate for a circular binary:
    //
    //   dr/dt = - (64/5) * G^3 * m1 * m2 * (m1 + m2) / (c^5 * r^3)
    //
    const double dr_dt = -(64.0 / 5.0) *
        std::pow(G_NEWTON, 3.0) * m1_ * m2_ * M /
        (std::pow(C_LIGHT, 5.0) * r * r * r);

    // Orbital phase advances at the Keplerian angular frequency:
    //   dphi/dt = omega = sqrt(G M / r^3)
    const double dphi_dt = std::sqrt(G_NEWTON * M / (r * r * r));

    return { dr_dt, dphi_dt };
}

void BinarySystem::step(double dt) {
    State<2> y = { r_, phi_ };
    Derivative<2> f = [this](double t, const State<2>& s) {
        return this->derivatives(t, s);
    };
    State<2> y_next = rk4_step<2>(f, t_, y, dt);
    r_   = y_next[0];
    phi_ = y_next[1];
    t_  += dt;
}

void BinarySystem::strain(double distance,
                          double inclination,
                          double& h_plus,
                          double& h_cross) const
{
    // Restricted-PN / quadrupole strain for a circular binary
    // (Maggiore Vol. 1, eq. 4.30):
    //
    //   h+ = (4/D) (G M_c / c^2)^(5/3) (pi f_gw / c)^(2/3)
    //          * (1 + cos^2 iota)/2 * cos(2 phi)
    //   hx = (4/D) (G M_c / c^2)^(5/3) (pi f_gw / c)^(2/3)
    //          * cos(iota)         * sin(2 phi)
    //
    // Note: f_gw and phi here are the *gravitational-wave* frequency and
    //       phase, but since f_gw = 2 f_orb the GW phase is just 2 * phi_.

    const double Mc   = chirp_mass();
    const double f_gw = gw_frequency();

    const double prefactor = (4.0 / distance)
        * std::pow(G_NEWTON * Mc / (C_LIGHT * C_LIGHT), 5.0 / 3.0)
        * std::pow(M_PI * f_gw / C_LIGHT,             2.0 / 3.0);

    const double cos_i = std::cos(inclination);
    const double plus_factor  = 0.5 * (1.0 + cos_i * cos_i);
    const double cross_factor = cos_i;

    h_plus  = prefactor * plus_factor  * std::cos(2.0 * phi_);
    h_cross = prefactor * cross_factor * std::sin(2.0 * phi_);
}

} // namespace gw
