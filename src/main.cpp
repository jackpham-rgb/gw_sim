// main.cpp -- run a binary-black-hole inspiral and write out the waveform.
//
// Default parameters are GW150914-like (the first LIGO detection, 2015).
// Change them below to match GW250114 or any other event you care about.

#include "binary_system.hpp"
#include "wav_writer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

int main() {
    using namespace gw;

    // ---- Source parameters --------------------------------------------------
    // GW150914 (first direct GW detection) -- LIGO/Virgo inferred:
    //   m1 ~ 36 M_sun, m2 ~ 29 M_sun, D_L ~ 410 Mpc.
    // For GW250114 you'd plug in the values from the LVK paper (the masses
    // are similar; check arXiv:2601.13173 for the most recent fit).
    const double m1_solar    = 36.0;
    const double m2_solar    = 29.0;
    const double f_gw_start  = 35.0;                 // Hz, lower edge of LIGO band
    const double distance    = 410.0 * MEGAPARSEC;   // m
    const double inclination = 0.0;                  // rad, 0 = face-on

    // ---- Integrator settings ------------------------------------------------
    const int    sample_rate = 8192;                 // Hz (a real LIGO data rate)
    const double dt          = 1.0 / sample_rate;    // s
    const int    max_steps   = 60 * sample_rate;     // hard cap: 60 s

    BinarySystem bbh(m1_solar, m2_solar, f_gw_start);

    std::cout << "==== Binary black hole inspiral ====\n";
    std::cout << "  m1            = " << m1_solar << " M_sun\n";
    std::cout << "  m2            = " << m2_solar << " M_sun\n";
    std::cout << "  chirp mass    = " << bbh.chirp_mass() / M_SUN << " M_sun\n";
    std::cout << "  start f_gw    = " << f_gw_start << " Hz\n";
    std::cout << "  start sep.    = " << bbh.separation() / 1000.0 << " km\n";
    std::cout << "  ISCO sep.     = " << bbh.schwarzschild_isco() / 1000.0 << " km\n";
    std::cout << "  distance      = " << distance / MEGAPARSEC << " Mpc\n";
    std::cout << "  sample rate   = " << sample_rate << " Hz\n";
    std::cout << "  time step     = " << dt << " s\n\n";

    // ---- Integrate the inspiral --------------------------------------------
    std::vector<double> hp_samples;
    std::vector<double> hc_samples;
    hp_samples.reserve(max_steps);
    hc_samples.reserve(max_steps);

    std::ofstream csv("output/waveform.csv");
    if (!csv) {
        std::cerr << "Could not open output/waveform.csv for writing\n";
        return EXIT_FAILURE;
    }
    csv << "t_s,h_plus,h_cross,f_gw_Hz,separation_km\n";
    csv.precision(10);

    int steps = 0;
    while (!bbh.merged() && steps < max_steps) {
        double hp, hc;
        bbh.strain(distance, inclination, hp, hc);

        hp_samples.push_back(hp);
        hc_samples.push_back(hc);

        csv << bbh.time()                << ','
            << hp                        << ','
            << hc                        << ','
            << bbh.gw_frequency()        << ','
            << bbh.separation() / 1000.0 << '\n';

        bbh.step(dt);
        ++steps;
    }
    csv.close();

    std::cout << "Inspiral finished after " << steps << " steps ("
              << bbh.time() << " s)\n";
    std::cout << "  final f_gw    = " << bbh.gw_frequency() << " Hz\n";
    std::cout << "  final sep.    = " << bbh.separation() / 1000.0 << " km\n";
    std::cout << "  reached ISCO? = " << (bbh.merged() ? "yes" : "no") << "\n\n";

    // ---- Write a WAV of the chirp ------------------------------------------
    // Strain values are ~1e-21 -- way below audible. We just normalize to
    // peak amplitude 1.0 so the chirp comes out clearly.
    double peak = 0.0;
    for (double v : hp_samples) peak = std::max(peak, std::fabs(v));

    if (peak > 0.0) {
        std::vector<double> audio(hp_samples.size());
        for (std::size_t i = 0; i < hp_samples.size(); ++i) {
            audio[i] = hp_samples[i] / peak;
        }
        if (write_wav("output/chirp.wav", audio, sample_rate)) {
            std::cout << "Wrote output/chirp.wav  ("
                      << audio.size() << " samples @ "
                      << sample_rate << " Hz)\n";
        } else {
            std::cerr << "Failed to write output/chirp.wav\n";
        }
    }

    std::cout << "Wrote output/waveform.csv\n";
    return EXIT_SUCCESS;
}
