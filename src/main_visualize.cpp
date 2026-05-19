// main_visualize.cpp -- terminal ASCII animation of the inspiral.
//
// Renders both black holes orbiting in real time using ANSI escape codes,
// leaves a trail showing their spiral path, and plots the gravitational-
// wave strain h+(t) as an ASCII chart at the end.
//
// Tested on Linux, macOS, Windows Terminal, and the VS Code integrated
// terminal. Will look broken on legacy Windows cmd.exe.

#include "binary_system.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace {

using namespace gw;

// Display dimensions in characters.
constexpr int COLS = 70;
constexpr int ROWS = 26;

// ANSI escape codes for terminal control.
constexpr const char* CLEAR_SCREEN = "\033[2J";
constexpr const char* CURSOR_HOME  = "\033[H";
constexpr const char* HIDE_CURSOR  = "\033[?25l";
constexpr const char* SHOW_CURSOR  = "\033[?25h";

// 2D character grid that maps physical (x, y) coordinates to character
// cells, accounting for the ~2:1 aspect ratio of terminal characters so
// circles render as circles.
class AsciiCanvas {
public:
    explicit AsciiCanvas(double half_width_phys)
        : grid_(ROWS, std::string(COLS, ' ')),
          half_width_(half_width_phys),
          // Terminal chars are ~2x taller than wide; correct so that a
          // physical circle looks like a circle on screen.
          half_height_(2.0 * half_width_phys * ROWS / COLS)
    {}

    // Plot character `c` at physical coordinate (x, y). Off-canvas
    // coordinates are silently dropped.
    void plot(double x, double y, char c) {
        int col = static_cast<int>(std::round(
            (x / half_width_) * (COLS / 2.0) + COLS / 2.0));
        int row = static_cast<int>(std::round(
            -(y / half_height_) * (ROWS / 2.0) + ROWS / 2.0));
        if (col >= 0 && col < COLS && row >= 0 && row < ROWS) {
            grid_[row][col] = c;
        }
    }

    void render(std::ostream& os) const {
        for (const auto& r : grid_) os << "  " << r << '\n';
    }

private:
    std::vector<std::string> grid_;
    double half_width_;
    double half_height_;
};

// ASCII chart of h+(t) printed once the simulation finishes.
void plot_strain_chart(const std::vector<double>& times,
                       const std::vector<double>& strain) {
    if (strain.empty()) return;

    constexpr int W = 70;
    constexpr int H = 14;

    // Sub-sample the strain timeseries down to W columns.
    std::vector<double> cols(W);
    for (int i = 0; i < W; ++i) {
        std::size_t idx = static_cast<std::size_t>(
            static_cast<double>(i) / W * strain.size());
        cols[i] = strain[std::min(idx, strain.size() - 1)];
    }

    double peak = 0.0;
    for (double s : cols) peak = std::max(peak, std::fabs(s));
    if (peak == 0.0) return;

    std::vector<std::string> chart(H, std::string(W, ' '));
    const int mid = H / 2;
    for (int c = 0; c < W; ++c) chart[mid][c] = '-';   // zero line

    for (int c = 0; c < W; ++c) {
        double n = cols[c] / peak;
        int row = mid - static_cast<int>(std::round(n * (H / 2 - 1)));
        row = std::max(0, std::min(H - 1, row));
        chart[row][c] = '*';
    }

    std::cout << "\n  Gravitational-wave strain h+(t)   (peak-normalized)\n";
    std::cout << "  " << std::string(W, '-') << '\n';
    for (const auto& row : chart) std::cout << "  " << row << '\n';
    std::cout << "  " << std::string(W, '-') << '\n';
    std::cout << "  t = 0 s"
              << std::string(W - 6 - 16, ' ')
              << "t = " << std::fixed << std::setprecision(4)
              << times.back() << " s\n";
}

} // namespace

int main() {
    using namespace gw;

    BinarySystem bbh(36.0, 29.0, 35.0);  // GW150914-like

    // Choose viewport so that both bodies stay on screen comfortably.
    // The body furthest from the center of mass is at (max(m1,m2)/M) * r.
    const double max_body_dist =
        std::max(bbh.m1(), bbh.m2()) / bbh.total_mass() * bbh.separation();
    const double half_view = max_body_dist * 1.6;
    AsciiCanvas trail(half_view);

    const double sim_dt          = 1.0 / 8192.0;   // 8192 Hz integrator step
    const int    steps_per_frame = 4;              // ~0.49 ms simulated per frame
    const int    frame_delay_ms  = 30;             // ~33 fps wall-clock

    std::vector<double> times, hp_samples;

    std::cout << HIDE_CURSOR << CLEAR_SCREEN;

    while (!bbh.merged()) {
        const double r   = bbh.separation();
        const double phi = bbh.phase();
        const double M   = bbh.total_mass();

        // Body positions relative to the center of mass:
        //   r1 = (m2/M) * r * (cos phi, sin phi)
        //   r2 = -(m1/M) * r * (cos phi, sin phi)
        const double x1 =  (bbh.m2() / M) * r * std::cos(phi);
        const double y1 =  (bbh.m2() / M) * r * std::sin(phi);
        const double x2 = -(bbh.m1() / M) * r * std::cos(phi);
        const double y2 = -(bbh.m1() / M) * r * std::sin(phi);

        // Lay down a permanent trail point for each body.
        trail.plot(x1, y1, '.');
        trail.plot(x2, y2, '.');

        // Compose this frame: trail + center-of-mass marker + current bodies.
        AsciiCanvas frame = trail;
        frame.plot(0.0, 0.0, '+');
        frame.plot(x1, y1, '@');   // body 1 (heavier)
        frame.plot(x2, y2, 'O');   // body 2 (lighter)

        // Repaint the screen by moving the cursor home (no clear -> no flicker).
        std::cout << CURSOR_HOME;
        std::cout << "  ====== Binary Black Hole Inspiral (GW150914-like) ======\n";
        frame.render(std::cout);
        std::cout << "  ========================================================\n";
        std::cout << "   t = "       << std::fixed << std::setprecision(4) << bbh.time()
                  << " s    f_gw = " << std::setprecision(1) << bbh.gw_frequency()
                  << " Hz    sep = " << std::setprecision(0) << bbh.separation()/1000.0
                  << " km        \n";
        std::cout << "   '@' "  << std::setprecision(0) << bbh.m1()/M_SUN << " M_sun"
                  << "   'O' "  << bbh.m2()/M_SUN << " M_sun"
                  << "   '+' COM   '.' trail        \n";
        std::cout.flush();

        // Record strain for the post-merger chart.
        double hp, hc;
        bbh.strain(410.0 * MEGAPARSEC, 0.0, hp, hc);
        times.push_back(bbh.time());
        hp_samples.push_back(hp);

        std::this_thread::sleep_for(std::chrono::milliseconds(frame_delay_ms));

        for (int i = 0; i < steps_per_frame && !bbh.merged(); ++i) {
            bbh.step(sim_dt);
        }
    }

    std::cout << SHOW_CURSOR;
    std::cout << "\n  *** MERGED at t = " << std::fixed << std::setprecision(4)
              << bbh.time() << " s  --  separation reached ISCO ***\n";

    plot_strain_chart(times, hp_samples);
    std::cout << "\n";
    return 0;
}
