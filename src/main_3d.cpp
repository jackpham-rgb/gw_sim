// main_3d.cpp -- 3D real-time visualization of the N-body simulation
// using raylib. Switch between scenarios with number keys; orbit the
// camera with the left mouse button; zoom with the scroll wheel.
//
// Build: requires raylib (sudo apt install libraylib-dev on Ubuntu).
//   make 3d
// Run:
//   ./gw_3d

#include "nbody.hpp"
#include "scenarios.hpp"

#include <raylib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <deque>
#include <string>
#include <vector>

using namespace gw;

namespace {

// ----- Display constants -----
constexpr int   WINDOW_W       = 1280;
constexpr int   WINDOW_H       = 720;
constexpr float SCENE_RADIUS   = 5.0f;     // bodies normalized to ~|x|<5 in 3D world
constexpr int   MAX_TRAIL      = 800;
constexpr int   STRAIN_SAMPLES = 600;

// Per-body display colors (cycles if more than 4 bodies).
constexpr std::array<Color, 4> BODY_COLORS = {{
    Color{255, 100, 100, 255},  // red
    Color{100, 180, 255, 255},  // blue
    Color{255, 215,  90, 255},  // gold
    Color{160, 220, 130, 255},  // green
}};

// Convert a Vec3 in physical units to raylib Vector3 in scene units.
inline Vector3 to_scene(const Vec3& p, double view_scale) {
    return Vector3{
        static_cast<float>(p.x / view_scale * SCENE_RADIUS),
        static_cast<float>(p.z / view_scale * SCENE_RADIUS),  // z up in 3D
        static_cast<float>(p.y / view_scale * SCENE_RADIUS),
    };
}

// Simple manual orbital camera.
struct OrbitCamera {
    float theta  = 0.7f;        // azimuth (rad)
    float phi    = 0.4f;        // elevation (rad)
    float radius = 14.0f;
    Vector3 target = {0, 0, 0};

    Camera3D camera() const {
        Camera3D c;
        c.position = Vector3{
            target.x + radius * std::cos(phi) * std::sin(theta),
            target.y + radius * std::sin(phi),
            target.z + radius * std::cos(phi) * std::cos(theta),
        };
        c.target = target;
        c.up = Vector3{0, 1, 0};
        c.fovy = 50.0f;
        c.projection = CAMERA_PERSPECTIVE;
        return c;
    }

    void update() {
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            Vector2 d = GetMouseDelta();
            theta -= d.x * 0.006f;
            phi   -= d.y * 0.006f;
            phi = std::max(-1.5f, std::min(1.5f, phi));
        }
        radius -= GetMouseWheelMove() * 1.0f;
        radius = std::max(2.0f, std::min(80.0f, radius));
    }
};

// Per-body trail: a ring buffer of recent scene-space positions.
struct Trail {
    std::deque<Vector3> points;
    void push(Vector3 p) {
        points.push_back(p);
        if (points.size() > MAX_TRAIL) points.pop_front();
    }
    void clear() { points.clear(); }
    void draw(Color c) const {
        if (points.size() < 2) return;
        const std::size_t n = points.size();
        for (std::size_t i = 1; i < n; ++i) {
            float t = static_cast<float>(i) / static_cast<float>(n);
            Color cc = c;
            cc.a = static_cast<unsigned char>(255.0f * t * t);  // quadratic fade
            DrawLine3D(points[i - 1], points[i], cc);
        }
    }
};

// Approximate gravitational-wave strain h+ at distance D for the inner
// binary, using the quasi-circular quadrupole formula:
//   h+(t) = -(4 G^2 M mu / (c^4 D r)) * cos(2 phi)
// with phi the projected orbital angle in the binary's instantaneous
// plane. This is an approximation but tracks the chirp shape correctly.
double inner_strain_h_plus(const NBodySystem& sys, int i, int j) {
    if (i < 0 || j < 0) return 0.0;
    const Body& a = sys.bodies()[i];
    const Body& b = sys.bodies()[j];
    Vec3 r_rel = a.pos - b.pos;
    double r = r_rel.norm();
    if (r == 0.0) return 0.0;
    double phi = std::atan2(r_rel.y, r_rel.x);
    double M  = a.mass + b.mass;
    double mu = a.mass * b.mass / M;
    constexpr double D = 410.0 * MEGAPARSEC_3D;
    return -4.0 * G_NEWTON_3D * G_NEWTON_3D * M * mu /
           (std::pow(C_LIGHT_3D, 4.0) * D * r) *
           std::cos(2.0 * phi);
}

// Body display radius in scene units. We use a logarithmic mass scale
// so very different masses are still visible.
float body_radius(double mass_kg) {
    double m_solar = mass_kg / M_SUN_3D;
    return 0.10f + 0.05f * static_cast<float>(std::log10(1.0 + m_solar));
}

// One running visualization state.
struct Sim {
    Scenario              scenario;
    std::vector<Trail>    trails;
    std::deque<double>    strain_history;
    std::deque<double>    strain_times;
    int                   steps_per_frame = 5;
    bool                  paused = false;

    void load(Scenario s) {
        scenario = std::move(s);
        trails.assign(scenario.system.size(), Trail{});
        strain_history.clear();
        strain_times.clear();
        // Tune speed so a "natural" period takes ~5 s of wall clock.
        if (scenario.name.find("Figure-8") != std::string::npos) {
            steps_per_frame = 4;
        } else if (scenario.name.find("Triple") != std::string::npos) {
            steps_per_frame = 8;
        } else {
            steps_per_frame = 5;
        }
    }

    void step() {
        if (paused) return;
        for (int i = 0; i < steps_per_frame; ++i) {
            scenario.system.step_rk4(scenario.suggested_dt);
        }
        // Record strain from the inner binary (if any).
        // The strain itself is ~1e-21 -- we just store the raw value
        // and renormalize for plotting.
        // We take advantage of binary indices being stored on the system,
        // but they're private; just look at scenario type:
        if (scenario.name.find("Inspiral") != std::string::npos) {
            // bodies 0 and 1 are the inner binary
            double h = inner_strain_h_plus(scenario.system, 0, 1);
            strain_history.push_back(h);
            strain_times.push_back(scenario.system.time());
            if (strain_history.size() > STRAIN_SAMPLES) {
                strain_history.pop_front();
                strain_times.pop_front();
            }
        }
        // Update trails
        for (std::size_t b = 0; b < scenario.system.bodies().size(); ++b) {
            trails[b].push(to_scene(scenario.system.bodies()[b].pos,
                                    scenario.view_scale));
        }
    }

    void draw_3d() const {
        // Reference grid at z=0 (which is "y=0" in scene because we
        // mapped physical-z to scene-y above).
        DrawGrid(20, 0.5f);

        // Bodies + trails
        for (std::size_t b = 0; b < scenario.system.bodies().size(); ++b) {
            const Body& body = scenario.system.bodies()[b];
            Vector3 p = to_scene(body.pos, scenario.view_scale);
            Color c = BODY_COLORS[b % BODY_COLORS.size()];
            float r = body_radius(body.mass);
            // Draw trail behind the body
            trails[b].draw(c);
            // Body itself: draw a slightly darker outer "horizon" sphere
            // and a brighter core, so it reads as a black hole-ish thing.
            Color core = c;
            DrawSphereEx(p, r, 12, 12, core);
            // A faint "halo" ring
            DrawSphereWires(p, r * 1.3f, 6, 6,
                            Color{c.r, c.g, c.b, (unsigned char)80});
        }
    }

    void draw_hud() const {
        const int x = 14;
        int y = 12;
        const int line = 22;

        DrawRectangle(0, 0, WINDOW_W, 4 * line + 8,
                      Color{0, 0, 0, 160});

        DrawText(scenario.name.c_str(), x, y, 20, RAYWHITE); y += line;
        DrawText(scenario.blurb.c_str(), x, y, 16,
                 Color{200, 200, 200, 255}); y += line;

        char buf[256];
        std::snprintf(buf, sizeof(buf),
                      "t = %.4f s   bodies = %zu   2.5PN = %s   %s",
                      scenario.system.time(),
                      scenario.system.bodies().size(),
                      scenario.system.uses_pn25() ? "ON" : "off",
                      paused ? "[PAUSED]" : "");
        DrawText(buf, x, y, 16, RAYWHITE); y += line;

        if (scenario.system.bodies().size() >= 2 &&
            scenario.name.find("Figure-8") == std::string::npos)
        {
            std::snprintf(buf, sizeof(buf),
                          "inner sep = %.0f km   f_gw = %.2f Hz   ecc = %.4f",
                          scenario.system.inner_separation() / 1000.0,
                          scenario.system.inner_gw_frequency(),
                          scenario.system.inner_eccentricity());
            DrawText(buf, x, y, 16, Color{200, 220, 255, 255});
        }

        // Footer with controls
        DrawRectangle(0, WINDOW_H - 28, WINDOW_W, 28, Color{0, 0, 0, 160});
        DrawText("[1] binary inspiral   [2] hierarchical triple   [3] figure-8   "
                 "[SPACE] pause   [R] reset   [drag] orbit camera   [scroll] zoom",
                 14, WINDOW_H - 22, 14, Color{180, 180, 180, 255});

        // Strain plot (only for inspiral)
        draw_strain_plot();
    }

    void draw_strain_plot() const {
        if (strain_history.size() < 2) return;
        const int W = 320;
        const int H = 90;
        const int x0 = WINDOW_W - W - 14;
        const int y0 = 14;
        DrawRectangle(x0, y0, W, H, Color{0, 0, 0, 160});
        DrawRectangleLines(x0, y0, W, H, Color{120, 120, 120, 255});
        DrawText("h+(t)  inner-binary strain",
                 x0 + 8, y0 + 6, 12, Color{200, 220, 255, 255});

        double peak = 0.0;
        for (double h : strain_history) peak = std::max(peak, std::fabs(h));
        if (peak == 0.0) return;

        const int plot_y0 = y0 + 22;
        const int plot_h  = H - 28;
        const int mid_y   = plot_y0 + plot_h / 2;
        DrawLine(x0 + 4, mid_y, x0 + W - 4, mid_y,
                 Color{80, 80, 80, 255});

        const int N = static_cast<int>(strain_history.size());
        for (int i = 1; i < N; ++i) {
            float fx0 = x0 + 4 +
                static_cast<float>(i - 1) / N * (W - 8);
            float fx1 = x0 + 4 +
                static_cast<float>(i)     / N * (W - 8);
            float fy0 = mid_y -
                static_cast<float>(strain_history[i - 1] / peak) * (plot_h / 2 - 2);
            float fy1 = mid_y -
                static_cast<float>(strain_history[i]     / peak) * (plot_h / 2 - 2);
            DrawLine(static_cast<int>(fx0), static_cast<int>(fy0),
                     static_cast<int>(fx1), static_cast<int>(fy1),
                     Color{120, 200, 255, 255});
        }
    }
};

} // namespace

int main() {
    InitWindow(WINDOW_W, WINDOW_H, "GW: Binary Inspiral / Triple System (3D)");
    SetTargetFPS(60);

    Sim sim;
    sim.load(make_binary_inspiral());

    OrbitCamera cam;

    while (!WindowShouldClose()) {
        // ---- Input ----
        cam.update();
        if (IsKeyPressed(KEY_ONE))   sim.load(make_binary_inspiral());
        if (IsKeyPressed(KEY_TWO))   sim.load(make_hierarchical_triple());
        if (IsKeyPressed(KEY_THREE)) sim.load(make_figure8());
        if (IsKeyPressed(KEY_R)) {
            // Re-load same scenario by copying the type
            std::string name = sim.scenario.name;
            if (name.find("Inspiral")    != std::string::npos) sim.load(make_binary_inspiral());
            else if (name.find("Triple") != std::string::npos) sim.load(make_hierarchical_triple());
            else                                                sim.load(make_figure8());
        }
        if (IsKeyPressed(KEY_SPACE)) sim.paused = !sim.paused;

        // ---- Step physics ----
        sim.step();

        // ---- Draw ----
        BeginDrawing();
        ClearBackground(Color{8, 10, 16, 255});

        Camera3D c = cam.camera();
        BeginMode3D(c);
            sim.draw_3d();
        EndMode3D();

        sim.draw_hud();
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
