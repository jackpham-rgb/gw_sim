# gw_sim — Binary Inspirals & Three-Body Gravitation in C++

A from-scratch C++17 simulator for compact-binary inspirals, hierarchical
triples (Kozai–Lidov), and the Chenciner–Montgomery figure-8 three-body
orbit. Three different ways to look at the result:

| Target            | What it does                                        |
|-------------------|-----------------------------------------------------|
| `make run`        | Computes a single binary inspiral, writes a CSV     |
|                   | of strain h+(t), h×(t), f_gw(t), separation(t),     |
|                   | and a `.wav` of the chirp.                          |
| `make visualize`  | Real-time 2D ASCII animation of a single binary,    |
|                   | with an ASCII strain plot at the end.               |
| `make visualize3d`| Real-time 3D OpenGL window (raylib). Switch         |
|                   | between binary, hierarchical triple, and figure-8   |
|                   | with the number keys; orbit camera with the mouse.  |

The physics engine is shared and self-contained — no scientific
libraries used. All formulas are documented inline in the source.

---

## Two physics engines

The project contains **two** different physics modules, used by
different targets:

### `src/binary_system.{hpp,cpp}` — orbit-averaged 2-body

Used by `gw_sim` and `gw_visualize`. Tracks just the orbital separation
`r` and phase `phi` of a circular binary; integrates Peters' (1964)
orbit-averaged equation `dr/dt ∝ -1/r^3` with RK4. Very fast, very
accurate for the inspiral phase, but assumes the orbit stays exactly
circular.

### `src/nbody.{hpp,cpp}` — full 3D N-body with 2.5PN

Used by `gw_3d`. Integrates the actual 3D positions and velocities of
every body, with full Newtonian pairwise gravity plus optional 2.5
post-Newtonian radiation reaction on a designated "inner binary." This
is what real numerical-relativity codes use for the early inspiral, and
it lets us simulate triple systems where Newtonian three-body dynamics
matter on the same plot as the 2-body GW radiation.

---

## Physics in detail

### A. The orbit-averaged binary (engine 1)

Two point masses on a circular orbit shrink under quadrupole radiation
reaction (Peters 1964):

```
dr/dt  = -(64/5) G^3 m1 m2 (m1+m2) / (c^5 r^3)
dphi/dt = sqrt(G (m1+m2) / r^3)                    [Kepler's 3rd law]
```

Strain at the detector at distance D and inclination ι (Maggiore 2008,
eq. 4.30):

```
h+ = (4/D) (G M_c / c^2)^(5/3) (π f_gw / c)^(2/3)
       * (1 + cos²ι) / 2 * cos(2 phi)

hx = (4/D) (G M_c / c^2)^(5/3) (π f_gw / c)^(2/3)
       * cos(ι)            * sin(2 phi)

M_c = (m1 m2)^(3/5) / (m1 + m2)^(1/5)              [chirp mass]
```

Cutoff at the Schwarzschild ISCO `r = 6 G M / c²`.

### B. The 3D N-body engine (engine 2)

For each body i,

```
a_i = sum_{j≠i} G m_j (r_j - r_i) / |r_j - r_i|^3   [Newtonian]
       + a_i^(2.5PN)                                  [if enabled]
```

The 2.5PN piece, applied only between the two designated "inner
binary" bodies, is the relative-acceleration form (Blanchet 2014,
Living Rev. Relativity 17, eq. 203):

```
a_rel = (8/5) G² M² ν / (c⁵ r³)
        × { [3 v² + 17 G M / (3 r)] ṙ n̂  −  [v² + 3 G M / r] v }
```

with `M = m1+m2`, `ν = m1 m2 / M²`, `n̂ = r/|r|`, `ṙ = v · n̂`. For a
circular orbit (`ṙ = 0`, `v² = GM/r`) this collapses to

```
a_rel = -(32/5) G³ M³ ν / (c⁵ r⁴) v
```

which reproduces Peters' orbit-averaged energy loss exactly. The split
between the two bodies follows Newton's third law:
`a_self = (m_other / M) a_rel`, `a_other = -(m_self / M) a_rel`.

The integrator is classical RK4 on the 6N-dimensional phase space.
This is good enough for visualization but is not symplectic — for
production work, swap in leapfrog / Wisdom–Holman.

### C. Three-body: Kozai–Lidov

In a hierarchical triple (close inner binary + distant third body on an
inclined orbit), the third body torques the inner binary. When the
mutual inclination exceeds ~39.2°, the inner-binary eccentricity and
inclination undergo periodic exchange — the **Kozai–Lidov mechanism**
(Kozai 1962, Lidov 1962). This is the leading astrophysical channel for
forming the close binaries that LIGO detects: the third body pumps the
inner orbit to high eccentricity, so close pericenter passages bleed
energy to gravitational waves and the binary inspirals from a
configuration that would otherwise take longer than the age of the
universe.

The hierarchical-triple scenario in this code uses a much tighter outer
orbit than astrophysical reality (so the cycle is visible in seconds
rather than megayears) and runs with PN off, so what you see is the
**Newtonian** precession and eccentricity exchange that drives the
phenomenon.

### D. Three-body: figure-8

A genuine periodic solution to the Newtonian three-body problem.
Discovered numerically by Cris Moore (1993) and proven to exist by
Chenciner & Montgomery (Annals of Mathematics, 2000): three equal masses
chase each other around a single figure-eight curve, exchanging roles
every third of a period. It has zero total angular momentum.

The standard initial conditions in `G = m = 1` units are:

```
r1 = ( 0.97000436, -0.24308753, 0)
r2 = -r1
r3 = (0, 0, 0)
v3 = (-0.93240737, -0.86473146, 0)
v1 = v2 = -v3 / 2
```

Verified test in this code: the figure-8 conserves total energy to
~10⁻¹⁴ over three full periods.

---

## Build & run

### Prerequisites

- A C++17 compiler (`g++ ≥ 7` or `clang ≥ 5`)
- `make`
- For the 3D target only: **raylib**
  ```bash
  sudo apt install libraylib-dev    # Ubuntu / Debian / WSL
  brew install raylib                # macOS
  ```

### 2D targets (no extra dependencies)

```bash
make            # builds gw_sim and gw_visualize
make run        # runs gw_sim, writes output/{waveform.csv, chirp.wav}
make visualize  # runs gw_visualize, ASCII orbit + strain chart
```

### 3D target

```bash
make 3d            # builds gw_3d (needs libraylib-dev installed)
make visualize3d   # runs it
```

### Cleanup

```bash
make clean
```

---

## 3D viewer controls

| Key / mouse        | Action                                       |
|--------------------|----------------------------------------------|
| **1**              | Load Binary Inspiral (GW150914-like)         |
| **2**              | Load Hierarchical Triple (Kozai–Lidov)       |
| **3**              | Load Figure-8 Three-Body                     |
| **R**              | Reset current scenario                       |
| **SPACE**          | Pause / unpause                              |
| Drag (left mouse)  | Orbit the camera                             |
| Mouse wheel        | Zoom                                         |
| ESC                | Quit                                         |

The HUD shows simulation time, separation, GW frequency, and orbital
eccentricity for the inner binary; for the inspiral scenario, a live
plot of h+(t) appears in the top right.

---

## Project layout

```
gw_sim/
├── Makefile
├── README.md
├── src/
│   ├── vec3.hpp                # 3D vector type (header-only)
│   ├── rk4.hpp                 # generic RK4 template (engine 1)
│   ├── binary_system.{hpp,cpp} # orbit-averaged 2-body (engine 1)
│   ├── wav_writer.hpp          # tiny WAV-file writer
│   ├── nbody.{hpp,cpp}         # 3D N-body + 2.5PN (engine 2)
│   ├── scenarios.hpp           # binary / triple / figure-8 presets
│   ├── main.cpp                # gw_sim:        CSV + WAV
│   ├── main_visualize.cpp      # gw_visualize:  ASCII animation
│   ├── main_3d.cpp             # gw_3d:         raylib 3D viewer
│   └── test_nbody.cpp          # physics sanity checks
└── output/                     # created at run time
```

Build the test binary explicitly:

```bash
g++ -std=c++17 -O2 -Isrc src/nbody.cpp src/test_nbody.cpp -o test_nbody
./test_nbody
```

Expected output: Newtonian binary energy drift ~10⁻¹¹, 2.5PN energy-loss
rate within 0.5% of Peters' formula, figure-8 drift ~10⁻¹⁴.

---

## Trying other events

In `src/scenarios.hpp`, edit `make_binary_inspiral()` to plug in
different masses / starting frequency. Some real events:

| Event       | m1 (M☉) | m2 (M☉) | D (Mpc) | Notes                                   |
|-------------|---------|---------|---------|-----------------------------------------|
| GW150914    | 36      | 29      | 410     | first-ever GW detection (2015)          |
| GW170817    | 1.46    | 1.27    | 40      | binary neutron star (different physics) |
| GW190521    | 85      | 66      | ~5300   | most massive merger seen at the time    |
| GW250114    | ~33     | ~32     | ~330    | the cleanest signal to date — fit       |
|             |         |         |         | values from arXiv:2601.13173            |

For GW250114 specifically, see the LIGO–Virgo–KAGRA 2026 paper.

---

## What's *not* in here

This is a serious starting point, not a research code. Things you might
want to add:

1. **Higher-order PN.** Currently 2.5PN. Real templates go to 3.5PN or
   higher. See Blanchet (2014).
2. **Spin.** Each black hole has a spin; precession and final-spin
   effects matter for real signals.
3. **Merger + ringdown.** This code stops at ISCO. Stitching on a
   phenomenological merger (IMRPhenomD) and a Kerr quasi-normal-mode
   ringdown is the natural next step.
4. **Symplectic integrator.** Replace RK4 with leapfrog or Wisdom–
   Holman so long N-body runs don't drift.
5. **Detector response.** LIGO arm geometry, antenna patterns, noise
   from the public PSDs.
6. **Comparison against real data.** Pull GW150914 / GW250114 strain
   from the GW Open Science Center (gwosc.org) and overlay your
   template.

---

## References

- **Einstein, A.** (1916, 1918) — general relativity and the prediction
  of gravitational waves.
- **Peters, P. C. & Mathews, J.** (1963) and **Peters, P. C.** (1964) —
  the quadrupole-formula back-reaction. *Phys. Rev.* **131**, 435 and
  **136**, B1224.
- **Damour, T. & Deruelle, N.** (1981) — radiation-reaction force at
  2.5PN. *Phys. Lett. A* **87**, 81.
- **Iyer, B. R. & Will, C. M.** (1995) — radiation-reaction in
  arbitrary gauges. *Phys. Rev. D* **52**, 6882.
- **Blanchet, L.** (2014) — definitive review of PN methods.
  *Living Rev. Relativity* **17**, 2.
- **Maggiore, M.** (2008) — *Gravitational Waves, Vol. 1: Theory and
  Experiments*, Oxford. Used for the strain conventions in §A.
- **Kozai, Y.** (1962); **Lidov, M. L.** (1962) — original
  Kozai–Lidov papers. *AJ* **67**, 591 and *Planet. Space Sci.* **9**, 719.
- **Naoz, S.** (2016) — modern review of the eccentric Kozai–Lidov
  mechanism. *Annu. Rev. Astron. Astrophys.* **54**, 441.
- **Moore, C.** (1993) — numerical discovery of the figure-8 orbit.
  *Phys. Rev. Lett.* **70**, 3675.
- **Chenciner, A. & Montgomery, R.** (2000) — existence proof.
  *Annals of Mathematics* **152**, 881.
- **LIGO–Virgo–KAGRA Collaboration** — GW150914 (PRL 116, 061102, 2016)
  and GW250114 (PRL 2026, arXiv:2601.13173).
- **GW Open Science Center**, https://gwosc.org — public data archive.
- **raylib**, https://www.raylib.com/ — used for 3D rendering.

---

## License

Public domain / CC0. Use it however you like.
