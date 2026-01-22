# Project-ES313

# Frank–Kamenetskii Thermal Ignition Model (2D)

This repository contains a **numerical implementation of a 2-D reaction–diffusion model for thermal ignition**, inspired by **Frank–Kamenetskii (FK) theory**.
The project investigates **thermal runaway, critical temperature, and ignition time** using a cellular-automaton–style explicit time integration in Julia.

The model is intended for **theoretical and educational purposes**, focusing on **slow thermal ignition (cook-off)** rather than detonation or shock initiation.

---

## 🔬 Physical model

The temperature field (T(x,y,t)) is governed by the heat equation with a volumetric Arrhenius heat source:

[
\rho c_p \frac{\partial T}{\partial t}
======================================

\nabla \cdot (k \nabla T)
+
\rho Q A \exp!\left(-\frac{E}{RT}\right)
]

where:

* heat diffusion competes with
* exponentially temperature-dependent chemical heat release.

This equation forms the basis of **Frank–Kamenetskii thermal explosion theory**, which predicts the existence of a **critical condition separating decay from thermal runaway**.

---

## 📐 Frank–Kamenetskii framework

Frank–Kamenetskii theory shows that ignition is controlled by a single dimensionless parameter:

[
\delta
\sim
\frac{Q A L^2}{\rho c_p \alpha T^2}
\exp!\left(-\frac{E}{RT}\right)
]

For a 2-D infinite plane, ignition occurs when:

[
\delta > \delta_c \approx 1
]

In this project:

* FK theory is used to **interpret** results analytically
* the **critical temperature is determined numerically** via simulation

---

## 🧮 Numerical method

* **Spatial discretization:** 2-D Cartesian grid (finite differences, 5-point stencil)
* **Time integration:** explicit Euler (cellular-automaton–style update)
* **Boundary conditions:** zero-flux (infinite-domain approximation)
* **Initial condition:** localized hot spot (circular or line geometry)

The simulation captures:

* subcritical decay
* near-critical metastability
* supercritical thermal runaway

---

## 🔥 Ignition criterion

Ignition is **not** defined by a fixed temperature threshold.

Instead, ignition is detected as:

* **sustained accelerated growth of the maximum temperature**
* corresponding to the **absence of a steady solution**, consistent with FK theory

A numerical temperature cutoff is used only as a **safety stop condition**.

---

## 📊 Output and visualization

The code produces:

* real-time 2-D temperature field visualizations
* maximum temperature vs time plots
* optional animations (MP4 / GIF)

These outputs allow direct comparison with FK predictions such as:

* critical behavior
* divergence of ignition time near threshold

---

## ⚠️ Scope and limitations

This model:

* ✔ describes slow thermal ignition (cook-off)
* ✔ captures FK critical behavior
* ❌ does **not** model detonation, shock initiation, gas dynamics, or phase change
* ❌ uses a single-step Arrhenius approximation

All parameters are used for **numerical modeling only**.

---

## 🚀 How to run

Requirements:

* Julia ≥ 1.9
* Packages: `Plots`, `Statistics`

Run the main script:

```bash
julia inf_plane_chatgpt.jl
```

Adjust parameters such as grid size, time step, and initial hotspot temperature in the **INPUT** section of the script.

---

## 📚 References

* D. A. Frank-Kamenetskii, *Diffusion and Heat Transfer in Chemical Kinetics*, Plenum Press
* Law, C. K., *Combustion Physics*, Cambridge University Press
* Wikipedia: *Arrhenius equation*, *Heat equation*, *Thermal explosion theory*

---

## ✍️ Author

[Your Name]
[University / Course / Project name]

---

## 📄 License

This project is released for academic and educational use.
