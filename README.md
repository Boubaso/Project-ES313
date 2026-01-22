# Project-ES313

# Frank–Kamenetskii Thermal Ignition Model (2D)

This repository contains a **numerical implementation of a 2-D reaction–diffusion model for thermal ignition**, inspired by **Frank–Kamenetskii (FK) theory**.
The project investigates **thermal runaway, critical temperature, and ignition time** using a cellular-automaton–style explicit time integration in Julia.

The detail description, calculations, sources and background info are given in the Notion Website: 

https://prairie-aries-a8b.notion.site/ES313-Project-2025-b00002e455564949bd950bd71b22c402?source=copy_link 

---

## ⚠️ Scope and limitations

This model:

* ✔ describes slow thermal ignition (cook-off)
* ✔ captures FK critical behavior
* ❌ does **not** model detonation, shock initiation, gas dynamics, or phase change
* ❌ uses a single-step Arrhenius approximation

---

## 🚀 How to run

Requirements:

* Julia ≥ 1.9
* Packages: `Plots`, `Statistics`

Run the main script (main purpose):

```bash
ES313_PROJECT_ADJT_KBO_SOW_178POL\CA numerical MIOT ALL FIRE NO FIRE.jl
```

This gives the best calculations, the others are purely for single calculations and visualization/educational purposes.
Adjust parameters such as grid size, time step, and initial hotspot temperature in the **INPUT** section of the script.

---


## ✍️ With kind regards

[Sow Boubacar]
[Royal Military Academy / ES313]

---

