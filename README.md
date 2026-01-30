# WaveSpec

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://shagunTUD.github.io/WaveSpec.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://shagunTUD.github.io/WaveSpec.jl/)
[![Build Status](https://github.com/shagunTUD/WaveSpec.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/shagunTUD/WaveSpec.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/shagunTUD/WaveSpec.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/shagunTUD/WaveSpec.jl)

---

## Description

WaveSpec provides functions related to ocean-wave input for numerical models.
These include

| Function |  Description |
| ---- | ---- |
| `dispersionRel(h, T)` | Function for calculation the wave-length for a ocean wave with period `T` and in water-depth `h`,  |
| `jonswap(Hs, Tp)` | Function for calculation jonswap spectrum, for a given significant wave height `Hs` and peak time period `Tp`) |
| `waveAiry1D(sp, t, x, z)`  | Generate time-series for wave-elevation and particle velocity at location `(x,z)` at at time-instances `t` for a ocean-wave spectrum defined by `sp` |

---

## Tutorials

Please refer to the `test` folder for the tutorials

---

## Installation 
WaveSpec is a registered package in the official [Julia package registry](https://github.com/JuliaRegistries/General). Thus, the installation of WaveSpec is straight forward using the Julia's package manager. Open the Julia REPL, type `]` to enter package mode, and install as follows
```
pkg> add WaveSpec
```


### Guidelines for Using `Revise.jl` in Development  

[`Revise.jl`](https://timholy.github.io/Revise.jl/stable/) enables developers to modify dependencies and immediately reflect those changes in the active Julia REPL without requiring a restart.  

While `Revise.jl` is a powerful tool, it is **not recommended** to add it as a package dependency. Instead, you should configure it to load automatically in your local Julia REPL. You can follow the instructions [here](https://timholy.github.io/Revise.jl/stable/config/#Using-Revise-by-default) to set this up.  

This approach ensures that `Revise.jl` is available in your local environment, allowing you to take advantage of its features without adding it to the package itself.


---

## Contact

Please contact [Oriol Colomes](https://www.tudelft.nl/en/staff/j.o.colomesgene/?cHash=d85db1dfc98f5e255324852d31948ede), [Shagun Agarwal](https://shagun751.github.io/) for further questions.

---

