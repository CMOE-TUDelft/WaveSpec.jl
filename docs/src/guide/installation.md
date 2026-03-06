# Installation

WaveSpec is a registered package in the official [Julia package registry](https://github.com/JuliaRegistries/General). Thus, the installation of WaveSpec is straightforward using Julia's package manager.

## Installing via Pkg REPL

Open the Julia REPL, type `]` to enter package mode, and install as follows:

```julia
pkg> add WaveSpec
```

## Installing via Pkg API

Alternatively, you can use the `Pkg` API in your scripts:

```julia
using Pkg
Pkg.add("WaveSpec")
```

## Development Version

If you want to use the latest development version from the master branch, you can install it directly from the GitHub repository:

```julia
pkg> add https://github.com/CMOE-TUDelft/WaveSpec.jl
```
