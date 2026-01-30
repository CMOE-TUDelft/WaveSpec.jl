module PhysicalConstants

# We use 'const' for performance. 
# Changing these during a session requires a restart in Julia.

export g, ρ_water, ρ_air, ν_water, σ_water

# --- Physical Constants ---

"""
    g
Standard gravity (m/s²)
"""
const g = 9.80665

"""
    ρ_water
Density of seawater (kg/m³). 
Standard value at 15°C and 3.5% salinity.
"""
const ρ_water = 1025.0

"""
    ρ_air
Density of air (kg/m³) at sea level.
"""
const ρ_air = 1.225

"""
    ν_water
Kinematic viscosity of water (m²/s) at 15°C.
"""
const ν_water = 1.19e-6

"""
    σ_water
Surface tension of water-air interface (N/m).
"""
const σ_water = 0.0728

# --- Mathematical Constants (For convenience in 3D logic) ---

export RAD_TO_DEG, DEG_TO_RAD

const RAD_TO_DEG = 180.0 / π
const DEG_TO_RAD = π / 180.0

end # module