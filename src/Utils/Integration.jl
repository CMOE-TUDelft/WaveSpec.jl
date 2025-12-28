module Integration

"""
Integration
===========
Module providing numerical integration methods.
"""

export simpsonInteg1D, trapzInteg1D, gaussQuad1D, gaussQuad1D_4

# ---------------------Start---------------------
function simpsonInteg1D(y, dx)
  ind = 1:length(y)  
  w = ifelse.(iseven.(ind), 4, 2)  
  w[1] = 1
  w[end] = 1

  return sum(dx/3.0 * w .* y)
end


function trapzInteg1D(y, dx)
  ind = 1:length(y)  
  w = ifelse.(iseven.(ind), 2, 2)  
  w[1] = 1
  w[end] = 1

  return sum(dx/2.0 * w .* y)
end


function gaussQuad1D(y, dx)
  # Gauss-Legendre Quadrature 2 point
  q1dx = 0.5*(1 - 1/√(3))
  q2dx = 0.5*(1 + 1/√(3))

  ly = y[1:end-1]
  ry = y[2:end]

  yq1 = (ry - ly)*q1dx + ly
  yq2 = (ry - ly)*q2dx + ly

  return sum(yq1 + yq2)*dx/2
end

"""
    gauss_quad_4(f, a, b)

Integrates function `f` from `a` to `b` using a 4-point Gauss-Legendre rule.
The mapping used is: 
    x = (b-a)/2 * ξ + (a+b)/2
    dx = (b-a)/2 * dξ
"""
function gaussQuad1D_4(f::Function, a::Float64, b::Float64, N::Int64)
    nodes = (-0.8611363115940526, -0.3399810435848563,  0.3399810435848563,  0.8611363115940526)
    weights = (0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538)

    h_diff = (b - a) / 2.0
    h_sum  = (b + a) / 2.0
    
    bin_width = (b - a) / N
    total_integral = 0.0

    for i in 1:N
        # Define bounds for the current bin
        bin_a = a + (i - 1) * bin_width
        bin_b = a + i * bin_width
        
        # Mapping constants for the current bin
        h_diff = bin_width / 2.0
        h_sum  = (bin_a + bin_b) / 2.0
        
        # Apply 4-point Gauss rule to this bin
        bin_integral = 0.0
        for j in 1:4
            x = h_diff * nodes[j] + h_sum
            bin_integral += weights[j] * f(x)
        end
        total_integral += bin_integral * h_diff
    end
    
    return total_integral
end
# ----------------------End----------------------

end