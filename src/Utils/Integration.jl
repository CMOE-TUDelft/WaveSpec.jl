module Integration

"""
Integration
===========
Module providing numerical integration methods.
"""

export IntegrateTrapezoidal, IntegrateGaussQuad

"""
    IntegrateTrapezoidal(f; a, b, n)

Integrates function `f` from `a` to `b` using the Composite Trapezoidal Rule 
with `n` sub-intervals (panels).
"""
function IntegrateTrapezoidal(f::Function; a::Real=-1.0, b::Real=1.0, n::Int=100)
    # Check integration bounds
    if !isfinite(a) || !isfinite(b)
        throw(ArgumentError("Trapezoidal rule requires finite bounds. Received a=$a, b=$b"))
    end
    if n < 1
        throw(ArgumentError("Number of intervals n must be at least 1."))
    end

    bin_width = (b - a) / n
    total_integral = 0.0

    # Integration Loop
    for i in 1:n
        x_left  = a + (i - 1) * bin_width
        x_right = a + i * bin_width
        
        # Area of trapezoid: width * average height
        total_integral += (f(x_left) + f(x_right))
    end
    
    # We multiply by width/2 at the end for efficiency 
    return total_integral * (bin_width / 2.0)
end


"""
    IntegrateGaussQuad(f, order, a, b, n)

Integrates function `f` from `a` to `b` using a composite Gauss-Legendre rule of order 'order' and 
and `n` sub-intervals (panels).
The mapping used is: 
    x = (b-a)/2 * ξ + (a+b)/2
    dx = (b-a)/2 * dξ
"""
function IntegrateGaussQuad(f::Function; order::Int=2, a::Real=-1.0, b::Real=1.0, n::Int=100)
    # Check integration bounds
    if !isfinite(a) || !isfinite(b)
        throw(ArgumentError("Gauss-Legendre rule requires finite bounds. Received a=$a, b=$b"))
    end
    if n < 1
        throw(ArgumentError("Number of intervals n must be at least 1."))
    end

    # Order of Gauss Quadrature integration 
    if order == 1
      nodes = ( 0.0, )
      weights = ( 2.0, )
    elseif order == 2
      nodes = ( - 1/√(3), 1/√(3) )
      weights = ( 1.0, 1.0)
    elseif order == 3
      nodes = ( - √(3/5), 0.0,  √(3/5))
      weights = ( 5/9, 8/9, 5/9 )
    elseif order == 4
      nodes = (-0.8611363115940526, -0.3399810435848563,  0.3399810435848563,  0.8611363115940526)
      weights = (0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538)
    else 
      error("Gauss Quadrature for that order not implemented.")
    end

    bin_width = (b - a) / n
    total_integral = 0.0

    for i in 1:n
        # Define bounds for the current bin
        bin_a = a + (i - 1) * bin_width
        bin_b = a + i * bin_width
        
        # Mapping constants for the current bin
        h_diff = bin_width / 2.0
        h_sum  = (bin_a + bin_b) / 2.0
        
        # Apply 4-point Gauss rule to this bin
        bin_integral = 0.0
        for j in 1:order
            x = h_diff * nodes[j] + h_sum
            bin_integral += weights[j] * f(x)
        end
        total_integral += bin_integral * h_diff
    end
    
    return total_integral
end
# ----------------------End----------------------

end