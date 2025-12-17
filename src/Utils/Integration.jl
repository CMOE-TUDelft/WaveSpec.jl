module Integration

"""
Integration
===========
Module providing numerical integration methods.
"""

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
# ----------------------End----------------------

end