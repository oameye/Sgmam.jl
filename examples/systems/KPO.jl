#! format: off

const λ = 3 / 1.21 * 2 / 295
const ω0 = 1.000
const ω = 1.000
const γ = 1 / 295
const η = 0
const α = -1

fu(u, v) = -γ*v/2 - λ*u/(4*ω) + (1 - ω^2)*u/(2*ω) + 3*α*u*(u^2 + v^2)/(8*ω) - η*v*(u^2 + v^2)/8
fv(u, v) = -γ*u/2 - λ*v/(4*ω) - (1 - ω^2)*v/(2*ω) - 3*α*v*(u^2 + v^2)/(8*ω) - η*u*(u^2 + v^2)/8
dfudv(u, v) = -(4*γ*ω - 6*α*u*v + η*ω*(u^2 + 3*v^2))/(8*ω)
dfudu(u, v) = (4 - 2*λ - 4*ω^2 + 9*α*u^2 - 2*η*ω*u*v + 3*α*v^2)/(8*ω)
dfvdv(u, v) = -(4 + 2*λ - 4*ω^2 + 3*α*u^2 + 2*η*ω*u*v + 9*α*v^2)/(8*ω)
dfvdu(u, v) = -(4*γ*ω + 6*α*u*v + η*ω*(3*u^2 + v^2))/(8*ω)


function H(x, p)
  u, v = eachrow(x)
  pv, pu = eachrow(p)

  return @. (pv^2 + pu^2)/2 + pu*fu(u, v) + pv*fv(u, v)
end
function H_x(x, p)
  u, v = eachrow(x)
  pv, pu = eachrow(p)

  H_u = @. pu*dfudu(u, v) + pv*dfvdu(u, v)
  H_v = @. pv*dfvdv(u, v) + pu*dfudv(u, v)
  return Matrix([H_u H_v]')
end

function H_p(x, p)
  u, v = eachrow(x)
  pv, pu = eachrow(p)

  H_pv = @. pv + fv(u, v)
  H_pu = @. pu + fu(u, v)
  return Matrix([H_pv H_pu]')
end
#! format: on
