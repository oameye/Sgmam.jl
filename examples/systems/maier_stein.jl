#! format: off

const β = 10
const α = 1.0

fu(u, v) = u - u^3 - β*u*v^2
fv(u, v) = -α*(1 + u^2)*v
dfudv(u, v) = - 2*β*u*v
dfudu(u, v) = 1 - 3*u^2 - β*v^2
dfvdv(u, v) = -α*(1+u^2)
dfvdu(u, v) = -2*α*u*v

include("2d_hamiltonian.jl")
#! format: on
