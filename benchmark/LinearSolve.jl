# 2x in speed when using LinearSolve

using Sgmam, Plots
using NLsolve: nlsolve
using BenchmarkTools

const λ = 3 / 1.21 * 2 / 295
const ω0 = 1.000
const ω = 1.000
const γ = 1 / 295
const η = 0
const α = -1

function get_steady_state(sys, x0)
    steady_state(x) = sys.H_p(x, zeros(size(x)))
    return Xs = nlsolve(steady_state, x0).zero
end

include("systems/KPO.jl")
sys = System(H, H_x, H_p)

# setup
Nt = 500  # number of discrete time steps
s = collect(range(0; stop=1, length=Nt))

xa = get_steady_state(sys, [1.0, 1.0])
xb = -xa
xsaddle = [0.0, 0.0]

# Initial trajectory
xx = @. (xb[1] - xa[1]) * s + xa[1] + 4 * s * (1 - s) * xsaddle[1]
yy = @. (xb[2] - xa[2]) * s + xa[2] + 4 * s * (1 - s) * xsaddle[2] + 0.01 * sin(2π * s)
x_initial = Matrix([xx yy]')

x_min, S_min, lambda, p, xdot = sgmam(
    sys, x_initial; iterations=1_000, ϵ=10e2, show_progress=false
)
@show S_min;
# plot(x_initial[1, :], x_initial[2, :], label = "init", lw=3, c=:black)
# plot!(x_min[1, :], x_min[2, :], label = "MLP", lw=3, c=:red)

@benchmark $sgmam($sys, $x_initial, iterations=100, ϵ=10e2, show_progress=false)
@profview sgmam(sys, x_initial, iterations=100, ϵ=10e2, show_progress=false)

# @benchmark sgmam($params, $x_initial, iterations=1e3, ϵ=10e2)
# @profile sgmam(params, x_initial, iterations=1e3, ϵ=10e2, reltol=1e-10)
# owntime()
# owntime(stackframe_filter=filecontains(srcdir("sgmam.jl")))
# totaltime()

# @profview sgmam_opt(params, x_initial, iterations=1e3, ϵ=10e2, reltol=1e-10)
ϵ = 10e2
iH = sys.invH_pp(x_min, ())
dof = 1
x = x_min
idxc = 2:(Nt - 1)
xdotdot = zeros(size(x))
Sgmam.central_diff!(xdotdot, xdot)
pdot = zeros(size(x))
Sgmam.central_diff!(pdot, p)
Hx = sys.H_x(x, p)

rhs = @. (
    x_min[dof, idxc] +
    ϵ * (
        lambda[idxc] * pdot[dof, idxc] + Hx[dof, idxc] -
        iH[dof, idxc] * lambda[idxc]^2 * xdotdot[dof, idxc]
    )
)
rhs[1] += ϵ * iH[dof, 2] * lambda[2]^2 * xa[dof]
rhs[end] += ϵ * iH[dof, end - 1] * lambda[end - 1]^2 * xb[dof]

A = spdiagm( # spdiagm makes it significantly faster
    0 => 1 .+ 2 .* ϵ .* iH[1, 2:(end - 1)] .* lambda[2:(end - 1)] .^ 2,
    1 => -ϵ .* iH[1, 2:(end - 2)] .* lambda[2:(end - 2)] .^ 2,
    -1 => -ϵ .* iH[1, 3:(end - 1)] .* lambda[3:(end - 1)] .^ 2,
)
using LinearSolve
cond(Array(A), 2)

prob = LinearProblem(A, rhs)
@btime $solve($prob);
@btime $solve($prob, assumptions=OperatorAssumptions(true));
@btime $solve($prob, KrylovJL_GMRES());
@btime $solve($prob, UMFPACKFactorization());
@btime $solve($prob, KLUFactorization());
@btime $A \ $rhs;
