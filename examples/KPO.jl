using Sgmam, Plots
using NLsolve: nlsolve
using BenchmarkTools

const λ = 3/1.21*2/295
const ω0 = 1.000
const ω = 1.000
const γ = 1/295
const η = 0
const α = -1

function get_steady_state(sys, x0)
    steady_state(x) = sys.H_p(x, zeros(size(x)))
    Xs = nlsolve(steady_state, x0).zero
end

include("systems/KPO.jl")
sys = System(H, H_x, H_p)

# setup
Nt = 2^10  # number of discrete time steps
s = range(0, stop=1, length=Nt) |> collect

xa = get_steady_state(sys, [1.0, 1.0])
xb = -xa
xsaddle = [0.0,0.0]

# Initial trajectory
xx = @. (xb[1] - xa[1])*s + xa[1] + 4*s*(1 - s)*xsaddle[1]
yy = @. (xb[2] - xa[2])*s + xa[2] + 4*s*(1 - s)*xsaddle[2] + 0.01*sin(2π*s)
x_initial = Matrix([xx yy]')

x_min, S_min = sgmam(sys, x_initial, iterations=1000, ϵ=10e2, show_progress=false)
@show S_min;
plot(x_initial[1, :], x_initial[2, :], label = "init", lw=3, c=:black)
plot!(x_min[1, :], x_min[2, :], label = "MLP", lw=3, c=:red)

@benchmark sgmam(sys, x_initial, iterations=100, ϵ=10e2, show_progress=false)
@profview sgmam(sys, x_initial, iterations=100, ϵ=10e2, show_progress=false)

# @benchmark sgmam($params, $x_initial, iterations=1e3, ϵ=10e2)
# @profile sgmam(params, x_initial, iterations=1e3, ϵ=10e2, reltol=1e-10)
# owntime()
# owntime(stackframe_filter=filecontains(srcdir("sgmam.jl")))
# totaltime()

# @profview sgmam_opt(params, x_initial, iterations=1e3, ϵ=10e2, reltol=1e-10)
