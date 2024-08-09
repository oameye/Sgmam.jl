using Sgmam, Plots, NLsolve
using BenchmarkTools

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
plot(x_initial[1, :], x_initial[2, :], label = "init", lw=3, c=:black)
plot!(x_min[1, :], x_min[2, :], label = "MLP", lw=3, c=:red)

@benchmark $sgmam($sys, $x_initial, iterations=100, ϵ=10e2, show_progress=false)
# @profview sgmam(sys, x_initial, iterations=100, ϵ=10e2, show_progress=false)
# v_range = range(-.15, .15, 100); u_range = range(-.05, .05, 100)
# iter = Iterators.product(u_range, v_range) |> collect
# heatmap(u_range,v_range,map( x-> H([x...], ones(2))[1], iter)')
# plot!(x_initial[1, :], x_initial[2, :], label = "init", lw=3, c=:black)
# plot!(x_min[1, :], x_min[2, :], label = "MLP", lw=3, c=:red)
