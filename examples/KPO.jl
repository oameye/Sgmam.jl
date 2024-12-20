using Sgmam, Plots, NLsolve
using LinearSolve
using BenchmarkTools

function get_steady_state(sys, x0)
    steady_state(x) = sys.H_p(x, zeros(size(x)))
    return Xs = nlsolve(steady_state, x0).zero
end

include("systems/KPO.jl")
sys = System(H_x, H_p)

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

@btime $sgmam($sys, $x_initial, iterations=100, ϵ=10e2, show_progress=false) # 25.803 ms (29024 allocations: 105.69 MiB)
@profview sgmam(sys, x_initial, iterations=100, ϵ=10e2, show_progress=false)

# The bottleneck is atm at the LinearSolve call to update the x in the new iteration. So the more improve, one needs to write it own LU factorization.
