module Sgmam

export sgmam

using DataStructures: CircularBuffer
using ProgressMeter: Progress, next!
using Dierckx: ParametricSpline
# using ProgressBars: tqdm

using LinearAlgebra, SparseArrays

using DispatchDoctor: @stable

@stable begin # enforces type_stability


function central_diff!(xdot, x)
    # ̇xₙ = 0.5(xₙ₊₁ - xₙ₋₁) central finite difference
    xdot[:, 2:(end - 1)] = 0.5 * (x[:, 3:end] - x[:, 1:(end - 2)])
    nothing
end

FW_action(xdot, p) = sum(sum(xdot .* p; dims = 1)) / 2

function sgmam(sys, x_initial;
        xsaddle = [],
        ϵ = 1e-1,
        iterations = 1000,
        show_progress = false,
        save_info = false,
        reltol = NaN)

    H_p = sys["H_p"] # H_p(x,p) = b(x) + a*p
    H_x = sys["H_x"] # H_x(x, p) = ForwardDiff.jacobian(b, x)'*p
    invH_pp = sys["invH_pp"] # invH_pp(x, p) = a

    Nx, Nt = size(x_initial)
    s = range(0, stop = 1, length = Nt)

    # preallocate
    x = deepcopy(x_initial)
    p = zeros(size(x))
    pdot = zeros(size(x))
    xdot = zeros(size(x))
    lambda = zeros(1, Nt) # Lagrange multiplier
    alpha = zeros(Nt)

    S = CircularBuffer{Float64}(2)
    fill!(S, Inf)

    progress = Progress(iterations; enabled=show_progress)
    for i in 1:iterations
        central_diff!(xdot, x)

        update_p!(p, lambda, x, xdot, H_p, invH_pp)

        central_diff!(pdot, p)
        Hx = H_x(x, p)

        # explicit update (deactivated, using implicit instead)
        # x = x + ϵ*(lambda .* pdot + Hx);

        # implicit update
        iH = invH_pp(x, p)
        xdotdot = zeros(size(xdot))
        central_diff!(xdotdot, xdot)

        # each dof has same lambda, but possibly different H_pp^{-1}
        update_x!(x, lambda, pdot, xdotdot, Hx, iH, ϵ) # 63% of time

        # reparametrize to arclength
        interpolate_path!(x, alpha, iH, s) # 14% of time
        push!(S, FW_action(xdot, p))

        tol = abs(S[end] - S[1]) / S[end]
        if tol < reltol
            @info "Converged after $i iterations with $tol"
            break
        end
        next!(progress)
    end
    return save_info ? (x, S[end], lambda, p, xdot) : (x, S[end])
end

function interpolate_path!(x, alpha, iH, s) # 14% of time
    alpha[2:end] .= vec(sqrt.(sum(iH[:, 2:end] .* diff(x, dims = 2) .^ 2, dims = 1)))
    alpha .= cumsum(alpha, dims = 1)
    alpha .= alpha ./ last(alpha)
    interp = ParametricSpline(vec(alpha), x)
    x .= Matrix(interp(s))
    nothing
end

function update_x!(x, λ, p′, x′′, Hx, iH, ϵ) # 63% of time
    # each dof has same lambda, but possibly different H_pp^{-1}
    Nx, Nt = size(x)
    xa = x[:, 1]
    xb = x[:, end]

    # rhs = zeros(Nt)
    idxc = 2:(Nt - 1)
    for dof in 1:Nx
        rhs = @. (x[dof, idxc] +
                  ϵ * (λ[idxc] * p′[dof, idxc] + Hx[dof, idxc] -
                   iH[dof, idxc] * λ[idxc]^2 * x′′[dof, idxc]))
        rhs[1] += ϵ * iH[dof, 2] * λ[2]^2 * xa[dof]
        rhs[end] += ϵ * iH[dof, end - 1] * λ[end - 1]^2 * xb[dof]

        A = spdiagm(
            0 => 1 .+ 2 .* ϵ .* iH[dof, 2:(end - 1)] .* λ[2:(end - 1)] .^ 2,
            1 => -ϵ .* iH[dof, 2:(end - 2)] .* λ[2:(end - 2)] .^ 2,
            -1 => -ϵ .* iH[dof, 3:(end - 1)] .* λ[3:(end - 1)] .^ 2
        )
        x[dof, 2:(end - 1)] .= A \ rhs
    end
end

function update_p!(p, lambda, x, xdot, H_p, invH_pp)
    # Alternative: Direct computation, only correct for quadratic Hamiltonian in p,
    # where H_pp does not depend on p
    b_ = H_p(x, 0 * x)
    invHpp_ = invH_pp(x, 0 * x)
    lambda .= sqrt.(sum(b_ .^ 2 .* invHpp_, dims = 1) ./
                    sum(xdot .^ 2 .* invHpp_, dims = 1))
    lambda[1] = 0
    lambda[end] = 0
    p .= invHpp_ .* (lambda .* xdot .- b_)
    nothing
end

end # @stable
end
