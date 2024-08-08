function sgmam_Hpp(sys::System, x_initial;
        ϵ::Float64 = 1e-1,
        iterations::Int64 = 1000,
        show_progress::Bool = false,
        save_info::Bool = false,
        reltol::Float64 = NaN,
        has_invH_pp::Bool = false)
    if has_invH_pp
        @unpack H_p, H_x, invH_pp = sys
    else
        @unpack H_p, H_x = sys
    end

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

    progress = Progress(iterations; enabled = show_progress)
    for i in 1:iterations
        if has_invH_pp
            update!(x, xdot, p, pdot, lambda, H_x, H_p, invH_pp, ϵ)
        else
            update!(x, xdot, p, pdot, lambda, H_x, H_p, ϵ)
        end

        # reparametrize to arclength
        interpolate_path!(x, alpha, s) # 14% of time
        push!(S, FW_action(xdot, p))

        tol = abs(S[end] - S[1]) / S[end]
        if tol < reltol
            @info "Converged after $i iterations with $tol"
            break
        end
        next!(progress)
    end
    return (x, S[end], lambda, p, xdot)
end

function interpolate_path!(path, α, iH, s) # 14% of time
    α[2:end] .= vec(sqrt.(sum(iH[:, 2:end] .* diff(path; dims = 2) .^ 2, dims = 1)))
    α .= cumsum(α; dims = 1)
    α .= α ./ last(α)
    path[1, :] .= LinearInterpolation(α, path[1, :])(s)
    path[2, :] .= LinearInterpolation(α, path[2, :])(s)
    nothing
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
    return nothing
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

        A = spdiagm( # spdiagm makes it significantly faster
            0 => 1 .+ 2 .* ϵ .* iH[dof, 2:(end - 1)] .* λ[2:(end - 1)] .^ 2,
            1 => -ϵ .* iH[dof, 2:(end - 2)] .* λ[2:(end - 2)] .^ 2,
            -1 => -ϵ .* iH[dof, 3:(end - 1)] .* λ[3:(end - 1)] .^ 2
        )
        prob = LinearProblem(A, rhs)
        x[dof, 2:(end - 1)] .= solve(prob, KLUFactorization()).u
    end
    return nothing
end

function update!(x, xdot, p, pdot, lambda, H_x, H_p, invH_pp, ϵ)
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
end
