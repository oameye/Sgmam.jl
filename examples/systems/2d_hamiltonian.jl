function H(x, p) # ℜ² → ℜ
    u, v = eachrow(x)
    pu, pv = eachrow(p)
    return @. (pu^2 + pv^2)/2 + pv*fu(u, v) + pu*fv(u, v)
end
function H_x(x, p) # ℜ² → ℜ²
    u, v = eachrow(x)
    pu, pv = eachrow(p)

    H_u = @. pu*dfudu(u, v) + pv*dfvdu(u, v)
    H_v = @. pu*dfudv(u, v) + pv*dfvdv(u, v)
    return Matrix{eltype(x)}([H_u H_v]')
end
function H_p(x, p) # ℜ² → ℜ²
    u, v = eachrow(x)
    pu, pv = eachrow(p)

    H_pu = @. pu + fu(u, v)
    H_pv = @. pv + fv(u, v)
    return Matrix{eltype(x)}([H_pu H_pv]')
end
