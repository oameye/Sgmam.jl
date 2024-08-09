using BenchmarkTools
using SplitApplyCombine
using LoopVectorization
using Strided

λ = 0.0168
ω0 = 1.000
ω = 1.000
γ = 1 / 295
η = 0
α = -1
function H_p(x)
  u, v = eachrow(x)

  H_p1 = @. -γ * u / 2 + ω * v / 2 + (-1 / 2 - λ / 4 - 0.375 * α * (u^2 + v^2)) * v / ω
  H_p2 = @. -γ * v / 2 - ω * u / 2 + (1 / 2 - λ / 4 + 0.375 * α * (u^2 + v^2)) * u / ω
  return Matrix([H_p1 H_p2]')
end

function KPO(x) # out-of-place
  u, v = x
  du = -γ * u / 2 - (3 * α * (u^2 + v^2) / 8 + (ω0^2 - ω^2) / 2 + λ / 4) * v / ω
  dv = -γ * v / 2 + (3 * α * (u^2 + v^2) / 8 + (ω0^2 - ω^2) / 2 - λ / 4) * u / ω
  return [du, dv]
end

# f(x) = KPO(x)
# g(x) = H_p(x)
function h(x)
  new = zeros(size(x))
  for i in 1:size(x, 2)
    new[:, i] = KPO(view(x, :, i))
  end
  return new
end
function i(x)
  x_vv = splitdimsview(x)
  return combinedimsview(KPO.(x_vv))
end

x = rand(2, 1000)
@benchmark reduce(hcat, KPO.(eachcol($x))) # 753.249 μs
@benchmark mapslices(KPO, $x, dims=1) # 1.06 ms
@benchmark h($x) #  815.349 μs
@benchmark i($x) # 760.491 μs
@benchmark H_p($x) # 11.9 μs
