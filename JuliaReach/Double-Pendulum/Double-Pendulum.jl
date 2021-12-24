module DoublePendulum

using NeuralNetworkAnalysis, LaTeXStrings
using NeuralNetworkAnalysis: Specification
import Plots, DifferentialEquations
using Plots: plot, plot!, savefig, font, Measures.mm

# problem

m = 0.5
L = 0.5
c = 0.0
g = 1.0
gL = g/L
mL = m*L^2

@taylorize function double_pendulum!(dx, x, p, t)
    x₁, x₂, x₃, x₄, T₁, T₂ = x

    Δ12 = x₁ - x₂
    ★ = cos(Δ12)
    x3sin12 = x₃^2 * sin(Δ12)
    x4sin12 = x₄^2 * sin(Δ12) / 2
    gLsin1 = gL * sin(x₁)
    gLsin2 = gL * sin(x₂)
    T1_frac = (T₁ - c * x₃) / (2 * mL)
    T2_frac = (T₂ - c * x₄) / mL
    bignum = x3sin12 - ★ * (gLsin1 - x4sin12 + T1_frac) + gLsin2 + T2_frac
    denom = ★^2 / 2 - 1

    dx[1] = x₃
    dx[2] = x₄
    dx[3] = ★ * bignum / (2 * denom) - x4sin12 + gLsin1 + T1_frac
    dx[4] = - bignum / denom
end;

controller = read_nnet_sherlock(@modelpath "" "controller_Double-Pendulum_mr_sherlock")

X₀ = BallInf(fill(1.15, 4), 0.15)
X₀_small = BallInf(fill(1.000005, 4), 0.000005)
U₀ = ZeroSet(2)
vars_idx = Dict(:state_vars=>1:4, :control_vars=>5:6)
ivp = @ivp(x' = double_pendulum!(x), dim: 6, x(0) ∈ X₀ × U₀)
ivp_small = @ivp(x' = double_pendulum!(x), dim: 6, x(0) ∈ X₀_small × U₀)

period = 0.02
k = 20
T = k * period
T_warmup = 2 * period

prob = ControlledPlant(ivp, controller, vars_idx, period)
prob_small = ControlledPlant(ivp_small, controller, vars_idx, period)

# simulation

println("simulation")
res = @timed simulate(prob, T=T, trajectories=10, include_vertices=true)
sim = res.value
print_timed(res)

# reachability analysis

alg = TMJets20(abstol=1e-9, orderT=8, orderQ=1)
alg_nn = Ai2()

function benchmark(; T=T, silent::Bool=false)
    silent || println("flowpipe construction")
    res = @timed solve(prob, T=T, alg_nn=alg_nn, alg=alg)
    sol = res.value
    silent || print_timed(res)

    silent || println("flowpipe construction (small initial states)")
    res = @timed solve(prob_small, T=T, alg_nn=alg_nn, alg=alg)
    sol_small = res.value
    silent || print_timed(res)

    return sol, sol_small
end

benchmark(T=T_warmup, silent=true)
res = @timed benchmark()
sol, sol_small = res.value

# plotting

vars = (1, 2)
fig = plot(xlab=L"\theta", ylab=L"\theta'", leg=:topleft,
		   legendfontsize=5,
		   tickfont=font(25, "Times"),
		   guidefontsize=25,
		   xguidefont=font(25, "Times"),
		   yguidefont=font(25, "Times"),
		   xticks=([1, 1.25, 1.5, 1.75, 2], ["1.0", "1.25", "1.5", "1.75", "2.0"]),
		   bottom_margin=0mm, left_margin=0mm, right_margin=4mm, top_margin=0mm, size=(500, 450))
plot!(fig, sol, vars=vars, color=:orange, lab="")
plot!(project(X₀, vars), c=:blue, lab=L"X_0", alpha=0.6)
plot!(fig, sol_small, vars=vars, color=:red, lab="")
plot!(project(X₀_small, vars), c=:cyan, lab=L"X_0'", alpha=1)
plot_simulation!(fig, sim; vars=vars, color=:black, lab="")
savefig("doublependulum_more_robust_x1x2_juliareach.png")

end
nothing
