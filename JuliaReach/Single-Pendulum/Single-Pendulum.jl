module SinglePendulum

using NeuralNetworkAnalysis, LaTeXStrings
import Plots, DifferentialEquations
using Plots: plot, plot!, savefig, font, Measures.mm

# problem

m = 0.5
L = 0.5
c = 0.0
g = 1.0
gL = g/L
mL = 1/(m*L^2)

@taylorize function single_pendulum!(dx, x, p, t)
    dx[1] = x[2]
    dx[2] = gL * sin(x[1]) + mL*(x[3] - c*x[2])
    dx[3] = zero(x[3])
    return dx
end

controller = read_nnet_sherlock(@modelpath "" "controller_Single-Pendulum_sherlock")

X0 = Hyperrectangle([1.1, 0.1], [0.1, 0.1])
U0 = ZeroSet(1)
ivp = @ivp(x' = single_pendulum!(x), dim: 3, x(0) ∈ X0 × U0)
vars_idx = Dict(:state_vars=>1:2, :control_vars=>3)

period = 0.05
k = 20
T = k * period
T_warmup = 2 * period

prob = ControlledPlant(ivp, controller, vars_idx, period)

# simulation

println("simulation")
res = @timed simulate(prob, T=T, trajectories=10, include_vertices=true)
sim = res.value
print_timed(res)

# reachability analysis

alg = TMJets(abstol=1e-7, orderT=4, orderQ=1)
alg_nn = Ai2()

function benchmark(; T=T, silent::Bool=false)
    silent || println("flowpipe construction")
    res = @timed solve(prob, T=T, alg_nn=alg_nn, alg=alg)
    sol = res.value
    silent || print_timed(res)

    return sol
end

benchmark(T=T_warmup, silent=true)
res = @timed benchmark()
sol = res.value

# plotting

vars = (1, 2)
fig = plot(xlab=L"\theta", ylab=L"\theta'",
           xlims=(1, 2.6), ylims=(-0.5, 2.5),
           leg=:topleft,
           legendfontsize=5,
           tickfont=font(25, "Times"),
           guidefontsize=25,
           xguidefont=font(25, "Times"),
           yguidefont=font(25, "Times"),
           bottom_margin=0mm, left_margin=0mm, right_margin=0mm, top_margin=2mm, size=(500, 450))
plot!(fig, sol, vars=vars, color=:orange, lab="")
plot!(project(X0, vars), c=:blue, lab=L"X_0", alpha=0.6)
plot_simulation!(fig, sim; vars=vars, color=:black, lab="")
savefig("singlependulum_juliareach.png")

end
nothing
