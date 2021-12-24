module TORA

using NeuralNetworkAnalysis, LaTeXStrings
using NeuralNetworkAnalysis: UniformAdditivePostprocessing
import Plots, DifferentialEquations
using Plots: plot, plot!, xlims!, ylims!, savefig, font, Measures.mm

# problem

@taylorize function TORA!(dx, x, p, t)
    x₁, x₂, x₃, x₄, u = x

    aux = 0.1 * sin(x₃)
    dx[1] = x₂
    dx[2] = -x₁ + aux
    dx[3] = x₄
    dx[4] = u
    dx[5] = zero(u)
    return dx
end

controller = read_nnet_sherlock(@modelpath "" "controller_TORA_sherlock")

X₀ = Hyperrectangle(low=[0.6, -0.7, -0.4, 0.5], high=[0.7, -0.6, -0.3, 0.6])
U = ZeroSet(1)
vars_idx = Dict(:state_vars=>1:4, :control_vars=>5)
ivp = @ivp(x' = TORA!(x), dim: 5, x(0) ∈ X₀ × U)
control_postprocessing = UniformAdditivePostprocessing(-10.0)

period = 1.0
k = 20
T = k * period
T_warmup = 2 * period

prob = ControlledPlant(ivp, controller, vars_idx, period;
                       postprocessing=control_postprocessing)

safe_states = cartesian_product(BallInf(zeros(4), 2.0), Universe(1))
predicate = X -> X ⊆ safe_states

# simulation

println("simulation")
res = @timed simulate(prob, T=T, trajectories=10, include_vertices=true)
sim = res.value
print_timed(res)

# reachability analysis

alg = TMJets(abstol=1e-10, orderT=8, orderQ=3)
alg_nn = Ai2()
splitter = BoxSplitter([4, 4, 3, 5])

function benchmark(; T=T, silent::Bool=false)
    silent || println("flowpipe construction")
    res = @timed solve(prob, T=T, alg_nn=alg_nn, alg=alg, splitter=splitter)
    sol = res.value
    silent || print_timed(res)

    silent || println("property checking")
    solh = overapproximate(sol, Hyperrectangle)
    res = @timed predicate(solh)
    silent || print_timed(res)
    if res.value
        silent || println("The property is satisfied.")
    else
        silent || println("The property may be violated.")
    end

    return solh
end

benchmark(T=T_warmup, silent=true)
res = @timed benchmark()
sol = res.value

# plotting

vars = (1, 2)
fig = plot(xlab=L"x_1", ylab=L"x_2",
           legendfontsize=15,
           tickfont=font(25, "Times"),
           guidefontsize=25,
           xguidefont=font(25, "Times"),
           yguidefont=font(25, "Times"),
           bottom_margin=0mm, left_margin=0mm, right_margin=3mm, top_margin=2mm,
           size=(500, 450))
plot!(project(X₀, vars), c=:blue, lab=L"X_0", alpha=0.6)
plot!(fig, sol, vars=vars, color=:orange, lab="")
plot!(project(X₀, vars), c=:blue, lab="", alpha=0.6)
plot_simulation!(fig, sim; vars=vars, color=:black, lab="")
xlims!(-1, 1)
ylims!(-1.5, 1)
savefig("tora_x1x2_juliareach.png")

vars = (3, 4)
fig = plot(xlab=L"x_3", ylab=L"x_4",
           legendfontsize=15,
           tickfont=font(25, "Times"),
           guidefontsize=25,
           xguidefont=font(25, "Times"),
           yguidefont=font(25, "Times"),
           bottom_margin=0mm, left_margin=0mm, right_margin=0mm, top_margin=2mm,
           size=(500, 450))
# plot!(project(X₀, vars), c=:blue, lab=L"X_0", alpha=0.6)
plot!(fig, sol, vars=vars, color=:orange, lab="")
plot!(project(X₀, vars), c=:blue, lab="", alpha=0.6)
plot_simulation!(fig, sim; vars=vars, color=:black, lab="")
xlims!(-2, 2)
ylims!(-2, 2)
savefig("tora_x3x4_juliareach.png")

end
nothing
