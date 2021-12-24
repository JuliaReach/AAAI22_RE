module ACC

using NeuralNetworkAnalysis, LaTeXStrings
using NeuralNetworkAnalysis: FunctionPreprocessing
import Plots, DifferentialEquations
using Plots: plot, plot!, lens!, bbox, savefig, font, Measures.mm, annotations

# problem

u = 0.0001
a_lead = -2.0

@taylorize function ACC!(dx, x, p, t)
    v_lead = x[2]
    γ_lead = x[3]
    v_ego = x[5]
    γ_ego = x[6]
    a_ego = x[7]

    dx[1] = v_lead
    dx[2] = γ_lead
    dx[3] = 2 * (a_lead - γ_lead) - u * v_lead^2
    dx[4] = v_ego
    dx[5] = γ_ego
    dx[6] = 2 * (a_ego - γ_ego) - u * v_ego^2
    dx[7] = zero(a_ego)
    return dx
end

controller = read_nnet_sherlock(@modelpath "" "controller_ACC_sherlock")

X₀ = Hyperrectangle(low= [ 90, 32,   0, 10, 30,   0],
                    high=[110, 32.2, 0, 11, 30.2, 0])
U₀ = ZeroSet(1)
vars_idx = Dict(:state_vars=>1:6, :control_vars=>7)
ivp = @ivp(x' = ACC!(x), dim: 7, x(0) ∈ X₀ × U₀)

v_set = 30.0
T_gap = 1.4
M = zeros(3, 6)
M[1, 5] = 1.0
M[2, 1] = 1.0
M[2, 4] = -1.0
M[3, 2] = 1.0
M[3, 5] = -1.0
function preprocess(X::LazySet)
    Y1 = Singleton([v_set, T_gap])
    Y2 = linear_map(M, X)
    return cartesian_product(Y1, Y2)
end
function preprocess(X::AbstractVector)
    Y1 = [v_set, T_gap]
    Y2 = M * X
    return vcat(Y1, Y2)
end
control_preprocessing = FunctionPreprocessing(preprocess)

period = 0.1
k = 50
T = k * period
T_warmup = 2 * period

prob = ControlledPlant(ivp, controller, vars_idx, period;
                       preprocessing=control_preprocessing)

# simulation

println("simulation")
res = @timed simulate(prob, T=T, trajectories=10, include_vertices=true)
sim = res.value
print_timed(res)

# reachability analysis

alg = TMJets(abstol=1e-6, orderT=6, orderQ=1)
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

vars = (1, 4)
fig = plot(xlab=L"x_{lead}", ylab=L"x_{ego}",
           xlims=(50, 300), ylims=(0, 200),
           leg=:topleft,
           legendfontsize=25,
           tickfont=font(25, "Times"),
           guidefontsize=25,
           xguidefont=font(25, "Times"),
           yguidefont=font(25, "Times"),
           bottom_margin=0mm, left_margin=0mm, right_margin=5mm, top_margin=2mm, size=(500, 450))
plot!(fig, sol, vars=vars, color=:orange, lab="")
plot!(project(X₀, vars), c=:blue, lab=L"X_0", alpha=0.6)
plot_simulation!(fig, sim; vars=vars, color=:black, lab="")
lens!([89, 111], [9.7, 11.5], inset=(1, bbox(0.65, 0.6, 0.3, 0.25)), xtick=90:10:110, ytick=10:1:11, tickfontsize=15, subplot=2)
savefig("acc_juliareach.png")

end
nothing
