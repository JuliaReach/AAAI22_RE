module Unicycle

using NeuralNetworkAnalysis, LaTeXStrings
using NeuralNetworkAnalysis: UniformAdditivePostprocessing
import Plots, DifferentialEquations
using Plots: plot, plot!, scatter!, lens!, bbox, savefig, font, Measures.mm

# problem

@taylorize function unicycle!(dx, x, p, t)
    x₁, x₂, x₃, x₄, w, u₁, u₂ = x

    dx[1] = x₄ * cos(x₃)
    dx[2] = x₄ * sin(x₃)
    dx[3] = u₂
    dx[4] = u₁ + w
    dx[5] = zero(x[5])
    dx[6] = zero(x[6])
    dx[7] = zero(x[7])
    return dx
end

controller = read_nnet_sherlock(@modelpath "" "controller_Unicycle_sherlock")

X₀ = Hyperrectangle(low=[ 9.5,  -4.5,  2.1,  1.5, -1e-4],
                    high=[9.55, -4.45, 2.11, 1.51, 1e-4])
U₀ = ZeroSet(2)
vars_idx = Dict(:state_vars=>1:4, :input_vars=>[5], :control_vars=>6:7)
ivp = @ivp(x' = unicycle!(x), dim: 7, x(0) ∈ X₀ × U₀)

control_postprocessing = UniformAdditivePostprocessing(-20.0)

period = 0.2
T = 10.0
T_warmup = 2 * period

prob = ControlledPlant(ivp, controller, vars_idx, period;
                       postprocessing=control_postprocessing)

target_states = cartesian_product(Hyperrectangle(zeros(4), [0.6, 0.2, 0.06, 0.3]),
                                  Universe(3))

predicate_sol_suff = sol -> sol[end][end] ⊆ target_states

# simulation

println("simulation")

res = @timed simulate(prob, T=T; trajectories=1, include_vertices=false)
sim1 = res.value
print_timed(res)

res = @timed simulate(prob, T=T; trajectories=10, include_vertices=true)
sim = res.value
print_timed(res)

# reachability analysis

alg = TMJets(abstol=1e-15, orderT=10, orderQ=1)
alg_nn = Ai2()
splitter = BoxSplitter([3, 1, 8, 1])

function benchmark(; T=T, silent::Bool=false)
    silent || println("flowpipe construction")
    res = @timed solve(prob, T=T, alg_nn=alg_nn, alg=alg,
                                 splitter=splitter)
    sol = res.value
    silent || print_timed(res)

    silent || println("property checking")
    res = @timed predicate_sol_suff(sol)
    silent || print_timed(res)
    if res.value
        silent || println("The property is satisfied.")
    else
        silent || println("The property may be violated.")
    end
    return sol
end

benchmark(T=T_warmup, silent=true)
res = @timed benchmark()
sol = res.value
println("total analysis time")
print_timed(res)

# plotting

function plot_sim!(p, sim, vars)
    plot_simulation!(p, sim; vars=vars, color=:black, lab="")
end

function plot_sol!(p, sol, vars)
    plot!(p, sol, vars=vars, c=:orange, alpha=0.1)
    plot!(p, sol[end][end], vars=vars, c=:red, alpha=0.3)
end

function plot_zoom!(p, vars)
    if vars == [1, 2]
        lens!(p, [9.45, 9.6], [-4.55, -4.4], inset=(1, bbox(0.7, 0.2, 0.3, 0.25)), xticks=[9.5, 9.55], yticks=[-4.5, -4.45], subplot=2)
        lens!(p, [0.3, 0.7], [-0.3, 0.3], inset=(1, bbox(0.1, 0.3, 0.3, 0.3)), xticks=[0.4, 0.6], subplot=3)
    elseif vars == [3, 4]
        plot!(xlims=[-0.2, 2.7], ylims=[-1.5, 3], leg=:topleft)
        lens!(p, [2.095, 2.12], [1.49, 1.52], inset=(1, bbox(0.4, 0.1, 0.3, 0.3)), xticks=[2.10, 2.11, 2.12], subplot=2)
        lens!(p, [-0.1, 0.1], [-0.5, -0.1], inset=(1, bbox(0.6, 0.6, 0.3, 0.3)), xticks=[-0.1, 0, 0.1], subplot=3)
    end
end

function plot_X0!(p, vars; lab=L"X_0")
    plot!(p, project(X₀, vars), c=:blue, lab=lab, alpha=0.6)
end

function label(var)
    var == 1 && return L"x"
    var == 2 && return L"y"
    var == 3 && return L"\theta"
    var == 4 && return L"v"
    throw(ArgumentError(var))
end

function plot_all(sol, sim, vars)
    p = plot(xlab=label(vars[1]), ylab=label(vars[2]),
             tickfont=font(15, "Times"),
             guidefontsize=15,
             xguidefont=font(15, "Times"),
             yguidefont=font(15, "Times"),
             bottom_margin=1mm, left_margin=0mm, right_margin=0mm, top_margin=0mm)
    plot_X0!(p, vars)
    plot!(p, project(target_states, vars), c=:green, lab=L"X_T")

    !isnothing(sol) && plot_sol!(p, sol, vars)
    !isnothing(sol) && plot_X0!(p, vars, lab="")
    plot_sim!(p, sim, vars)
    plot_zoom!(p, vars)
    return p
end

plot_all(nothing, sim, [1, 2])
savefig("unicycle_sim_x1x2.png")

plot_all(nothing, sim, [3, 4])
savefig("unicycle_sim_x3x4.png")

plot_all(sol, sim1, [1, 2])
savefig("unicycle_reach_x1x2.png")

plot_all(sol, sim1, [3, 4])
savefig("unicycle_reach_x3x4.png")

tdom = range(0, T - period, step=period)
fig = plot(xlab=L"t", ylab=L"u_1(t) / u_2(t)",
           leg=:topleft,
           legendfontsize=5,
           tickfont=font(15, "Times"),
           guidefontsize=15,
           xguidefont=font(15, "Times"),
           yguidefont=font(15, "Times"),
           bottom_margin=2mm, left_margin=0mm, right_margin=5mm, top_margin=1mm)
[plot!(fig, tdom, [c[1] for c in controls(sim, i)], lab=(i == 1 ? L"u_1" : ""), c=:red, seriestype=:steppre) for i in 1:length(sim)]
[plot!(fig, tdom, [c[2] for c in controls(sim, i)], lab=(i == 1 ? L"u_2" : ""), c=:blue, seriestype=:steppre) for i in 1:length(sim)]
plot!()
savefig("unicycle_controls.png")

end
nothing
