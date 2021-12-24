module Airplane

using NeuralNetworkAnalysis, LaTeXStrings
import Plots, DifferentialEquations
using Plots: plot, plot!, savefig, font, Measures.mm

# problem

Tψ = ψ -> [cos(ψ)  -sin(ψ)  0;
           sin(ψ)   cos(ψ)  0;
                0        0  1]

Tθ = θ -> [ cos(θ)  0  sin(θ);
                 0  1       0;
           -sin(θ)  0  cos(θ)]

Tϕ = ϕ -> [1       0        0;
           0  cos(ϕ)  -sin(ϕ);
           0  sin(ϕ)   cos(ϕ)]

Rϕθ = (ϕ, θ) -> [1  tan(θ) * sin(ϕ)  tan(θ) * cos(ϕ);
                 0           cos(θ)          -sin(ϕ);
                 0  sec(θ) * sin(ϕ)  sec(θ) * cos(ϕ)]

m = 1.0
g = 1.0
Ix = 1.0
Iy = 1.0
Iz = 1.0
Ixz = 0.0

@taylorize function airplane!(dx, x, p, t)
    _x, y, z, u, v, w, ϕ, θ, ψ, r, _p, q, Fx, Fy, Fz, Mx, My, Mz = x

    T_ψ = Tψ(ψ)
    T_θ = Tθ(θ)
    T_ϕ = Tϕ(ϕ)
    mat_1 = T_ψ * T_θ * T_ϕ
    xyz = mat_1 * vcat(u, v, w)

    mat_2 = Rϕθ(ϕ, θ)
    ϕθψ = mat_2 * vcat(_p, q, r)

    dx[1] = xyz[1]
    dx[2] = xyz[2]
    dx[3] = xyz[3]
    dx[4] = -g * sin(θ) + Fx / m - q * w + r * v
    dx[5] = g * cos(θ) * sin(ϕ) + Fy / m - r * u + _p * w
    dx[6] = g * cos(θ) * cos(ϕ) + Fz / m - _p * v + q * u
    dx[7] = ϕθψ[1]
    dx[8] = ϕθψ[2]
    dx[9] = ϕθψ[3]
    dx[10] = Mz
    dx[11] = My
    dx[12] = Mx
    dx[13] = zero(Fx)
    dx[14] = zero(Fy)
    dx[15] = zero(Fz)
    dx[16] = zero(Mx)
    dx[17] = zero(My)
    dx[18] = zero(Mz)
end

controller = read_nnet_sherlock(@modelpath "" "controller_Airplane_sherlock")

X₀ = Hyperrectangle(low=[ 0,       0,       0,       0.99999, 0.99999, 0.99999, 0.99999, 0.99999, 0.99999, 0,       0,       0],
                    high=[0.00001, 0.00001, 0.00001, 1,       1,       1,       1,       1,       1,       0.00001, 0.00001, 0.00001])
U₀ = ZeroSet(6)
vars_idx = Dict(:state_vars=>1:12, :control_vars=>13:18)
ivp = @ivp(x' = airplane!(x), dim: 18, x(0) ∈ X₀ × U₀)

period = 0.1
k = 20
T = k * period
T_warmup = 2 * period

prob = ControlledPlant(ivp, controller, vars_idx, period)

# simulation

println("simulation")
res = @timed simulate(prob, T=T, trajectories=3, include_vertices=false)
sim = res.value
print_timed(res)

# reachability analysis

alg = TMJets(abstol=1e-10, orderT=7, orderQ=1)
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

vars = (2, 7)
fig = plot(xlab=L"y", ylab=L"\phi", leg=:topleft,
           legendfontsize=5,
           tickfont=font(25, "Times"),
           guidefontsize=25,
           xguidefont=font(25, "Times"),
           yguidefont=font(25, "Times"),
           bottom_margin=0mm, left_margin=0mm, right_margin=0mm, top_margin=0mm, size=(500, 450))
plot!(fig, sol, vars=vars, color=:orange, lab="")
plot!(project(X₀, vars), c=:blue, lab=L"X_0", alpha=0.6)
plot_simulation!(fig, sim; vars=vars, color=:black, lab="")
savefig("airplane_yphi_juliareach.png")

end
nothing
