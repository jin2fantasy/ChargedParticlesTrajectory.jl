module ChargedParticlesTrajectory

# package code goes here
using DifferentialEquations
using StaticArrays
using Plots


export solve_electrostatics
export track_particles

abstract type AbstractSuperParticle end

# super electron.
struct SuperElectron{T}
    pos::MVector{3,T}
    charge
    mass
end

# Δϕ = ρ/ϵ₀
function solve_electrostatics(source_func, boundary_cond, domain, dx)
    # dx = 0.05
    # mesh = notime_squaremesh([0 1 0 1], dx, :dirichlet)
    # source_func(x) = 0.0
    # boundary_cond(x) = (x[:,1].≈0) .* 0 .+ (x[:,1].≈1) .* 1
    mesh = notime_squaremesh(domain, dx, :dirichlet)
    prob = PoissonProblem(source_func, mesh, gD=boundary_cond)
    sol = solve(prob)
    plot(sol)
end

function track_particles()

end

end # module
