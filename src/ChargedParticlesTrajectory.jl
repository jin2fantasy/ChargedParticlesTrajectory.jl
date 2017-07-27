module ChargedParticlesTrajectory

# package code goes here
using DifferentialEquations,
      StaticArrays,
      Plots


export solve_electrostatics,
       track_particles

abstract type AbstractSuperParticle end

# super electron.
struct SuperElectron{T<:AbstractFloat}
    position::MVector{3,T}
    momentum::MVector{3,T}
    charge::T
    mass::T
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
    return sol
end

# solve the equation of motion
# and update the position and momentum
# ṗ = E⃗
# ẋ = pc² / √( (m₀c²)² + (pc)² )
#   = 1 / √( (m₀/p)² + (1/c)² )
function track_particles()
    1
end

end # module
potential = solve_electrostatics(x->0, x-> (x[:,1] .≈ 0) .+ (x[:,1] .≈ 1), [0 1 0 1], 0.05)
plot(potential)
potential.u[1][220]

function ϕ_interpolated(x, y, ϕ, domain)
    m = n = round(Int, sqrt(length(ϕ[1])))
    xmin, xmax, ymin, ymax = domain[1], domain[2], domain[3], domain[4]
    xlength, ylength = xmax-xmin, ymax-ymin
    xstep, ystep = xlength/m, ylength/n
    x_lowind = floor(Int, x / xstep)
    x_highind = ceil(Int, x / xstep)
    y_lowind = floor(Int, y / ystep)
    y_highind = ceil(Int, y / ystep)
    x_low = xstep * x_lowind
    x_high = xstep * x_highind
    y_low = ystep * y_lowind
    y_high = ystep * y_highind
    cell_area = (x_high - x_low)*(y_high - y_low)
    ϕ_lowleft = (x - x_low)*(y - y_low)/cell_area * ϕ[1][sub2ind((m,n), x_lowind, y_lowind)]
    ϕ_upleft = (x - x_low)*(y_high - y)/cell_area * ϕ[1][sub2ind((m,n), x_lowind, y_highind)]
    ϕ_lowright = (x_high - x)*(y - y_low)/cell_area * ϕ[1][sub2ind((m,n), x_highind, y_lowind)]
    ϕ_upright = (x_high - x)*(y_high - y)/cell_area * ϕ[1][sub2ind((m,n), x_highind, y_lowind)]
    return ϕ_upright + ϕ_lowright + ϕ_upleft + ϕ_lowleft
end
Efield(x) = -ForwardDiff.gradient(x -> ϕ_interpolated(x[1], x[2], potential, [0 1 0 1]), x)
Efield([0.01, 0.5])
