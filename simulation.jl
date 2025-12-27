#!/usr/bin/env julia
using Oceananigans
using Oceananigans.Units
using ClimaOcean
using SeawaterPolynomials.TEOS10
using FjordSim
using ArgParse

include("Oxydep.jl")
using .OXYDEPModel

const FT = Oceananigans.defaults.FloatType

# ----------------------------------------------------------
# Command-line argument parser
# ----------------------------------------------------------
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--grid_path"
        help = "Path to the bathymetry NetCDF file."
        arg_type = String
        default = joinpath(homedir(), "FjordSim_data", "oslofjord", "bathymetry_drammen1.nc") #"bathymetry_105to232.nc")

        "--forcing_path"
        help = "Path to the forcing NetCDF file."
        arg_type = String
        default = joinpath(homedir(), "FjordSim_data", "oslofjord", "forcing_drammen.nc") #"forcing_105to232.nc")

        "--atmospheric_forcing_path"
        help = "Path to the atmospheric JRA55 forcing directory."
        arg_type = String
        default = joinpath(homedir(), "FjordSim_data", "JRA55")

        "--results_path"
        help = "Directory where results are stored."
        arg_type = String
        default = joinpath(homedir(), "FjordSim_results", "oslofjord")
    end
    return parse_args(s)
end

# ----------------------------------------------------------
# Main simulation setup
# ----------------------------------------------------------
function main()
    args = parse_commandline()
    println("Running simulation with:")
    println("  grid_path = $(args["grid_path"])")
    println("  forcing_path = $(args["forcing_path"])")
    println("  atmospheric_forcing_path = $(args["atmospheric_forcing_path"])")
    println("  results_path = $(args["results_path"])")

    arch = GPU()
    grid = ImmersedBoundaryGrid(args["grid_path"], arch, (7, 7, 7))
    buoyancy = SeawaterBuoyancy(FT, equation_of_state=TEOS10EquationOfState(FT))
    closure = (
        TKEDissipationVerticalDiffusivity(minimum_tke=7e-6),
        Oceananigans.TurbulenceClosures.HorizontalScalarBiharmonicDiffusivity(ν=15, κ=10),
    )
    # tracer_advection = (T=WENO(), S=WENO(), e=nothing, ϵ=nothing)
    tracer_advection = (
        T=WENO(),
        S=WENO(),
        C=WENO(),
        e=nothing,
        ϵ=nothing,
        NUT=WENO(),
        P=WENO(),
        HET=WENO(),
        POM=WENO(),
        DOM=WENO(),
        O₂=WENO(),
    )
    momentum_advection = WENOVectorInvariant(FT)
    # tracers = (:T, :S, :e, :ϵ)
    tracers = (:T, :S, :e, :ϵ, :C, :NUT, :P, :HET, :POM, :DOM, :O₂)
    # initial_conditions = (T=5.0, S=33.0)
    initial_conditions = (T = 5.0, S = 32.0, C = 0.0, NUT = 10.0, P = 0.05, HET = 0.01, O₂ = 150.0, DOM = 1.0)
    free_surface = SplitExplicitFreeSurface(grid, cfl=0.7)
    coriolis = HydrostaticSphericalCoriolis(FT)
    forcing = forcing_from_file(;
        grid=grid,
        filepath=args["forcing_path"],
        tracers=tracers,
    )
    tbbc = top_bottom_boundary_conditions(;
        grid=grid,
        bottom_drag_coefficient=0.003,
    )
    sobc = (v=(south=OpenBoundaryCondition(nothing),),)
    boundary_conditions = map(x -> FieldBoundaryConditions(; x...), recursive_merge(tbbc, sobc))
    # biogeochemistry = nothing
    biogeochemistry = OXYDEP(grid)
    boundary_conditions = merge(boundary_conditions, bgh_oxydep_boundary_conditions(biogeochemistry, grid.Nz))
    atmosphere = JRA55PrescribedAtmosphere(arch, FT;
        latitude=(58.98, 59.94),
        longitude=(10.18, 11.03),
        dir=args["atmospheric_forcing_path"],
    )
    downwelling_radiation = Radiation(arch, FT;
        ocean_emissivity=0.96,
        ocean_albedo=0.1
    )
    sea_ice = FreezingLimitedOceanTemperature()
    results_dir = args["results_path"]
    stop_time = 365days

    simulation = coupled_hydrostatic_simulation(
        grid,
        buoyancy,
        closure,
        tracer_advection,
        momentum_advection,
        tracers,
        initial_conditions,
        free_surface,
        coriolis,
        forcing,
        boundary_conditions,
        atmosphere,
        downwelling_radiation,
        sea_ice,
        biogeochemistry;
        results_dir,
        stop_time,
    )

    simulation.callbacks[:progress] = Callback(progress, TimeInterval(12hours))

    ocean_sim = simulation.model.ocean
    ocean_model = ocean_sim.model

    prefix = joinpath(results_dir, "snapshots_ocean")
    ocean_sim.output_writers[:ocean] = NetCDFWriter(
        ocean_model,
        (
            T=ocean_model.tracers.T,
            S=ocean_model.tracers.S,
            NUT=ocean_model.tracers.NUT,
            P=ocean_model.tracers.P,
            HET=ocean_model.tracers.HET,
            POM=ocean_model.tracers.POM,
            DOM=ocean_model.tracers.DOM,
            O₂=ocean_model.tracers.O₂,
            u=ocean_model.velocities.u,
            v=ocean_model.velocities.v,
        );
        filename="$prefix",
        schedule=TimeInterval(6hours),
        overwrite_existing=true,
    )

    conjure_time_step_wizard!(simulation; cfl=0.1, max_Δt=3minutes, max_change=1.01)
    run!(simulation)
end

# ----------------------------------------------------------
# Run script
# ----------------------------------------------------------
main()
