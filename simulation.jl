#!/usr/bin/env julia
using Oceananigans
using Oceananigans.Units
using ClimaOcean
using SeawaterPolynomials.TEOS10
using FjordSim
using FjordSim.FDatasets
using ArgParse
using Dates: DateTime

if !@isdefined(Ci_)
    const Ci_ = false  # true: include Ci fields and calculations; false: skip them
else
    global Ci_ = false
end

include("Oxydep.jl")
using .OXYDEPModel

const _uniform_IC = false # if false, read initial conditions from "snapshots_ocean.nc"

const FT = Oceananigans.defaults.FloatType

function select_arch(arch_str::AbstractString)
    arch = lowercase(strip(arch_str))

    if arch == "cpu"
        return CPU()
    elseif arch == "gpu"
        return GPU()  # throws if no CUDA GPU
    else
        error("Invalid --arch value: $(repr(arch_str)). Use 'cpu' or 'gpu'.")
    end
end


# ----------------------------------------------------------
# Command-line argument parser
# ----------------------------------------------------------
function default_paths()
    base_dir = @__DIR__
    grid     = joinpath(base_dir, "data", "input", "bathymetry_19to67.nc")
    forcing  = joinpath(base_dir, "data", "input", "forcing_drammen6.nc")
    atm      = joinpath(base_dir, "data", "input", "JRA55")
    results  = joinpath(base_dir, "data", "output", "drammensfjord")

    return (; grid, forcing, atm, results)
end

function parse_commandline()
    s = ArgParseSettings()
    paths = default_paths()

    @add_arg_table! s begin
        "--grid_path"
        help = "Path to the bathymetry NetCDF file."
        arg_type = String
        default = paths.grid

        "--forcing_path"
        help = "Path to the forcing NetCDF file."
        arg_type = String
        default = paths.forcing

        "--atmospheric_forcing_path"
        help = "Path to the atmospheric forcing directory."
        arg_type = String
        default = paths.atm

        "--results_path"
        help = "Directory where results are stored."
        arg_type = String
        default = paths.results

        "--arch"
        help = "Compute architecture: cpu or gpu."
        arg_type = String
        default = "gpu"

    end
    return parse_args(s)
end


# ----------------------------------------------------------
# Main simulation setup
# ----------------------------------------------------------
function main()
    args = parse_commandline()
    arch = select_arch(args["arch"])

    println("Running simulation with:")
    println("  arch = $(args["arch"])")
    println("  selected_architecture = $(typeof(arch))")
    println("  grid_path = $(args["grid_path"])")
    println("  forcing_path = $(args["forcing_path"])")
    println("  atmospheric_forcing_path = $(args["atmospheric_forcing_path"])")
    println("  results_path = $(args["results_path"])")

    grid = ImmersedBoundaryGrid(args["grid_path"], arch, (7, 7, 7))
    buoyancy = SeawaterBuoyancy(FT, equation_of_state=TEOS10EquationOfState(FT))
#¤    closure = (
#¤        TKEDissipationVerticalDiffusivity(minimum_tke=7e-6),
#¤        Oceananigans.TurbulenceClosures.HorizontalScalarBiharmonicDiffusivity(ν=15, κ=10),
#¤   )
#¤    closure = VerticalScalarDiffusivity(ν=1e-3, κ=6e-7) # κ is for tracers, ν is for momentum
#¤    tracer_advection = WENO()
    closure = (
        CATKEVerticalDiffusivity(minimum_tke = 7e-6), #% (minimum_tke = 7e-6)
        Oceananigans.TurbulenceClosures.HorizontalScalarBiharmonicDiffusivity(ν = 15, κ = 10),
     )    
#¤     tracer_advection = ()
    tracer_advection = if Ci_
        (
            T=WENO(),
            S=WENO(),
            e=nothing,
            NUT=WENO(),
            P=WENO(),
            HET=WENO(),
            POM=WENO(),
            DOM=WENO(),
            O₂=WENO(),
            Ci_free=WENO(),   
            Ci_PHY=WENO(),        
            Ci_HET=WENO(),   
            Ci_POM=WENO(),   
            Ci_DOM=WENO(), 
        )
    else
        (
            T=WENO(),
            S=WENO(),
            e=nothing,
            NUT=WENO(),
            P=WENO(),
            HET=WENO(),
            POM=WENO(),
            DOM=WENO(),
            O₂=WENO(),
        )
    end
    momentum_advection = WENOVectorInvariant(FT)

    tracers = if Ci_
        (
            :T, :S, :e, :NUT, :P, :HET, :POM, :DOM, :O₂,
            :Ci_free, :Ci_PHY, :Ci_HET, :Ci_POM, :Ci_DOM,
        )
    else
        (:T, :S, :e, :NUT, :P, :HET, :POM, :DOM, :O₂)
    end
    
    if !_uniform_IC
        restart_paths = default_paths()
        dataset = DSResults(
            "snapshots_ocean_6.nc",
            restart_paths.results;
            start_date_time = DateTime(2025, 1, 1),
        )
        initial_conditions_base = (
            T = Metadatum(:temperature; dataset, date = last_date(dataset, :temperature)),
            S = Metadatum(:salinity; dataset, date = last_date(dataset, :salinity)),
            e = Metadatum(:e; dataset, date = last_date(dataset, :e)),
            NUT = Metadatum(:NUT; dataset, date = last_date(dataset, :NUT)),
            P = Metadatum(:P; dataset, date = last_date(dataset, :P)),
            HET = Metadatum(:HET; dataset, date = last_date(dataset, :HET)),
            POM = Metadatum(:POM; dataset, date = last_date(dataset, :POM)),
            DOM = Metadatum(:DOM; dataset, date = last_date(dataset, :DOM)),
            O₂ = Metadatum(:O₂; dataset, date = last_date(dataset, :O₂)),
            u = Metadatum(
                :u_velocity;
                dataset,
                date = last_date(dataset, :u_velocity),
            ),
            v = Metadatum(
                :v_velocity;
                dataset,
                date = last_date(dataset, :v_velocity),
            ),
        )
        if Ci_
            initial_conditions = merge(initial_conditions_base, (
                Ci_free = Metadatum(:Ci_free; dataset, date = last_date(dataset, :Ci_free)),
                Ci_PHY = Metadatum(:Ci_PHY; dataset, date = last_date(dataset, :Ci_PHY)),
                Ci_HET = Metadatum(:Ci_HET; dataset, date = last_date(dataset, :Ci_HET)),
                Ci_POM = Metadatum(:Ci_POM; dataset, date = last_date(dataset, :Ci_POM)),
                Ci_DOM = Metadatum(:Ci_DOM; dataset, date = last_date(dataset, :Ci_DOM)),
            ))
        else
            initial_conditions = initial_conditions_base
        end
    else
        initial_conditions_base = (
            T = 1.0, S = 29.3, NUT = 18., P = 0.01, HET = 0.01, O₂ = 150., DOM = 0.05, POM = 0.01,
        )
        if Ci_
            initial_conditions = merge(initial_conditions_base, (
                Ci_free = 0.0, Ci_PHY = 0.0, Ci_HET = 0.0, Ci_POM = 0.0, Ci_DOM = 0.0,
            ))
        else
            initial_conditions = initial_conditions_base
        end
    end
    free_surface = SplitExplicitFreeSurface(grid, cfl = 0.7)
    coriolis = HydrostaticSphericalCoriolis(FT)
    forcing = forcing_from_file(;
        grid=grid,
        filepath=args["forcing_path"],
        tracers=tracers,
    )
    tbbc = top_bottom_boundary_conditions(;
        grid = grid,
        bottom_drag_coefficient = 0.003,
    )
    sobc = (v = (south = OpenBoundaryCondition(nothing),),)
    boundary_conditions = map(x -> FieldBoundaryConditions(; x...), recursive_merge(tbbc, sobc))
    # biogeochemistry disabled
    biogeochemistry = OXYDEPModel.OXYDEP(grid)
    boundary_conditions = merge(boundary_conditions, OXYDEPModel.bgh_oxydep_boundary_conditions(biogeochemistry, grid.Nz))    
    atmosphere = JRA55PrescribedAtmosphere(arch, FT;
         latitude=(58.98, 59.94),
         longitude=(10.18, 11.03),
         dir=args["atmospheric_forcing_path"],
     )
#    atmosphere = NORA3PrescribedAtmosphere(arch)
    downwelling_radiation = Radiation(arch, FT;
        ocean_emissivity=0.96,
        ocean_albedo=0.1
    )
    sea_ice = FreezingLimitedOceanTemperature()
    results_dir = args["results_path"]
    stop_time = 365days
    #¤biogeochemistry = nothing

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

    simulation.callbacks[:progress] = Callback(progress, TimeInterval(24hours))

    ocean_sim = simulation.model.ocean
    ocean_model = ocean_sim.model

    prefix = joinpath(results_dir, "snapshots_ocean")
    output_fields_base = (
        T=ocean_model.tracers.T,
        S=ocean_model.tracers.S,
        NUT=ocean_model.tracers.NUT,
        P=ocean_model.tracers.P,
        HET=ocean_model.tracers.HET,
        POM=ocean_model.tracers.POM,
        DOM=ocean_model.tracers.DOM,
        O₂=ocean_model.tracers.O₂,
        e=ocean_model.tracers.e,
        u=ocean_model.velocities.u,
        v=ocean_model.velocities.v,
    )
    output_fields = if Ci_
        merge(output_fields_base, (
            Ci_free=ocean_model.tracers.Ci_free,
            Ci_PHY=ocean_model.tracers.Ci_PHY,
            Ci_HET=ocean_model.tracers.Ci_HET,
            Ci_POM=ocean_model.tracers.Ci_POM,
            Ci_DOM=ocean_model.tracers.Ci_DOM,
        ))
    else
        output_fields_base
    end
    ocean_sim.output_writers[:ocean] = NetCDFWriter(
        ocean_model,
        output_fields;
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
