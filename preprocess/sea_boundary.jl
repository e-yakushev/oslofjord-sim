using XLSX
using CSV
using DataFrames
using Dates
using Statistics
using Tables
using Dierckx
using Interpolations
using Impute
using CairoMakie

# ---------------------------
# Config
# ---------------------------
base_dir = dirname(@__DIR__)
folder = joinpath(base_dir, "data", "input", "Sea_boundary")
#filename = joinpath(folder, "VT10_for_Oxydep_2000_before.xlsx")
filename = joinpath(folder, "Im2_for_OxyDep_old_13.xlsx") # worked 260404
#filename = joinpath(folder, "Drammen", "Drammensfjord_inside.xlsx") # 
interp_method = "linear"  # Options: "spline", "linear", "fill"

# ---------------------------
# Load Excel
# ---------------------------
sheet = XLSX.gettable(XLSX.readxlsx(filename)["WaterChemistry"]) |> DataFrame

rename!(sheet, Dict(
    :date   => :date,
    :depth1 => :depth,
    :o2     => :o2,
    :no3    => :no3,
    :nh4    => :nh4,
    :temp   => :temp,
    :salt   => :salt,
    :chla   => :chla,
))

if eltype(sheet.date) <: AbstractString
    sheet.date = Date.(sheet.date, "dd.mm.yyyy")
end

function row_has_lt(row, cols)
    @inbounds for c in cols
        v = row[c]
        if v isa AbstractString && startswith(v, "<")
            return true
        end
    end
    return false
end

#param_cols = [:no3, :nh4, :chla]
#param_cols = [:o2, :no3, :temp, :salt, :chla]
param_cols = [:o2, :no3, :nh4, :temp, :salt, :chla]
#param_cols = [:o2, :no3, :nh4]
mask = [!row_has_lt(row, param_cols) for row in eachrow(sheet)]
sheet = sheet[mask, :]
sheet = dropmissing(sheet, [:depth, :date])

function to_date(datestr)
    if datestr isa Date
        return datestr
    elseif datestr isa DateTime
        return Date(datestr)
    elseif datestr isa AbstractString && length(datestr) >= 10
        return Date(datestr[1:10], dateformat"yyyy-mm-dd")
    else
        error("Unrecognized date format: $datestr")
    end
end

dates = [to_date(d) for d in sheet.date]

println("Parsed dates: ", dates[1:5], " ... ", dates[end-4:end])
# ---------------------------
# Build time axis from first to last observation (day numbers)
# ---------------------------
#first_date = minimum(sheet.date)
#last_date  = maximum(sheet.date)
first_date = minimum(dates)
last_date  = maximum(dates)
ndays      = Int(Dates.value(last_date - first_date)) + 1           # inclusive span
dayindex(d::Date) = Int(Dates.value(d - first_date)) + 1            # 1..ndays for any observed date
dayindex(d) = dayindex(Date(d))

# ---------------------------
# Fill if needed
# ---------------------------
function apply_forward_backward_fill!(df::DataFrame, cols::Vector{Symbol})
    for col in cols
        col_data = df[!, col]
        if eltype(col_data) !== Union{Missing, Float64}
            try
                col_data = convert(Vector{Union{Missing, Float64}}, col_data)
            catch e
                @warn "Skipping $col: cannot convert to Float64 with missing." exception=(e, catch_backtrace())
                continue
            end
        end
        col_data = Impute.locf(col_data)
        col_data = Impute.nocb(col_data)
        df[!, col] = col_data
    end
end

if interp_method == "fill"
    apply_forward_backward_fill!(sheet, param_cols)
end

# ---------------------------
# Transformations
# --------------------------- here names are as in input file
param_transforms = Dict(
    :o2   => (x -> x * 44.88), # mol O2/L -> mg O2/m^3 (example factor from earlier)
    :no3  => (x -> x / 14),
    :nh4  => (x -> x / 14),
    :temp => identity,
    :salt => identity,
    :phy  => (x -> x / 2.) 
)

param_source_column = Dict(
    :o2    => :o2,
    :no3   => :no3,
    :nh4   => :nh4,
    :temp  => :temp,
    :salt  => :salt,
    :phy   => :chla
)

param_names = Dict(
    :o2   => "O2",
    :no3  => "NUT",
    :nh4  => "DOM",
    :temp => "TEMP",
    :salt => "SALT",
    :phy => "P",
)

parameters = collect(keys(param_transforms))

# ---------------------------
# Interpolation at specific depths
# ---------------------------
function interpolate_param(df::DataFrame, source_col::Symbol, convert_fn::Function,
                           ndays::Int, dayindex::Function)
#function interpolate_param(df::DataFrame, source_col::Symbol, convert_fn::Function)
                          days = 1:ndays # days = 1:365 #
# target_depths = [1, 1.5, 2.5, 4, 6.25, 8.75, 12.5, 20, 37.5, 62.5, 87.5, 125]
# target_depths = [0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 15.0, 20.0, 40.0, 60.0, 80.0, 90.0, 100.0, 110.0]
# target_depths = [0.5, 1.5, 2.5, 4., 6.25, 8.75, 12.5, 17.5 , 30., 50., 70., 85., 95., 105., 115, 125, 135, 150, 170, 195]     
  target_depths = [0.5, 1.5, 2.5, 4., 6.25, 8.75, 12.5, 17.5, 25, 35., 45., 55., 65., 75., 85., 95, 107.5] #Drammensfjord     
# target_depths = [2, 4, 8, 12, 16, 20, 30, 40, 50, 60, 80, 100, 125, 150, 195]
    interpolated = DataFrame(depth=Float64[], day=Int[], value=Float64[])

    for d in target_depths
        sub = df[df.depth .== d, :]
        raw_vec = sub[!, source_col]

        vec = Union{Missing, Float64}[x isa Missing ? missing :
            x isa Float64 ? x :
            x isa Int ? Float64(x) :
            x isa AbstractString && tryparse(Float64, x) !== nothing ? tryparse(Float64, x) :
            missing for x in raw_vec]

        valid_vals = collect(skipmissing(vec))
        if length(valid_vals) < 2
            # Try nearby depths (wider tolerance for deeper targets)
            tol = max(5.0, 0.15 * d)
            neighbor_df = df[abs.(df.depth .- d) .<= tol, :]
            sub = neighbor_df
            raw_vec = sub[!, source_col]
            vec = Union{Missing, Float64}[x isa Missing ? missing :
                x isa Float64 ? x :
                x isa Int ? Float64(x) :
                x isa AbstractString && tryparse(Float64, x) !== nothing ? tryparse(Float64, x) :
                missing for x in raw_vec]
            valid_vals = collect(skipmissing(vec))
            if length(valid_vals) < 2
                continue
            end
        end

        # Use day indices from first observation instead of day-of-year
        day_idx = dayindex.(sub.date)
        gd = groupby(DataFrame(day=day_idx, val=vec), :day)
        day_vals, val_means = Int[], Float64[]
#=
        doy = dayofyear.(sub.date)
        gd = groupby(DataFrame(doy=doy, val=vec), :doy)
        doy_vals, val_means = Int[], Float64[]
=#        
        for g in gd
            v = skipmissing(g.val)
            if !isempty(v)
                push!(day_vals, first(g.day))
                ##push!(doy_vals, first(g.doy))
                push!(val_means, mean(v))
            end
        end

        if length(day_vals) < 2
        ##if length(doy_vals) < 2
            continue
        end

        idx = sortperm(day_vals)
        xgrid = Float64.(day_vals[idx])
        ##idx = sortperm(doy_vals)
        ##xgrid = Float64.(doy_vals[idx])
        ygrid = val_means[idx]

        yhat = if interp_method == "spline"
            k = min(3, length(xgrid)-1)
            Spline1D(xgrid, ygrid; k=k, bc="extrapolate").(days)
        elseif interp_method == "linear"
            LinearInterpolation(xgrid, ygrid, extrapolation_bc=Line()).(days)
        elseif interp_method == "fill"
            fill(last(ygrid), length(days))
        else
            error("Invalid interp_method = $interp_method")
        end

        yconv = convert_fn.(yhat)
##        append!(interpolated, DataFrame(depth=fill(d, length(days)), day=collect(days), value=yconv))
        append!(interpolated, DataFrame(depth=fill(d, length(days)),
                                        day=collect(days), value=yconv))
    end

    return interpolated
end

# ---------------------------
# Run interpolation
# ---------------------------
interpolated_results = Dict{Symbol, DataFrame}()


println("will interpolate...")
#readline()

for param in parameters
    println("🔄 Interpolating $param with $interp_method ...")
    #readline()
    src_col   = param_source_column[param]
    transform = param_transforms[param]
    interpolated_results[param] = interpolate_param(sheet, src_col, transform, ndays, dayindex)
##    interpolated_results[param] = interpolate_param(sheet, src_col, transform)
end

# ---------------------------
# Compute 365-day averaged arrays
# ---------------------------
annual_results = Dict{Symbol, DataFrame}()

for param in parameters
    df = interpolated_results[param]

    # Convert to wide format: rows = day, columns = depth
    df_wide = unstack(df, :day, :depth, :value)
    sort!(df_wide, :day)

    # Depths as column names
    depth_labels = names(df_wide)[2:end]
    depths = parse.(Float64, string.(depth_labels))

    # Extract numeric matrix: rows = days, cols = depths
    Z = Matrix(df_wide[:, Not(:day)])

    ndays_total = size(Z, 1)
    newZ = Array{Float64}(undef, 365, size(Z, 2))

    # For each "day of year" (1–365), average values every 365 days starting from 366
    for i in 1:365
        inds = collect(i+365:365:ndays_total)   # 366+i, 731+i, ...
        if !isempty(inds)
            newZ[i, :] = vec(mean(Z[inds, :], dims=1))
        else
            newZ[i, :] .= NaN
        end
    end

    # Store as DataFrame
    annual_df = DataFrame(DayOfYear = 1:365)
    for (j, d) in enumerate(depths)
        annual_df[!, "depth_$(d)"] = newZ[:, j]
    end

    annual_results[param] = annual_df
end

# ---------------------------
# Save annual results to Excel
# ---------------------------
#=
annual_output_file = joinpath(folder, "annual_means_365.xlsx")
XLSX.openxlsx(annual_output_file, mode="w") do xf
    for param in parameters
        df = annual_results[param]
        ws = XLSX.addsheet!(xf, uppercase(string(param)))
        XLSX.writetable!(ws, Tables.columntable(df); write_columnnames=true)
    end
end
println("💾 Saved annual means to Excel: $annual_output_file")
=#
# ---------------------------
# Copy depth 1.5 values to depth 0.5
# ---------------------------
for param in parameters
    df = annual_results[param]
    if hasproperty(df, "depth_1.5") && hasproperty(df, "depth_0.5")
        df[!, "depth_0.5"] .= df[!, "depth_1.5"]
    end
end

# ---------------------------
# Save annual results to CSV
# ---------------------------

for param in parameters
    df = annual_results[param]
    csv_file = joinpath(folder, "$(param_names[param]).csv")
    CSV.write(csv_file, df; delim=';')
    println("💾 Saved CSV: $csv_file")
end

# ---------------------------
# Plot annual results (day of year vs depth)
# ---------------------------
for param in parameters
    df = annual_results[param]
    days = df.DayOfYear
    depth_cols = [n for n in names(df) if startswith(n, "depth_")]
    depths = sort(parse.(Float64, replace.(depth_cols, "depth_" => "")))
    sorted_cols = ["depth_$(d)" for d in depths]
    Z = Matrix(df[:, sorted_cols])'  # depths × days

    local fig = CairoMakie.Figure(size=(900, 500))
    ax = CairoMakie.Axis(fig[1, 1];
        xlabel="Day of Year", ylabel="Depth (m)",
        title="$(param_names[param]) ($(param))",
        yreversed=true)
    hm = CairoMakie.heatmap!(ax, days, depths, Z'; colormap=CairoMakie.Reverse(:cherry))
    CairoMakie.Colorbar(fig[1, 2], hm; label=string(param))
    CairoMakie.save(joinpath(folder, "$(param_names[param])_annual.png"), fig)
    display(fig)
    println("📊 Saved plot: $(param_names[param])_observed.png")
end


#=
# ---------------------------
# Excel output
# ---------------------------
output_file = joinpath(folder, "interpolated_1.xlsx")

XLSX.openxlsx(output_file, mode="w") do xf
    for param in parameters
        df = interpolated_results[param]
        df_out = rename(df, :value => param)
        ws = XLSX.addsheet!(xf, uppercase(string(param)))
        XLSX.writetable!(ws, Tables.columntable(df_out); write_columnnames=true)
    end
end

println("💾 Saved to Excel: $output_file")

# ---------------------------
# Save to CSV 
# ---------------------------

for param in parameters
    df = interpolated_results[param]
    df_out = rename(df, :value => param)

    # CSV export
    csv_path = joinpath(folder, "$(param)_interpolated.csv")
    CSV.write(csv_path, df_out)
    println("💾 Saved CSV: $csv_path")
 
end
=#