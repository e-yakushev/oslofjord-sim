using NCDatasets
using XLSX
using Statistics
using Dates

function export_river_data(file_path::String, output_file::String, selected_year::Int,
                           river_index::Int, s_rho_index::Int, param_name::String = "river_N3_n")
    ds = NCDataset(file_path, "r")
    
    param = ds[param_name]  # dimensions: (river, s_rho, river_time)
    time_var = ds["river_time"]
    water_discharge = ds["river_transport"]
   
    river_times_raw = time_var[:]                 # assumed to be DateTime array
    river_values = param[river_index, s_rho_index, :]
    river_discharge = water_discharge[river_index, :]

# to check names of 2D (or 3D) variables in .nc file
#    println("2D Variables in $file_path:")
#    for (varname, var) in ds
#        if ndims(var) == 2
#            println("Variable: $varname")
#            println("  Dimensions: ", dimnames(var))
#            println("  Sizes: ", size(var))
#        end
#    end

    # Pre-allocate filtered data
    selected_days = Int[]
    selected_times = typeof(river_times_raw[1])[]
    selected_values = Float64[]
    selected_discharge = Float64[]

    for (i, t) in enumerate(river_times_raw)
        d = Date(t)  # strip time part for year/day filtering
        if year(d) == selected_year
            push!(selected_days, dayofyear(d))
            push!(selected_times, t)
            push!(selected_values, river_values[i])
            push!(selected_discharge, river_discharge[i])
        end
    end

     # Output diagnostics
    println("📏 Dimensions of 'array': ", dimnames(param))
    println("📏 Sizes of 'array': ", size(param))
    println("📅 Selected year: $selected_year → ", length(selected_values), " entries")
    println("📏 ",param_name," values — min: ", minimum(selected_values), 
        ", mean: ", mean(selected_values), 
        ", max: ", maximum(selected_values))
  #  println("⏸ Press Enter to continue...")
  #  readline()

    close(ds)
#=
    # Write to Excel
    XLSX.openxlsx(output_file * ".xlsx", mode="w") do xf
        sheet = xf[1]
        sheet["A1"] = "DayOfYear"
        sheet["B1"] = "flux,m3/s"
        sheet["C1"] = "con,mmol/m3" 
        sheet["D1"] = "dis,mmol/s" 
        for i in 1:length(selected_values)
            sheet["A$(i+1)"] = selected_days[i]
            sheet["B$(i+1)"] = selected_discharge[i]
            sheet["C$(i+1)"] = selected_values[i]
            sheet["D$(i+1)"] = (selected_values[i] * selected_discharge[i])
        end
    end
=#
    open(output_file * ".txt", "w") do io
        # Write header
        println(io, "DayOfYear\tflux_m3/s\tc_mmol/m3\tdis_mmol/s")
        
        # Write rows
        for i in 1:length(selected_values)
            println(io, "$(selected_days[i])\t$(selected_discharge[i])\t$(selected_values[i])\t$(selected_values[i] * selected_discharge[i])")
        end
    end


    println("📤 Data for river $selected_river for year $selected_year written to: $output_file")
end

# Example usage
selected_year = 2020
selected_river = 18
selected_depth = 1
param_names = ["river_N3_n","river_DON0","river_PON0"]

folder = joinpath(homedir(), "FjordsSim_results", "oslofjord")
file = joinpath(folder, "of800_rivers_13_22.nc")

for ip in 1:length(param_names)
    output = joinpath(folder, "$(param_names[ip])_from_$(selected_river)_year_$(selected_year)")
    export_river_data(file, output, selected_year, selected_river, selected_depth, param_names[ip])
end


