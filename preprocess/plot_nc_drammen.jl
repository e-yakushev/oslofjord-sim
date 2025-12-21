using Oceananigans
using JLD2
using NCDatasets
using NetCDF
using Printf
using Oceananigans.Units
using Oceananigans.Utils: maybe_int
using CairoMakie 
#using Colors
include("/home/eya/src/oslofjord-sim/Oxydep.jl")
#: 
      #Auto, Axis, Figure, GridLayout, Colorbar, 
      #rowgap!, colgap!,GridLayout, Relative, scatter!, lines!,  
      #Observable, Reverse, record, heatmap!, contour!, @lift
#using FjordsSim: 
#    plot_1d_phys,
#    extract_z_faces,
#    record_vertical_tracer,
#    record_surface_speed,
#    plot_ztime,   
#    record_bottom_tracer,
#    BGCModel,
#    oxygen_saturation  
#include("/home/eya/src/FjordsSim.jl/src/BGCModels/BGCModels.jl")
#using .BGCModels
#include("/home/eya/src/FjordsSim.jl/src/BGCModels/boundary_conditions.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
function prettydays(time_seconds)
    days = round(Int, time_seconds / (24 * 3600))
    return "day $days"
end
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#map_axis_kwargs = (xlabel = "Grid points, East", ylabel = "Grid points, North")
#transect_axis_kwargs = (xlabel = "Grid points, East", ylabel = "z (m)")
#framerate = 12

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# replace zeros with NaN in 4D array slice
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function replace_zeros_with_NaN!(A, depth_index, day_index)
    slice = Float64.(view(A, :, :, depth_index, day_index))
    @. slice = ifelse(slice == 0, NaN, slice)
    return slice
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# vertical distribution changes at a point (i,j)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helper function to safely get the "interior" data
get_interior(A, inds...) = 
    hasmethod(interior, Tuple{typeof(A), Vararg{Any}}) ? interior(A, inds...) : A[inds...]
# ~~~~~~~~~~~~~~~~  Seasonal changes at a point (i,j) ~~~~~~~~~~~~~~~~
function plot_ztime(NUT, O₂, O₂_sat, PHY, HET, T, DOM, POM, S, i, j, times, z, folder)
#using CairoMakie  # ensures Makie symbols are available here
    fig = Figure(size = (1500, 1000), fontsize = 20)

    axis_kwargs = (
        xlabel = "Time (days)",
        ylabel = "z (m)",
        xticks = (0:50:times[end]),
        xtickformat = "{:.0f}",
    )

    axNUT = Axis(fig[1, 1]; title = "NUT [μM N]", axis_kwargs...)
    hmNUT = heatmap!(times / days, z, get_interior(NUT, i, j, :, :)', colormap = Reverse(:cherry))
    Colorbar(fig[1, 2], hmNUT)
    
    axOXY = Axis(fig[1, 3]; title = "O₂ [μM]", axis_kwargs...)
#    hmOXY = heatmap!(times / days, z, get_interior(O₂, i, j, :, :)', colormap = :turbo)
#    Colorbar(fig[1, 4], hmOXY)
    O₂_slice = get_interior(O₂, i, j, :, :)'   # transpose so z is vertical
    hmOXY = heatmap!(times / days, z, O₂_slice, colormap = :turbo)
# --- Add isoline O₂ = 67 uM (i.e. 1.5 ml/l as hypoxia threshold) ---
    contour!(times / days, z, O₂_slice; levels=[67], color=:white, linewidth=4, linestyle = :dot)
# --- Add manual label ---
#    text!(times[end] / (2 * days), 20;  # (x, z) position of the label
#      text = "67", color = :white, align = (:center, :center), fontsize = 18, font = "sans")
    Colorbar(fig[1, 4], hmOXY)

    axOXY_rel = Axis(fig[1, 5]; title = "O₂ saturation [%]", axis_kwargs...)

    O₂_sat_slice = get_interior(O₂_sat, i, j, :, :)'   # transpose so z is vertical
    hmOXY_rel = heatmap!(times / days, z, O₂_sat_slice, colormap = :gist_stern) 
    # --- Add isoline O₂_sat = 100 ---
    contour!(times / days, z, O₂_sat_slice; levels=[100], color=:white, 
                                            linewidth=4, linestyle = :dot)
    Colorbar(fig[1, 6], hmOXY_rel)

    axPHY = Axis(fig[2, 1]; title = "PHY [μM N]", axis_kwargs...)
    hmPHY = heatmap!(times / days, z, get_interior(PHY, i, j, :, :)', colormap = Reverse(:cubehelix))
    Colorbar(fig[2, 2], hmPHY)

    axHET = Axis(fig[2, 3]; title = "HET [μM N]", axis_kwargs...)
    hmHET = heatmap!(times / days, z, get_interior(HET, i, j, :, :)', colormap = Reverse(:afmhot))
    Colorbar(fig[2, 4], hmHET)

    axT = Axis(fig[2, 5]; title = "T [°C]", axis_kwargs...)
    hmT = heatmap!(times / days, z, get_interior(T, i, j, :, :)', colormap = Reverse(:RdYlBu))
    Colorbar(fig[2, 6], hmT)

    axDOM = Axis(fig[3, 1]; title = "DOM [μM N]", axis_kwargs...)
    hmDOM = heatmap!(times / days, z, get_interior(DOM, i, j, :, :)', colormap = Reverse(:CMRmap))
    Colorbar(fig[3, 2], hmDOM)

    axPOM = Axis(fig[3, 3]; title = "POM [μM N]", axis_kwargs...)
    hmPOM = heatmap!(times / days, z, get_interior(POM, i, j, :, :)', colormap = Reverse(:greenbrownterrain))
    Colorbar(fig[3, 4], hmPOM)

    axS = Axis(fig[3, 5]; title = "S [psu]", axis_kwargs...)
    hmS = heatmap!(times / days, z, get_interior(S, i, j, :, :)', colormap = :viridis)
    Colorbar(fig[3, 6], hmS)

    save(joinpath(folder, "ztime_$(i)_$(j).png"), fig)
    @info "Saved ztime_$(i)_$(j) plot in $folder"
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot a transect along the deeppest line of the fjord
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_param_transect(Par_transect_slice, title_str, depth, transect, folder; 
                             colormap=:turbo, day_index=1, whiteline=1.0)
    
    # --- Ensure depth is a 1-D vector ---
    depth = vec(depth)
    nmax = length(transect)
    
    # --- Compute cumulative distance along transect ---
    dist = zeros(Float64, nmax)
    for n in 2:nmax
        (_, i1, j1, _) = transect[n-1]
        (_, i2, j2, _) = transect[n]
        dist[n] = dist[n-1] + hypot(i2 - i1, j2 - j1) * 0.2
    end

    # --- Create figure ---
    fig = Figure(size = (900, 450), fontsize = 24)
    ax = Axis(fig[1, 1];
        xlabel = "Distance along transect (km)",
        ylabel = "Depth (m)",
        title  = "$title_str transect for day: $day_index",
        yreversed = false,
    )
    # ============================================
    # BASIC METHOD: First draw the gray background, then the data
    # ============================================
    
    # Create a copy of the data with NaN replacement
    data_plot = copy(Par_transect_slice)
    nan_mask = isnan.(data_plot)
    
    # Find the min/max of valid values
    if any(.!nan_mask)
        valid_data = data_plot[.!nan_mask]
        vmin, vmax = extrema(valid_data)
        
        # Replace NaN with a value BELOW the minimum
        offset = 0.1 * abs(vmax - vmin)
        data_plot[nan_mask] .= vmin - offset
          
        # Drawing a heatmap
        hm = heatmap!(ax, dist, depth, data_plot;
            colorrange = (vmin, vmax),
            # colormap = Reverse(:greenbrownterrain), #custom_colors,  # CairoMakieustom map
            #lowclip = custom_colors[1],  # Gray for values ​​below vmin
            #highclip = custom_colors[end],  # Last color for values ​​above vmax
            interpolate = false,
            nan_color = :silver,  
        )
    else
        # All NaN values ​​- just draw gray
        data_plot .= 0.0
        hm = heatmap!(ax, dist, depth, data_plot;
            colorrange = (0.0, 1.0),
            colormap = [:silver],
            interpolate = false,
        )
    end
    
    # --- Contour line ---
    if whiteline != 0.0
        contour_data = copy(Par_transect_slice)
        if any(nan_mask)
            # Replace NaN with a large value
            contour_data[nan_mask] .= whiteline + 1000.0
        end
        
        contour!(ax, dist, depth, contour_data; 
                 levels=[whiteline], color=:white, linewidth=4, linestyle=:dot)
    end
    
    # --- Colorbar ---
    if any(.!nan_mask)
        Colorbar(fig[1, 2], hm, label="O₂ [μM]", width=25)
    end
    
    # --- Save ---
    title_short = title_str[1:min(2, length(title_str))]
    filename = joinpath(folder, "transect_$(title_short)_day_$(day_index).png")
    save(filename, fig, px_per_unit=2)
    
    @info "Saved plot to $filename"
    println("NaN pixels $(sum(nan_mask)) of $(length(nan_mask))")
    
    return fig
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot a map of bottom depth indices or physical depth (m)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_bottom_depth_map!(fig, pos, bottom_z::AbstractMatrix{<:Integer}, z_vals::AbstractVector;
                                title_str="Bottom depth [m]", use_abs=true, colormap=:viridis, whiteline=0.0)
    # build 2D array of physical depths from index map (preserves shape)
    depth_map = [ z_vals[ bottom_z[i,j] ] for i in 1:size(bottom_z,1), j in 1:size(bottom_z,2) ]
    # convert to absolute positive depth if requested (z often negative)
    if use_abs
        depth_map = abs.(depth_map)
    end
    # mark invalid/zero indices as NaN
    @. depth_map = ifelse(bottom_z == 0, NaN, depth_map)
    # determine sensible colorrange from finite values
    finite_vals = depth_map[isfinite.(depth_map)]
    colorrange = isempty(finite_vals) ? (0.0, 1.0) : (minimum(finite_vals), maximum(finite_vals))
    plot_tracer_subplot!(fig, pos, depth_map, title_str; colorrange=colorrange, colormap=colormap, whiteline=whiteline)
end
# ====================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN CODE STARTS HERE: open file and extract data to "ds"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
folder = joinpath(homedir(), "FjordSim_results", "oslofjord") #joinpath(homedir(), "FjordsSim_results", "oslofjord")
filename = joinpath(folder, "snapshots_ocean")
ds = NCDataset("$filename.nc", "r")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
println("1D Variables in $filename.nc:")
    for (varname, var) in ds
        if ndims(var) == 1
            println("$varname, Dimensions: ", dimnames(var), ", Sizes: ", size(var))
        end
    end
println("2D Variables in $filename.nc:")
    for (varname, var) in ds
        if ndims(var) == 2
            println("$varname, Dimensions: ", dimnames(var),", Sizes: ", size(var))
        end
    end
println("3D Variables in $filename.nc:")
    for (varname, var) in ds
        if ndims(var) == 3
            println("$varname, Dimensions: ", dimnames(var),", Sizes: ", size(var))
        end
    end
println("4D Variables in $filename.nc:")
    for (varname, var) in ds
        if ndims(var) == 4
            println("$varname, Dimensions: ", dimnames(var), ", Sizes: ", size(var))
        end
    end
println("Groups in file:")        
for (name, grp) in ds.group
    println(" - ", name)
end
grid_group = ds.group["grid_attributes"]
println("Variables in grid_attributes:")
for (name, var) in grid_group
    println(" - $name  dims=$(dimnames(var)) size=$(size(var))")
end
# Get grid dimensions and properties
Nx = Int(grid_group.attrib["Nx"])
Ny = Int(grid_group.attrib["Ny"])
Nz = Int(grid_group.attrib["Nz"])
 
println("Grid dimensions: Nx=$Nx,", typeof(Nx),", Ny=$Ny,", typeof(Ny),", Nz=$Nz,", typeof(Nz))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract from NetCDF dataset 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
times = ds["time"][:]         # time in seconds, float
depth = ds["z_aac"][:]        # depth in m, starting from the deepest, negative, float
println("Time stats (in seconds) — min: ", minimum(times), ", max: ", maximum(times))
# --- Extract  T from NetCDF dataset --- 
T = ds["T"][:,:,:,:]               
println("T stats — min: ", minimum(T), ", max: ", maximum(T))
# --- Extract  S from NetCDF dataset --- 
S = ds["S"][:,:,:,:]                
println("S stats — min: ", minimum(S), ", max: ", maximum(S))
# --- Extract  P from NetCDF dataset --- 
P = ds["P"][:,:,:,:]               
println("P stats — min: ", minimum(P),  ", max: ", maximum(P))
# --- Extract  NUT from NetCDF dataset --- 
HET = ds["HET"][:,:,:,:]                 
println("HET stats — min: ", minimum(HET),  ", max: ", maximum(HET))
# --- Extract  NUT from NetCDF dataset --- 
NUT = ds["NUT"][:,:,:,:]                 
println("NUT stats — min: ", minimum(NUT),  ", max: ", maximum(NUT))
# --- Extract  DOM from NetCDF dataset --- 
POM = ds["POM"][:,:,:,:]                 
println("POM stats — min: ", minimum(POM),  ", max: ", maximum(POM))
# --- Extract  DOM from NetCDF dataset --- 
DOM = ds["DOM"][:,:,:,:]                 
println("DOM stats — min: ", minimum(DOM),  ", max: ", maximum(DOM))
# --- Extract  O₂ from NetCDF dataset --- 
O₂ = ds["O₂"][:,:,:,:]                 
println("O₂ stats — min: ", minimum(O₂),  ", max: ", maximum(O₂))

@info "BGH arrays loaded"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# subplot function for tracer plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Example for the Oslofjord (approximate boundaries)
#oslo_fjord_lon_min = 10.45   # East longitude
#oslo_fjord_lon_max = 10.8
#oslo_fjord_lat_min = 59.6  # Northern latitude
#oslo_fjord_lat_max = 59.95


# Example for the Drammensfjord (approximate boundaries)
oslo_fjord_lon_min = 10.23   # East longitude
oslo_fjord_lon_max = 10.45
oslo_fjord_lat_min = 59.585  # Northern latitude
oslo_fjord_lat_max = 59.755

# Create realistic coordinates
real_lon = range(oslo_fjord_lon_min, oslo_fjord_lon_max, length=Nx)
real_lat = range(oslo_fjord_lat_min, oslo_fjord_lat_max, length=Ny)

# Using real coordinates in the plotting function
function plot_tracer_subplot!(fig, pos, data, title_str; 
    longitudes=real_lon, latitudes=real_lat,
    colorrange=(0,1), colormap=:viridis, whiteline=1.0)

    ax = Axis(fig[pos...]; title=title_str, 
            width = 180,
            height = 300,
            xlabel="Longitude (°E)", ylabel="Latitude (°N)")
    
    hm = heatmap!(ax, longitudes, latitudes, data; 
                  colorrange=colorrange, colormap=colormap, nan_color=:silver)
    
    if whiteline != 0.0
        contour!(ax, longitudes, latitudes, data; 
                 levels=[whiteline], color=:white, linewidth=4, linestyle=:dot)
    end
    Colorbar(fig[pos[1], pos[2]+1], hm, vertical=true)
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate additional fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# - - - - 
# Pressure is needed in atm; we calculate  ro*g*z  in Pa and convert to atm using 1 atm = 101325 Pa
# - - - - 
Pressure = similar(T[:,:,:,1])  
for k in 1:Nz
    Pressure[:,:,k] .= 1. + 0.0992*(-depth[k])
end
# - - - - 
# Oxygen saturation related
# - - - - 
O₂_sat_val = similar(T)  # Oxygen saturation concentrtaion
O₂_sat = similar(T)       # Oxygen %
ϵ = eps(Float32)            # small positive number needed in division to zero
# - - - - 
# Bottom depths index for the bottom maps
# - - - - 
bottom_z = ones(Int, size(O₂, 1), size(O₂, 2))
for i = 1:size(O₂, 1)
    for j = 1:size(O₂, 2)
        for k = 1:size(O₂, 3)
            if O₂[i, j, k, 1] .!= 0
                bottom_z[i, j] = k
                break
                if k == Nz
                    bottom_z[i, j] = Nz
                end
            end
            if k == Nz
                bottom_z[i, j] = Nz
            end
        end
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute transect array (list of i, j, max_depth_index)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function compute_transect(bottom_z)
    transect = Vector{Tuple{Int, Int, Int, Int}}()  # (num, i, j, max_depth_index)
    num = 1
    drammen = 1
  if drammen == 1
        println("Computing transect for Drammen fjord...")
        for j in 1:66
        max_depth_index = Nz
        for i in 1:17
            if max_depth_index < bottom_z[i+1, j]
                push!(transect, (num, i, j, max_depth_index))
                num += 1
                break
            end
            max_depth_index = bottom_z[i+1, j]
        end
    end
  else
        println("Computing transect for Oslo fjord...")
    # First part: λ_caa < 27 → φ_aca increases
    for j in 1:71
        max_depth_index = 12
        for i in 27:-1:14
            if max_depth_index < bottom_z[i-1, j]
                push!(transect, (num, i, j, max_depth_index))
                num += 1
                break
            end
            max_depth_index = bottom_z[i-1, j]
        end
    end
    # Second part: λ_caa ≥ 27 → φ_aca decreases
    for j in 71:-1:34
        max_depth_index = 12
        for i in 27:44
            if max_depth_index < bottom_z[i+1, j]
                push!(transect, (num, i, j, max_depth_index))
                num += 1
                break
            end
            max_depth_index = bottom_z[i+1, j]
        end
    end
  end
    println("✅ Total transect points found: ", length(transect))
    return transect
end

# Compute transect trajectory
transect = compute_transect(bottom_z)

# Extract (i, j) positions for plotting
i_vals = [t[2] for t in transect]
j_vals = [t[3] for t in transect]
# ──────────────────────────────────────────────
# Plot a map of bottom depth indices (physical depth in meters)
fig_depth_map0 = Figure(size=(1200, 1000))
plot_bottom_depth_map!(fig_depth_map0, (1, 1), bottom_z, depth; 
    title_str="Bottom depth (m)", use_abs=true, colormap=Reverse(:oslo25), whiteline=0.0)
save(joinpath(folder, "Bottom_depth_map.png"), fig_depth_map0)
println("Saved: Bottom_depth_map.png")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute slice for vertical transect of Parameter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function vert_transect_slice(Param, transect, t, Nz)
    nmax = length(transect)
    Param_slice = Array{Float64}(undef, nmax, Nz)

    for (n, (_, i, j, _)) in enumerate(transect)
        @inbounds Param_slice[n, :] = Param[i, j, 1:Nz, t]
    end

    return Param_slice
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot bottom depth map and overlay transect line
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
using CairoMakie
# Create figure and axis
fig_depth_map = Figure(size = (1200, 1000))

# Plot bottom map and capture both the axis and heatmap
ax_depth = Axis(fig_depth_map[1, 1], title = "Bottom depth (m)",
    xlabel = "Longitude (°E)",  
    ylabel = "Latitude (°N)")   

# Convert indices into coordinates for the trajectory
lon_vals = real_lon[i_vals]
lat_vals = real_lat[j_vals]

# Create a depth map in coordinates
depth_map = [ depth[bottom_z[i,j]] for i in 1:size(bottom_z,1), j in 1:size(bottom_z,2) ]
depth_map = abs.(depth_map)  # absolute depth values
@. depth_map = ifelse(bottom_z == 0, NaN, depth_map)

hm = heatmap!(ax_depth, real_lon, real_lat, depth_map; colormap = Reverse(:oslo25))
cb = Colorbar(fig_depth_map[1, 2], hm, label = "Depth (m)")

# Overlay trajectory line on the *axis* using coordinates
CairoMakie.lines!(ax_depth, lon_vals, lat_vals;
    color = :white,
    linewidth = 2.5,
    linestyle = :solid)

# Save
save(joinpath(folder, "Bottom_depth_and_transect_map.png"), fig_depth_map)
println("💾 Saved: Bottom_depth_and_transect_map.png ✅")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#""" Function to calculate oxygen saturation in seawater """
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function oxygen_saturation(T::Float64, S::Float64, P::Float64)::Float64
    T_kelvin = T + 273.15  # Convert temperature to Kelvin

    # Calculate the natural logarithm of oxygen saturation concentration
    # Coefficients from Garcia and Gordon (1992)
    ln_O2_sat =
        -173.4292 +
        249.6339 * (100 / T_kelvin) +
        143.3483 * log(T_kelvin / 100) +
        -21.8492 * T_kelvin / 100 +
        -0.033096 * (T_kelvin / 100)^2 +
        0.014259 * (T_kelvin / 100)^3 +
        S * (-0.035274 + 0.001429 * (T_kelvin / 100) + -0.00007292 * (T_kelvin / 100)^2) +
        0.0000826 * S^2

    # Oxygen saturation concentration in µmol/kg
    O2_sat = exp(ln_O2_sat) * 44.66

    # Pressure correction factor (Weiss, 1970) for pressure in atm
    P_corr = 1.0 + P * (5.6e-6 + 2.0e-11 * P)

    # Adjusted oxygen saturation with pressure correction
    return (O2_sat * P_corr)
end

# - - - - 
# Fill oxygen saturation and percentage arrays
# - - - - 
for i = 1:size(O₂, 1)
    for j = 1:size(O₂, 2)
        for k = 1:size(O₂, 3)
            for it = 1:size(O₂, 4)
                O₂_sat_val[i, j, k, it] = oxygen_saturation(
                    Float64(T[i, j, k, it]),
                    Float64(S[i, j, k, it]),
                    Float64(Pressure[i, j, k])
                )
                denom = O₂_sat_val[i, j, k, it] == 0f0 ? ϵ : O₂_sat_val[i, j, k, it]
                O₂_sat[i, j, k, it] = 100f0 * O₂[i, j, k, it] / denom
            end
        end
    end
end
println("O₂_sat_val stats — min: ", minimum(O₂_sat_val),  ", max: ", maximum(O₂_sat_val))
println("O₂_sat % stats — min: ", minimum(O₂_sat),  ", max: ", maximum(O₂_sat))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  vertical distributions changes in a point
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plot_ztime(NUT, O₂, O₂_sat, P, HET, T, DOM, POM, S, 15, 52, times, depth, folder) # Vestfjorden
#plot_ztime(NUT, O₂, O₂_sat, P, HET, T, DOM, POM, S, 35, 50, times, depth, folder) # Bunnefjorden
plot_ztime(NUT, O₂, O₂_sat, P, HET, T, DOM, POM, S, 16, 22, times, depth, folder) # Bunnefjorden
plot_ztime(NUT, O₂, O₂_sat, P, HET, T, DOM, POM, S, 10, 44, times, depth, folder) # Bunnefjorden
#plot_ztime(NUT, O₂, O₂_sat, P, HET, T, DOM, POM, S, 6, 41, times, depth, folder) # Bunnefjorden

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot parameters MAPs as subplots at given days(plot_dates), depths(depth_indexes) and bottom
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_dates = [36, 72, 108, 144, 180, 226, 262, 298, 334]
bottom_layer = Nz
depth_indexes = [Nz] # [16]  # surface slice index Oslo:[12] 
fig_width =  950         # figure width
fig_height = 1150         # figure height

for plot_day in plot_dates
    day_index = plot_day * round(Int, length(times)/365) 
    for depth_index in depth_indexes
        println("Plotting full map figure for day $plot_day ...") 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot vertical slices at prescribed transect for the day_index
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
O2_slice = vert_transect_slice(O₂, transect, day_index, Nz)
fig1 = plot_param_transect(O2_slice, "O₂ [μM]", depth, transect, folder; 
        colormap=:turbo, day_index=plot_day, whiteline=67.0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot maps of horizontal slices at given depth_index and day_index
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        # Extract horizontal slices
        T_slice = replace_zeros_with_NaN!(T, depth_index, day_index)
        S_slice = replace_zeros_with_NaN!(S, depth_index, day_index)
        O₂_slice = replace_zeros_with_NaN!(O₂, depth_index, day_index)
        NUT_slice = replace_zeros_with_NaN!(NUT, depth_index, day_index)
        P_slice = replace_zeros_with_NaN!(P, depth_index, day_index)
        HET_slice = replace_zeros_with_NaN!(HET, depth_index, day_index)
        DOM_slice = replace_zeros_with_NaN!(DOM, depth_index, day_index)
        POM_slice  = replace_zeros_with_NaN!(POM, depth_index, day_index)
        O₂_sat_slice = replace_zeros_with_NaN!(O₂_sat, depth_index, day_index)
        println("Creating figure for day $plot_day at depth index $depth_index ...")

        # Create 3×3 subplot figure
        fig = Figure(size=(fig_width, fig_height))

        plot_tracer_subplot!(fig, (1, 1), T_slice, "T [°C]";  colorrange=(0, 20), colormap=Reverse(:RdYlBu), whiteline=0.0)
        plot_tracer_subplot!(fig, (1, 3), S_slice, "S [psu]"; colorrange=(15, 35), colormap=:viridis, whiteline=0.0)
        plot_tracer_subplot!(fig, (1, 5), O₂_slice,"O₂ [μM]"; colorrange=(0, 350), colormap=:turbo, whiteline=67.0)

        plot_tracer_subplot!(fig, (2, 1), P_slice,     "PHY [μM N]";   colorrange=(0, 5), colormap=Reverse(:cubehelix), whiteline=0.0)
        plot_tracer_subplot!(fig, (2, 3), HET_slice,   "HET [μM N]";   colorrange=(0, 5), colormap=Reverse(:afmhot), whiteline=0.0)
        plot_tracer_subplot!(fig, (2, 5), O₂_sat_slice,"O₂ [%]";colorrange=(0, 150), colormap=:gist_stern, whiteline=100.0)

        plot_tracer_subplot!(fig, (3, 1), DOM_slice,   "DOM [μM N]"; colorrange=(0, 15), colormap=Reverse(:CMRmap), whiteline=0.0)
        plot_tracer_subplot!(fig, (3, 3), POM_slice,   "POM [μM N]"; colorrange=(0, 5), colormap=Reverse(:greenbrownterrain), whiteline=0.0)
        plot_tracer_subplot!(fig, (3, 5), NUT_slice,   "NUT [μM N]"; colorrange=(0, 40), colormap=Reverse(:cherry), whiteline=0.0)

        save(joinpath(folder, "map_iz_$(depth_index)_day_$(plot_day).png"), fig)
        @info "Saved: map_iz_$(depth_index)_day_$(plot_day).png"
    end
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot bottom maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Create 3×3 subplot figure
    fig_b = Figure(size=(fig_width, fig_height))

    # Preallocate result arrays (same horizontal dimensions as O₂)
    T_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    S_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    O₂_slice_bot = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))

    P_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    HET_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    O₂_sat_slice_bot = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))

    DOM_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    POM_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    NUT_slice_bot = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    # Fill them
    for i in 1:size(O₂, 1)
        for j in 1:size(O₂, 2)
            z = bottom_z[i, j]
            T_slice_bot[i, j]  = Float64(T[i, j, z, day_index])
            S_slice_bot[i, j]  = Float64(S[i, j, z, day_index])
            O₂_slice_bot[i, j] = Float64(O₂[i, j, z, day_index])

            P_slice_bot[i, j]  = Float64(P[i, j, z, day_index])
            HET_slice_bot[i, j] = Float64(HET[i, j, z, day_index])
            O₂_sat_slice_bot[i, j] = Float64(O₂_sat[i, j, z, day_index])
        
            DOM_slice_bot[i, j] = Float64(DOM[i, j, z, day_index])
            POM_slice_bot[i, j] = Float64(POM[i, j, z, day_index])
            NUT_slice_bot[i, j] = Float64(NUT[i, j, z, day_index])
                             
        end
    end
    # Replace zeros with NaN
    @. T_slice_bot  = ifelse(T_slice_bot == 0, NaN, T_slice_bot)
    @. S_slice_bot  = ifelse(S_slice_bot == 0, NaN, S_slice_bot)
    @. O₂_slice_bot = ifelse(O₂_slice_bot == 0, NaN, O₂_slice_bot)
    @. P_slice_bot  = ifelse(P_slice_bot == 0, NaN, P_slice_bot)
    @. HET_slice_bot  = ifelse(HET_slice_bot == 0, NaN, HET_slice_bot)
    @. O₂_sat_slice_bot = ifelse(O₂_sat_slice_bot == 0, NaN, O₂_sat_slice_bot)
    @. DOM_slice_bot  = ifelse(DOM_slice_bot == 0, NaN, DOM_slice_bot)
    @. POM_slice_bot  = ifelse(POM_slice_bot == 0, NaN, POM_slice_bot)
    @. NUT_slice_bot = ifelse(NUT_slice_bot == 0, NaN, NUT_slice_bot)

    # Plot bottom maps
    plot_tracer_subplot!(fig_b, (1, 1), T_slice_bot,  "T [°C]"; colorrange=(0, 20), colormap=Reverse(:RdYlBu), whiteline=0.0)
    plot_tracer_subplot!(fig_b, (1, 3), S_slice_bot, "S [psu]"; colorrange=(15, 35), colormap=:viridis, whiteline=0.0)
    plot_tracer_subplot!(fig_b, (1, 5), O₂_slice_bot,"O₂ [μM]"; colorrange=(0, 350), colormap=:turbo, whiteline=67.0)

    plot_tracer_subplot!(fig_b, (2, 1), P_slice_bot,     "PHY [μM N]";   colorrange=(0, 5), colormap=Reverse(:cubehelix), whiteline=0.0)
    plot_tracer_subplot!(fig_b, (2, 3), HET_slice_bot,   "HET [μM N]";   colorrange=(0, 5), colormap=Reverse(:afmhot), whiteline=0.0)
    plot_tracer_subplot!(fig_b, (2, 5), O₂_sat_slice_bot,"O₂ [%]";colorrange=(0, 150), colormap=:gist_stern, whiteline=100.0)

    plot_tracer_subplot!(fig_b, (3, 1), DOM_slice_bot,   "DOM [μM N]"; colorrange=(0, 15), colormap=Reverse(:CMRmap), whiteline=0.0)
    plot_tracer_subplot!(fig_b, (3, 3), POM_slice_bot,   "POM [μM N]"; colorrange=(0, 5), colormap=Reverse(:greenbrownterrain), whiteline=0.0)
    plot_tracer_subplot!(fig_b, (3, 5), NUT_slice_bot,   "NUT [μM N]"; colorrange=(0, 40), colormap=Reverse(:cherry), whiteline=0.0)
    # Save figure with bottom maps
            save(joinpath(folder, "map_bottom_day_$(plot_day).png"), fig_b)
            @info "Saved: map_bottom_day_$(plot_day).png"        
end 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# subplot function for tracer plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_tracer_subplot!(fig, pos, data, title_str; 
    colorrange=(0,1), colormap=:viridis, whiteline=1.0)
    ax = Axis(fig[pos...]; title=title_str, 
            width  = 200, #Auto(),  # Adaptive width
            height = 300, #Auto()  # Adaptive height
            xlabel="", ylabel="")
    hm = heatmap!(ax, data; colorrange=colorrange, colormap=colormap, nan_color=:silver)
    if whiteline != 0.0
        contour!(ax, data; levels=[whiteline], color=:white, linewidth=4, linestyle = :dot)
    end
    Colorbar(fig[pos[1], pos[2]+1], hm, vertical=true)
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make animated gif of changes at selected depth
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#function record_horizontal_tracer(tracer, times, folder, name, label;
#                                longitudes=real_lon, latitudes=real_lat,
#                                colorrange = (-1, 30), colormap = :magma, iz = 10, framerate = 12)
#    Nt = length(times)
#    iter = Observable(1)
#   Ti = @lift begin
#        if tracer isa AbstractArray
#            Ti = tracer[:, :, iz, $iter]
#        elseif tracer isa FieldTimeSeries
#            Ti = interior(tracer[$iter], :, :, iz)
#        else
#            error("Unsupported tracer type: $(typeof(tracer))")
#        end
#        Ti[Ti .== 0] .= NaN
#        Ti
#    end
#    title = @lift "$label at $(prettydays(times[$iter]))"
#    fig = Figure(size = (400, 550), fontsize = 20)
#    ax = Axis(fig[1, 1]; title = title, xlabel="Longitude (°E)", ylabel="Latitude (°N)")
#    hm = heatmap!(ax, longitudes, latitudes, Ti, colorrange = colorrange, colormap = colormap, nan_color = :silver)
#    cb = Colorbar(fig[0, 1], hm, vertical = false)
#    record(fig, joinpath(folder, "movie_$(name)_iz_$iz.gif"), 1:Nt, framerate = framerate) do i
#        iter[] = i
#    end
#    @info "movie_$(name)_iz_$iz record made"
#end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make animated gif of changes at the bottom
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#function record_bottom_tracer(variable, var_name, Nz, times, folder;
#    longitudes=real_lon, latitudes=real_lat,
#    colorrange = (-1, 350), colormap = :turbo, figsize = (1000, 400), framerate = 12)
    # bottom_z evaluation
#    bottom_z = ones(Int, size(variable, 1), size(variable, 2))
#    for i = 1:size(variable, 1)
#        for j = 1:size(variable, 2)
#            for k = size(variable, 3):-1:1  # Loop backwards to find the latest non-zero
#                if variable[i, j, k, 1] == 0
#                    bottom_z[i, j] = k
#                    if k != Nz
#                        bottom_z[i, j] = k + 1
#                    end
#                    break
#                end
#            end
#        end
#    end
#    iter = Observable(1)
#    f = @lift begin
#        x = [variable[i, j, bottom_z[i, j], $iter] for i = 1:size(variable, 1), j = 1:size(variable, 2)]
#        x[x.==0] .= NaN
#        x
#    end
#    title = @lift "bottom $(var_name), μM at " * prettydays(times[$iter])
#    fig = Figure(size = figsize)
#    ax = Axis(fig[1, 1]; title = title, xlabel="Longitude (°E)", ylabel="Latitude (°N)")
#    hm = heatmap!(ax, longitudes, latitudes, f, colorrange = colorrange, colormap = colormap)
#    cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(var_name), μM")
#    Nt = length(times)
#    record(fig, joinpath(folder, "movie_$(var_name).gif"), 1:Nt, framerate = framerate) do i
#        iter[] = i
#    end
#end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make movies at the bottom        
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#record_bottom_tracer(O₂, "O2_bottom", Nz, times, folder;
#    colorrange = (-1, 350), colormap = :turbo, figsize = (400, 550),)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make movies at the surface (iz = Nz) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
# record_horizontal_tracer(NUT, times, folder, "NUTsurf", "NUT [μM N]",
#                          colorrange = (0, 40),colormap = Reverse(:cherry),iz = Nz, )
 # ~~~~~~~~~~~~~~~~~~~~~~~~~
 #record_horizontal_tracer(O₂, times, folder, "O2surf", "O₂ [μM]",
#     colorrange = (0, 350), colormap = :turbo, iz = Nz, )
 # ~~~~~~~~~~~~~~~~~~~~~~~~~
# record_horizontal_tracer(P, times, folder, "PHYsurf", "PHY [μM N]",
 #    colorrange = (0, 5), colormap = Reverse(:cubehelix), iz = Nz, )     
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# record_horizontal_tracer(HET, times, folder, "HETsurf", "HET [μM N]",
#     colorrange = (0, 5), colormap = Reverse(:afmhot), iz = Nz, )
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# record_horizontal_tracer(POM, times, folder, "POMsurf", "POM [μM N]",
#     colorrange = (0, 5), colormap = Reverse(:greenbrownterrain), iz = Nz, )    
# ~~~~~~~~~~~~~~~~~~~~~~~~~
#   record_horizontal_tracer(DOM, times, folder, "DOMsurf", "DOM [μM N]",
#     colorrange = (0, 20), colormap = Reverse(:CMRmap), iz = Nz, )
#println("⏸ Press Enter to continue...")
#readline()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to convert seconds to days and find indices for 365 days
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function get_365_days_indices(times_seconds, total_days=365)
    # Converting seconds to days
    times_days = times_seconds ./ (24 * 3600)
    
    # Find the maximum time in days
    max_days = maximum(times_days)
    
    if max_days < total_days
        @warn "Available data only for $max_days days, but requested $total_days days"
        total_days = Int(floor(max_days))
    end
    
    # Find indexes for 365 days
    target_times = range(0, total_days, length=total_days)
    indices = Int[]
    
    for target_day in target_times
        # Find the nearest available time step
        idx = argmin(abs.(times_days .- target_day))
        push!(indices, idx)
    end
    
    # Remove duplicates
    unique_indices = unique(indices)
    
    @info "Selected $(length(unique_indices)) frames for $total_days days"
    @info "Time range: day 0 to day $total_days"
    
    return unique_indices
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to plot six animations on one sheet for 365 days
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_six_animations_365_days(tracers, times_seconds, folder, labels, longitudes, latitudes;
                                      colorranges, colormaps, iz=Nz, framerate=1, 
                                      figsize=(1800, 1200), fontsize=24)
    
    # Receive indexes for 365 days
    day_indices = get_365_days_indices(times_seconds, 365)
    Nt = length(day_indices)
    
    # Create an Observable for the current day
    day_iter = Observable(1)
    current_frame = Observable(1)
    
    # Convert time to days for display
    times_days = times_seconds ./ (24 * 3600)
    
    # Create the figure
    fig = Figure(size=figsize, fontsize=fontsize)
    
    # Create arrays to store axes and heatmaps
    axes_vec = []
    heatmaps_vec = []
    
    # Create 6 subplots in 2x3 grid
    for i in 1:6
        row = ((i-1) ÷ 3) + 1
        col = ((i-1) % 3) + 1
        
        # Create main plot axis
        ax = Axis(fig[row, col*2-1], 
                  xlabel=(row == 2 ? "Longitude (°E)" : ""),
                  ylabel=(col == 1 ? "Latitude (°N)" : ""),
                  title=labels[i])
        
        # Create observable for current tracer и текущего дня
        Ti = @lift begin
            tracer_data = tracers[i]
            frame_idx = day_indices[$day_iter]
            # Check the dimensions of the data
            if ndims(tracer_data) == 4
                Ti_val = tracer_data[:, :, iz, frame_idx]
            else
                error("Expected 4D array, got $(ndims(tracer_data))D")
            end
            # Convert to Float64 and replace zeros with NaN
            Ti_val = Float64.(Ti_val)
            Ti_val[Ti_val .== 0] .= NaN
            Ti_val
        end
        
        # Create heatmap
        hm = heatmap!(ax, longitudes, latitudes, Ti, 
                     colorrange=colorranges[i], 
                     colormap=colormaps[i], 
                     nan_color=:white)
        
        # Create colorbar next to the plot
        Colorbar(fig[row, col*2], hm, width=20, label="")
        
        push!(axes_vec, ax)
        push!(heatmaps_vec, hm)
    end
    
    # Add super title with time in days
    super_title = @lift begin
        day_idx = day_indices[$day_iter]
        current_day = round(Int, times_days[day_idx])
        "Day $current_day - Surface"
    end
    Label(fig[0, :], super_title, fontsize=24, font=:bold)
    
    # Add progress information
    #progress_info = @lift "Frame $($day_iter)/$Nt"
    #Label(fig[3, :], progress_info, fontsize=22, color=:gray)
    
    # Adjust layout
    rowgap!(fig.layout, 15)
    colgap!(fig.layout, 10)
    
    output_file = joinpath(folder, "six_animations_365days_iz_$iz.gif")
    
    @info "Starting 365-day animation recording with $Nt frames..."
    @info "Output file: $output_file"
    @info "Frame rate: $framerate fps"
    
    # Record animation for selected days
    record(fig, output_file, 1:Nt, framerate=framerate) do i
        day_iter[] = i
        current_frame[] = day_indices[i]
        
        if i % 30 == 0 || i <= 5 || i >= Nt-5
            day_idx = day_indices[i]
            current_day = round(Int, times_days[day_idx])
            @info "Processing: frame $i/$Nt (day $current_day)"
        end
    end
    
    @info "365-day animation completed: $output_file"
    
    return fig
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare and call the function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Prepare data for the six animations
tracers = [NUT, O₂, P, HET, POM, DOM]

labels = ["NUT [μM N]", "O₂ [μM N]", "PHY [μM N]", "HET [μM N]", "POM [μM N]", "DOM [μM N]"]

colorranges = [
    (0.0, 40.0),    # NUT
    (0.0, 350.0),   # O₂
    (0.0, 5.0),     # PHY
    (0.0, 5.0),     # HET  
    (0.0, 5.0),     # POM
    (0.0, 20.0)     # DOM
]

colormaps = [
    :viridis,        # NUT
    :turbo,          # O₂
    :plasma,         # PHY
    :hot,            # HET
    :rainbow,        # POM
    :jet             # DOM
]

# Check the available data
times_days = times ./ (24 * 3600)
max_days = maximum(times_days)
@info "Available data: $(length(times)) time steps, up to day $(round(max_days, digits=1))"

# Select the appropriate function depending on the data
if max_days >= 365
    @info "Enough data for 365-day animation"
    fig = plot_six_animations_365_days(tracers, times, folder, labels, real_lon, real_lat,
                                      colorranges=colorranges, colormaps=colormaps, 
                                      iz=Nz, framerate=1)
else
    @info "Using all available data ($(round(max_days, digits=1)) days)"
    fig = plot_six_animations_daily(tracers, times, folder, labels, real_lon, real_lat,
                                   colorranges=colorranges, colormaps=colormaps, 
                                   iz=Nz, framerate=1)
end

# Close the dataset
close(ds)

@info "Script completed!"
