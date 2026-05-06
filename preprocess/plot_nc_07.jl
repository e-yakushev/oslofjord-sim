using Oceananigans
using NCDatasets
using NetCDF
using Printf
using Oceananigans.Units
using Oceananigans.Utils: maybe_int
using CairoMakie
using Statistics

# Put the helper code in a separate file:
include("plot_nc_stuff.jl")
using .plot_nc_stuff
import .plot_nc_stuff: prettydays, replace_zeros_with_NaN!, get_interior,
       plot_ztime, plot_param_transect, plot_bottom_depth_map!,
       compute_transect, compute_pressure_from_depth,
       compute_oxygen_saturation_fields, compute_bottom_index_from_O2,
       vert_transect_slice, make_bottom_depth_and_transect_figure,
       plot_tracer_subplot_map!, bottom_slices_at_day, oxygen_saturation,
       record_horizontal_tracer, record_bottom_tracer,
       get_365_days_indices, plot_six_animations_365_days,
       prepare_six_animations

# ===================== MAIN CODE STARTS HERE =====================
base_dir = dirname(@__DIR__)
folder = joinpath(base_dir, "data", "output", "drammensfjord")
filename = joinpath(folder, "snapshots_ocean")
#filename = joinpath(folder, "snapshots_new10")
#filename = joinpath(folder, "snapshots_ocean")
ds       = NCDataset("$filename.nc", "r")


###################################################
# flags to control which plots to make:
BGC_             = true     # to skip BGC fields and related plots
Ci_              = false     # to skip Ci fields and related plots
maps_animations_ = true    # to prepare individual tracer animations 
six_animations_  = false    # to prepare 6-maps-tracer animations 
# general settings for plots:
plot_dates     = [40, 100, 160, 220, 280, 340] # for maps and transects 
#plot_dates     = [36, 72, 108, 144, 180, 226, 262, 298, 334]

###################################################
# ----- Explore the contents of the NetCDF file ---

println("1D Variables in $filename.nc:")
for (varname, var) in ds
    if ndims(var) == 1
        println("$varname, Dimensions: ", dimnames(var), ", Sizes: ", size(var))
    end
end

println("2D Variables in $filename.nc:")
for (varname, var) in ds
    if ndims(var) == 2
        println("$varname, Dimensions: ", dimnames(var), ", Sizes: ", size(var))
    end
end

println("3D Variables in $filename.nc:")
for (varname, var) in ds
    if ndims(var) == 3
        println("$varname, Dimensions: ", dimnames(var), ", Sizes: ", size(var))
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

Nx = size(ds["λ_caa"], 1)
Ny = size(ds["φ_aca"], 1)
Nz = size(ds["z_aac"], 1)
println("Grid dimensions: Nx=$Nx, ", typeof(Nx), ", Ny=$Ny, ", typeof(Ny), ", Nz=$Nz, ", typeof(Nz))

# ------------------ Extract fields from NetCDF ------------------

times = ds["time"][:]
depth = ds["z_aac"][:]
println("Time stats (in seconds) — min: ", minimum(times), ", max: ", maximum(times))

T   = ds["T"][:, :, :, :]
println("T stats — min: ", minimum(T), ", max: ", maximum(T))

S   = ds["S"][:, :, :, :]
println("S stats — min: ", minimum(S), ", max: ", maximum(S))

e = ds["e"][:, :, :, :]
println("e stats — min: ", minimum(e), ", max: ", maximum(e))

u = ds["u"][:, :, :, :]
println("u stats — min: ", minimum(u), ", max: ", maximum(u))

v = ds["v"][:, :, :, :]
println("v stats — min: ", minimum(v), ", max: ", maximum(v))

if BGC_
    P   = ds["P"][:, :, :, :]
    println("P stats — min: ", minimum(P), ", max: ", maximum(P))

    HET = ds["HET"][:, :, :, :]
    println("HET stats — min: ", minimum(HET), ", max: ", maximum(HET))

    NUT = ds["NUT"][:, :, :, :]
    println("NUT stats — min: ", minimum(NUT), ", max: ", maximum(NUT))

    POM = ds["POM"][:, :, :, :]
    println("POM stats — min: ", minimum(POM), ", max: ", maximum(POM))

    DOM = ds["DOM"][:, :, :, :]
    println("DOM stats — min: ", minimum(DOM), ", max: ", maximum(DOM))

    O₂  = ds["O₂"][:, :, :, :]
    println("O₂ stats — min: ", minimum(O₂), ", max: ", maximum(O₂))
end

if Ci_
    Ci_free = ds["Ci_free"][:, :, :, :]
    println("Ci_free stats — min: ", minimum(Ci_free), ", max: ", maximum(Ci_free))

    Ci_DOM  = ds["Ci_DOM"][:, :, :, :]
    println("Ci_DOM stats — min: ", minimum(Ci_DOM), ", max: ", maximum(Ci_DOM))

    Ci_PHY  = ds["Ci_PHY"][:, :, :, :]
    println("Ci_PHY stats — min: ", minimum(Ci_PHY), ", max: ", maximum(Ci_PHY))

    Ci_HET  = ds["Ci_HET"][:, :, :, :]
    println("Ci_HET stats — min: ", minimum(Ci_HET), ", max: ", maximum(Ci_HET))

    Ci_POM  = ds["Ci_POM"][:, :, :, :]
    println("Ci_POM stats — min: ", minimum(Ci_POM), ", max: ", maximum(Ci_POM))
end

@info "all arrays extracted from NetCDF file. Starting to make plots ..."

##########################################
# ------------------ Coordinates for plotting ------------------

oslo_fjord_lon_min = 10.23
oslo_fjord_lon_max = 10.45
oslo_fjord_lat_min = 59.585
oslo_fjord_lat_max = 59.755

# Use actual data dimensions (may differ from grid coordinate variables due to staggering)
data_Nx = size(T, 1)
data_Ny = size(T, 2)
real_lon = range(oslo_fjord_lon_min, oslo_fjord_lon_max, length=data_Nx)
real_lat = range(oslo_fjord_lat_min, oslo_fjord_lat_max, length=data_Ny)

##########################################
# -----------Compute additional variables ------------------
# oxygen saturation and bottom depth from O₂
if BGC_
    Pressure = compute_pressure_from_depth(T, depth, Nz)
    O₂_sat_val, O₂_sat = compute_oxygen_saturation_fields(T, S, Pressure, O₂)
    println("O₂_sat_val stats — min: ", minimum(O₂_sat_val), ", max: ", maximum(O₂_sat_val))
    println("O₂_sat % stats — min: ", minimum(O₂_sat), ", max: ", maximum(O₂_sat))

    bottom_z = compute_bottom_index_from_O2(O₂, Nz)
end

# compute coordinated of the transect line (i, j) points along the 
# deepest bottom, for each column in the horizontal plane
transect = compute_transect(bottom_z, Nz)
#transect = compute_transect(bottom_z)
i_vals   = [t[2] for t in transect]
j_vals   = [t[3] for t in transect]

##########################################
# Bottom depth map
fig_depth_map0 = Figure(size=(1200, 1000))
plot_bottom_depth_map!(fig_depth_map0, (1, 1), bottom_z, depth;
                       title_str="Bottom depth (m)", use_abs=true,
                       colormap=Reverse(:oslo25), whiteline=0.0)
save(joinpath(folder, "Bottom_depth_map.png"), fig_depth_map0)
println("Saved: Bottom_depth_map.png")

##########################################
# -- Bottom map with transect line -------
make_bottom_depth_and_transect_figure(folder, bottom_z, depth, real_lon, real_lat, i_vals, j_vals)

##########################################
# - Vertical-time plots of BGC variables -
plot_ztime(NUT, O₂, O₂_sat, P, HET, T, DOM, POM, S, 16, 22, times, depth, folder)
plot_ztime(NUT, O₂, O₂_sat, P, HET, T, DOM, POM, S, 10, 44, times, depth, folder)


##########################################
# AVERAGED VALUES OVER THE VOLUME
#    days = ds["times"][:] ./ 4.0   # seconds → days
   nt = length(times)
   dt_days = 6 / 24     # 0.25  
   times_days = (0:nt-1) .* dt_days
#    day_index = plot_day * round(Int, length(times) / 365)   

    # Allocate arrays
    Int_S   = zeros(nt)
    Int_T   = zeros(nt)
    Int_e   = zeros(nt)
    if BGC_    
        Int_NUT = zeros(nt)
        Int_P   = zeros(nt)
        Int_HET = zeros(nt)
        Int_POM = zeros(nt)
        Int_DOM = zeros(nt)
        Int_O2  = zeros(nt)
    end
    # Helper function: NaN-safe mean
    nanmean(A) = mean(filter(!isnan, vec(A)))    
    for t in 1:nt
        Int_S[t]  = mean(ds["S"][:, :, :, t])
        Int_T[t]  = mean(ds["T"][:, :, :, t])
        Int_e[t]  = mean(ds["e"][:, :, :, t])
        if BGC_
            Int_NUT[t] = mean(ds["NUT"][:, :, :, t])
            Int_P[t]   = mean(ds["P"][:, :, :, t])
            Int_HET[t] = mean(ds["HET"][:, :, :, t])
            Int_POM[t] = mean(ds["POM"][:, :, :, t])
            Int_DOM[t] = mean(ds["DOM"][:, :, :, t])
            Int_O2[t]  = mean(ds["O₂"][:, :, :, t])
        end
    end

# ---------------- Plot averaged values over time ---------------- 
fig = Figure(size = (800, 800))

ax1 = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "T")    
ax2 = Axis(fig[1, 2], xlabel = "Time (days)", ylabel = "S")
ax3 = Axis(fig[1, 3], xlabel = "Time (days)", ylabel = "e")
lines!(ax1, times_days, Int_T, linewidth = 2, color = :cyan)
lines!(ax2, times_days, Int_S, linewidth = 2, color = :magenta)
lines!(ax3, times_days, Int_e, linewidth = 2, color = :blue)

if BGC_
    ax4 = Axis(fig[2, 1], xlabel = "Time (days)", ylabel = "P")    
    ax5 = Axis(fig[2, 2], xlabel = "Time (days)", ylabel = "HET")
    ax6 = Axis(fig[2, 3], xlabel = "Time (days)", ylabel = "POM")
    ax7 = Axis(fig[3, 1], xlabel = "Time (days)", ylabel = "DOM")
    ax8 = Axis(fig[3, 2], xlabel = "Time (days)", ylabel = "NUT")
    ax9 = Axis(fig[3, 3], xlabel = "Time (days)", ylabel = "N-total")
lines!(ax4, times_days, Int_P,   linewidth = 2, color = :green)
lines!(ax5, times_days, Int_HET, linewidth = 2, color = :orange)
lines!(ax6, times_days, Int_POM, linewidth = 2, color = :blue)
lines!(ax7, times_days, Int_DOM, linewidth = 2, color = :purple)
lines!(ax8, times_days, Int_NUT, linewidth = 2, color = :red)
lines!(ax9, times_days, Int_NUT .+ Int_P .+ Int_HET .+ Int_POM .+ Int_DOM,  linewidth = 2, color = :black)
lines!(ax9, times_days, Int_NUT, linewidth = 2, color = :red)
lines!(ax9, times_days, Int_NUT .+ Int_P,  linewidth = 2, color = :green)
lines!(ax9, times_days, Int_NUT .+ Int_P .+ Int_HET,  linewidth = 2, color = :orange)
lines!(ax9, times_days, Int_NUT .+ Int_P .+ Int_HET .+ Int_POM,  linewidth = 2, color = :blue)
end

    CairoMakie.ylims!(ax1, nothing, nothing)
    CairoMakie.ylims!(ax2, nothing, nothing)
    CairoMakie.ylims!(ax3, 0, nothing)
if BGC_    
    CairoMakie.ylims!(ax4, 0, nothing)
    CairoMakie.ylims!(ax5, 0, nothing)
    CairoMakie.ylims!(ax6, 0, nothing)
    CairoMakie.ylims!(ax7, 0, nothing)
    CairoMakie.ylims!(ax8, 0, nothing)
    CairoMakie.ylims!(ax9, 0, nothing)
end
filename1 = joinpath(folder, "vol_averaged_timeseries.png")
save(filename1, fig)
println("averaged values figure saved to: ", filename1)

##########################################
# - Maps and transects at selected dates -

bottom_layer   = Nz
depth_indexes  = [Nz]
fig_width      = 950
fig_height     = 1150

for plot_day in plot_dates
    day_index = plot_day * round(Int, length(times) / 365)

    for depth_index in depth_indexes
        println("Plotting full map figure for day $plot_day ...")
        ####################
        # vertical TRANSECTS
        O2_slice = vert_transect_slice(O₂, transect, day_index, Nz)
        _ = plot_param_transect(O2_slice, "O₂ [μM]", depth, transect, folder;
                                colormap=:turbo, day_index=plot_day, whiteline=67.0)
        NUT_slice = vert_transect_slice(NUT, transect, day_index, Nz)
        _ = plot_param_transect(NUT_slice, "NUT [μM]", depth, transect, folder;
                                colormap=:viridis, day_index=plot_day, whiteline=10.0)
        e_slice = vert_transect_slice(e, transect, day_index, Nz)
        _ = plot_param_transect(e_slice, "e [m²/s²]", depth, transect, folder;
                                colormap=:viridis, day_index=plot_day, whiteline=1.0)
        S_slice = vert_transect_slice(S, transect, day_index, Nz)
        _ = plot_param_transect(S_slice, "S [psu]", depth, transect, folder;
                                colormap=:viridis, day_index=plot_day, whiteline=15.0)
        P_slice = vert_transect_slice(P, transect, day_index, Nz)
        _ = plot_param_transect(P_slice, "P [μM]", depth, transect, folder;
                                colormap=:viridis, day_index=plot_day, whiteline=3.0)
        POM_slice = vert_transect_slice(POM, transect, day_index, Nz)
        _ = plot_param_transect(POM_slice, "POM [μM]", depth, transect, folder;
                                colormap=:viridis, day_index=plot_day, whiteline=3.0)
        if Ci_
        Ci_free_slice = vert_transect_slice(Ci_free, transect, day_index, Nz)
        _ = plot_param_transect(Ci_free_slice, "Ci_free [μM]", depth, transect, folder;
                                colormap=:turbo, day_index=plot_day, whiteline=0.01)
        Ci_DOM_slice = vert_transect_slice(Ci_DOM, transect, day_index, Nz)
        _ = plot_param_transect(Ci_DOM_slice, "Ci_DOM [μM]", depth, transect, folder;
                                colormap=:turbo, day_index=plot_day, whiteline=0.01)
        Ci_diss_slice = Ci_free_slice .+ Ci_DOM_slice
        _ = plot_param_transect(Ci_diss_slice, "Ci_diss [μM]", depth, transect, folder;
                                colormap=:turbo, day_index=plot_day, whiteline=0.01)
        Ci_PHY_slice = vert_transect_slice(Ci_PHY, transect, day_index, Nz)
        _ = plot_param_transect(Ci_PHY_slice, "Ci_PHY [μM]", depth, transect, folder;
                                colormap=:turbo, day_index=plot_day, whiteline=0.01)
        Ci_HET_slice = vert_transect_slice(Ci_HET, transect, day_index, Nz)
        _ = plot_param_transect(Ci_HET_slice, "Ci_HET [μM]", depth, transect, folder;
                                colormap=:turbo, day_index=plot_day, whiteline=0.01)
        Ci_POM_slice = vert_transect_slice(Ci_POM, transect, day_index, Nz)
        _ = plot_param_transect(Ci_POM_slice, "Ci_POM [μM]", depth, transect, folder;
                                colormap=:turbo, day_index=plot_day, whiteline=0.01)
        Ci_part_slice = Ci_PHY_slice .+ Ci_HET_slice .+ Ci_POM_slice
        _ = plot_param_transect(Ci_part_slice, "Ci_part [μM]", depth, transect, folder;
                                colormap=:turbo, day_index=plot_day, whiteline=0.01)
        end

        ####################
        # horizontal slices for maps
        T_slice      = replace_zeros_with_NaN!(T, depth_index, day_index)
        S_slice      = replace_zeros_with_NaN!(S, depth_index, day_index)
        O₂_slice     = replace_zeros_with_NaN!(O₂, depth_index, day_index)
        NUT_slice    = replace_zeros_with_NaN!(NUT, depth_index, day_index)
        P_slice      = replace_zeros_with_NaN!(P, depth_index, day_index)
        HET_slice    = replace_zeros_with_NaN!(HET, depth_index, day_index)
        DOM_slice    = replace_zeros_with_NaN!(DOM, depth_index, day_index)
        POM_slice    = replace_zeros_with_NaN!(POM, depth_index, day_index)
        O₂_sat_slice = replace_zeros_with_NaN!(O₂_sat, depth_index, day_index)

        if Ci_
            Ci_free_slice = replace_zeros_with_NaN!(Ci_free, depth_index, day_index)
            Ci_DOM_slice  = replace_zeros_with_NaN!(Ci_DOM, depth_index, day_index)
            Ci_PHY_slice  = replace_zeros_with_NaN!(Ci_PHY, depth_index, day_index)
            Ci_POM_slice  = replace_zeros_with_NaN!(Ci_POM, depth_index, day_index)
            Ci_HET_slice  = replace_zeros_with_NaN!(Ci_HET, depth_index, day_index)
        end

        v_slice = replace_zeros_with_NaN!(v, depth_index, day_index)
        u_slice = replace_zeros_with_NaN!(u, depth_index, day_index)

        ####################
        # HORIZONTAL MAPS
        println("Creating maps for day $plot_day at depth index $depth_index ...")
        fig = Figure(size=(fig_width, fig_height))
        plot_tracer_subplot_map!(fig, (1, 1), T_slice,       "T [°C]",       real_lon, real_lat;
                                 colorrange=(0, 20), colormap=Reverse(:RdYlBu), whiteline=0.0)
        plot_tracer_subplot_map!(fig, (1, 3), S_slice,       "S [psu]",      real_lon, real_lat;
                                 colorrange=(15, 35), colormap=:viridis, whiteline=0.0)
        plot_tracer_subplot_map!(fig, (1, 5), O₂_slice,      "O₂ [μM]",      real_lon, real_lat;
                                 colorrange=(0, 350), colormap=:turbo, whiteline=67.0)

        plot_tracer_subplot_map!(fig, (2, 1), P_slice,       "PHY [μM N]",   real_lon, real_lat;
                                 colorrange=(0, 5), colormap=Reverse(:cubehelix), whiteline=0.0)
        plot_tracer_subplot_map!(fig, (2, 3), HET_slice,     "HET [μM N]",   real_lon, real_lat;
                                 colorrange=(0, 5), colormap=Reverse(:afmhot), whiteline=0.0)
        plot_tracer_subplot_map!(fig, (2, 5), O₂_sat_slice,  "O₂ [%]",       real_lon, real_lat;
                                 colorrange=(0, 150), colormap=:gist_stern, whiteline=100.0)

        plot_tracer_subplot_map!(fig, (3, 1), DOM_slice,     "DOM [μM N]",   real_lon, real_lat;
                                 colorrange=(0, 15), colormap=Reverse(:CMRmap), whiteline=0.0)
        plot_tracer_subplot_map!(fig, (3, 3), POM_slice,     "POM [μM N]",   real_lon, real_lat;
                                 colorrange=(0, 5), colormap=Reverse(:greenbrownterrain), whiteline=0.0)
        plot_tracer_subplot_map!(fig, (3, 5), NUT_slice,     "NUT [μM N]",   real_lon, real_lat;
                                 colorrange=(0, 40), colormap=Reverse(:cherry), whiteline=0.0)

        save(joinpath(folder, "map_iz_$(depth_index)_day_$(plot_day).png"), fig)
        @info "Saved: map_iz_$(depth_index)_day_$(plot_day).png"

        if Ci_
            fig2 = Figure(size=(fig_width, fig_height))

            plot_tracer_subplot_map!(fig2, (1, 1), Ci_free_slice, "Ci_free [μM]",  real_lon, real_lat;
                                     colorrange=(0, 1), colormap=:turbo, whiteline=0.0)
            plot_tracer_subplot_map!(fig2, (1, 3), Ci_DOM_slice,  "Ci_DOM [μM]",   real_lon, real_lat;
                                     colorrange=(0, 0.01), colormap=Reverse(:cherry), whiteline=0.0)
            Ci_part_slice = Ci_PHY_slice .+ Ci_HET_slice .+ Ci_POM_slice
            plot_tracer_subplot_map!(fig2, (1, 5), Ci_part_slice, "Ci_part [μM]",  real_lon, real_lat;
                                     colorrange=(0, 0.01), colormap=Reverse(:turbo), whiteline=0.0)
            plot_tracer_subplot_map!(fig2, (2, 1), Ci_PHY_slice,  "Ci_PHY [μM]",   real_lon, real_lat;
                                     colorrange=(0, 0.01), colormap=Reverse(:cubehelix), whiteline=0.0)
            plot_tracer_subplot_map!(fig2, (2, 3), Ci_HET_slice,  "Ci_HET [μM N]", real_lon, real_lat;
                                     colorrange=(0, 0.01), colormap=Reverse(:afmhot), whiteline=0.0)
            plot_tracer_subplot_map!(fig2, (2, 5), Ci_POM_slice,  "Ci_POM [μM]",   real_lon, real_lat;
                                     colorrange=(0, 0.01), colormap=Reverse(:gist_stern), whiteline=0.0)
            plot_tracer_subplot_map!(fig2, (3, 1), u_slice,       "u [m/s]",       real_lon, real_lat;
                                     colorrange=(-0.2, 0.2), colormap=Reverse(:seismic), whiteline=0.0)
            plot_tracer_subplot_map!(fig2, (3, 3), v_slice,       "v [m/s]",       real_lon, real_lat;
                                     colorrange=(-0.2, 0.2), colormap=Reverse(:seismic), whiteline=0.0)
             Ci_diss_slice = Ci_free_slice .+ Ci_DOM_slice
            plot_tracer_subplot_map!(fig2, (3, 5), Ci_diss_slice, "Ci_diss [μM]",  real_lon, real_lat;
                                     colorrange=(0, 1), colormap=Reverse(:turbo), whiteline=0.0)
            save(joinpath(folder, "map2_iz_$(depth_index)_day_$(plot_day).png"), fig2)
            @info "Saved: map2_iz_$(depth_index)_day_$(plot_day).png"
        end
    end

    # bottom maps for this day
    fig_b = Figure(size=(fig_width, fig_height))

    T_slice_bot, S_slice_bot, O₂_slice_bot,
    P_slice_bot, HET_slice_bot, O₂_sat_slice_bot,
    DOM_slice_bot, POM_slice_bot, NUT_slice_bot =
        bottom_slices_at_day(T, S, O₂, P, HET, O₂_sat, DOM, POM, NUT, bottom_z, day_index)

    plot_tracer_subplot_map!(fig_b, (1, 1), T_slice_bot, "T [°C]",      real_lon, real_lat;
                             colorrange=(0, 20), colormap=Reverse(:RdYlBu), whiteline=0.0)
    plot_tracer_subplot_map!(fig_b, (1, 3), S_slice_bot, "S [psu]",     real_lon, real_lat;
                             colorrange=(15, 35), colormap=:viridis, whiteline=0.0)
    plot_tracer_subplot_map!(fig_b, (1, 5), O₂_slice_bot,"O₂ [μM]",     real_lon, real_lat;
                             colorrange=(0, 350), colormap=:turbo, whiteline=67.0)

    plot_tracer_subplot_map!(fig_b, (2, 1), P_slice_bot, "PHY [μM N]",  real_lon, real_lat;
                             colorrange=(0, 5), colormap=Reverse(:cubehelix), whiteline=0.0)
    plot_tracer_subplot_map!(fig_b, (2, 3), HET_slice_bot,"HET [μM N]", real_lon, real_lat;
                             colorrange=(0, 5), colormap=Reverse(:afmhot), whiteline=0.0)
    plot_tracer_subplot_map!(fig_b, (2, 5), O₂_sat_slice_bot, "O₂ [%]", real_lon, real_lat;
                             colorrange=(0, 150), colormap=:gist_stern, whiteline=100.0)

    plot_tracer_subplot_map!(fig_b, (3, 1), DOM_slice_bot, "DOM [μM N]", real_lon, real_lat;
                             colorrange=(0, 15), colormap=Reverse(:CMRmap), whiteline=0.0)
    plot_tracer_subplot_map!(fig_b, (3, 3), POM_slice_bot, "POM [μM N]", real_lon, real_lat;
                             colorrange=(0, 5), colormap=Reverse(:greenbrownterrain), whiteline=0.0)
    plot_tracer_subplot_map!(fig_b, (3, 5), NUT_slice_bot, "NUT [μM N]", real_lon, real_lat;
                             colorrange=(0, 40), colormap=Reverse(:cherry), whiteline=0.0)

    save(joinpath(folder, "map_bottom_day_$(plot_day).png"), fig_b)
    @info "Saved: map_bottom_day_$(plot_day).png"
end

##########################################
#        record_horizontal_tracer(Ci_free, times, folder, "Ci_free", "Ci_free [μM N]",
#                                real_lon, real_lat;
#                                colorrange=(0, 1), colormap=Reverse(:cherry), iz=Nz)


if maps_animations_ 
# ------------------ Movies ------------------
    record_horizontal_tracer(S,   times, folder, "Ssurf",    "S [psu]",
                            real_lon, real_lat;
                            colorrange=(0, 35), colormap=:viridis, iz=Nz)
    record_horizontal_tracer(NUT, times, folder, "NUTsurf",  "NUT [μM N]",
                            real_lon, real_lat;
                            colorrange=(0, 40), colormap=Reverse(:cherry), iz=Nz)
    record_horizontal_tracer(DOM, times, folder, "DOMsurf",  "DOM [μM N]",
                            real_lon, real_lat;
                            colorrange=(0, 5), colormap=Reverse(:cherry), iz=Nz)
    record_horizontal_tracer(O₂,  times, folder, "O2surf",   "O₂ [μM]",
                            real_lon, real_lat;
                            colorrange=(0, 350), colormap=:turbo, iz=Nz)
    record_horizontal_tracer(P,   times, folder, "PHYsurf",  "PHY [μM N]",
                            real_lon, real_lat;
                            colorrange=(0, 5), colormap=Reverse(:cubehelix), iz=Nz)

    if Ci_
        record_horizontal_tracer(Ci_free, times, folder, "Ci_free", "Ci_free [μM N]",
                                real_lon, real_lat;
                                colorrange=(0, 1), colormap=Reverse(:cherry), iz=Nz)
        record_horizontal_tracer(Ci_PHY,  times, folder, "Ci_PHY",  "Ci_PHY [μM N]",
                                real_lon, real_lat;
                                colorrange=(0, 0.0001), colormap=Reverse(:cherry), iz=Nz)
        record_horizontal_tracer(Ci_HET,  times, folder, "Ci_HET",  "Ci_HET [μM N]",
                                real_lon, real_lat;
                                colorrange=(0, 0.00001), colormap=Reverse(:cherry), iz=Nz)
        record_horizontal_tracer(Ci_POM,  times, folder, "Ci_POM",  "Ci_POM [μM N]",
                                real_lon, real_lat;
                                colorrange=(0, 0.000001), colormap=Reverse(:cherry), iz=Nz)
        record_horizontal_tracer(Ci_DOM,  times, folder, "Ci_DOM",  "Ci_DOM [μM N]",
                                real_lon, real_lat;
                                colorrange=(0, 0.0000001), colormap=Reverse(:cherry), iz=Nz)
    end

    record_bottom_tracer(O₂, "O2_bottom", Nz, times, folder,
                        real_lon, real_lat;
                        colorrange=(-1, 350), colormap=:turbo, figsize=(400, 550))
end
##########################################
#-------Six animation in 1 figure---------
if six_animations_
    prepare_six_animations(tracers=[NUT, O₂, P, HET, POM, DOM],
                           times=times, folder=folder,
                           real_lon=real_lon, real_lat=real_lat,
                           Nz=Nz)
end

close(ds)
@info "Script completed!"
