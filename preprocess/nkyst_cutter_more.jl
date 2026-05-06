# nkyst_cutter.jl
using NCDatasets
using Dates

function subset_netcdf(infile::String, outfile::String; xlim, ylim)
    # ------------------------
    # 0. Ensure output directory exists
    # ------------------------
    output_dir = dirname(outfile)
    if !isdir(output_dir)
        println("Creating directory: $output_dir")
        mkpath(output_dir)
    end
    
    # Remove existing file if it exists
    if isfile(outfile)
        println("Removing existing file: $outfile")
        rm(outfile)
    end

    # ------------------------
    # 1. Open input dataset
    # ------------------------
    println("Opening input file: $infile")
    ds = NCDataset(infile, "r")

    # Read coordinate variables
    X = ds["X"][:]
    Y = ds["Y"][:]

    # Determine index ranges for subsetting
    ix = findall(xlim[1] .<= X .<= xlim[2])
    iy = findall(ylim[1] .<= Y .<= ylim[2])

    if isempty(ix) || isempty(iy)
        error("Subregion does not overlap with dataset domain.")
    end

    println("Cropping X indices: $(first(ix)):$(last(ix))")
    println("Cropping Y indices: $(first(iy)):$(last(iy))")

    # ------------------------
    # 2. Create output dataset
    # ------------------------
    println("Creating output file: $outfile")
    ds_out = NCDataset(outfile, "c")

    # Define dimensions, cropping X and Y FIRST
    defDim(ds_out, "X", length(ix))
    defDim(ds_out, "Y", length(iy))

    # Copy all other dims unchanged - BUT SKIP X and Y
    for (dname, dimlen) in ds.dim
        if dname != "X" && dname != "Y"
            defDim(ds_out, dname, dimlen)
        end
    end

    # ------------------------
    # 3. Copy coordinate variables
    # ------------------------
    x_out = defVar(ds_out, "X", Float64, ("X",))
    y_out = defVar(ds_out, "Y", Float64, ("Y",))

    # Copy attributes from input X variable
    for (att, v) in ds["X"].attrib
        x_out.attrib[att] = v
    end
    # Copy attributes from input Y variable
    for (att, v) in ds["Y"].attrib
        y_out.attrib[att] = v
    end

    # Write cropped coordinates
    x_out[:] = X[ix]
    y_out[:] = Y[iy]

    # ------------------------
    # 4. Copy and subset all remaining variables
    # ------------------------
    # Get all variable names
    varnames = keys(ds)
    
    for varname in varnames
        if varname in ("X", "Y")
            continue
        end
        
        var = ds[varname]
        
        # Get dimensions
        dims = dimnames(var)        
        
        println("Processing variable: $varname")
        println("  Dimensions: $dims")
        println("  eltype: $(eltype(var))")
        
        # SPECIAL HANDLING FOR TIME VARIABLE
        # Check if this is the time variable (by name or by DateTime type)
        if varname == "time" || eltype(var) == DateTime
            println("  Detected as time variable - special handling")
            
            # Time variables should be stored as Float64 in NetCDF
            # even though they display as DateTime
            v_out = defVar(ds_out, varname, Float64, dims)
            
            # Copy ALL attributes (especially units!)
            for (att, v) in var.attrib
                try
                    v_out.attrib[att] = v
                catch e
                    println("Warning: Could not copy attribute '$att' for time variable: $e")
                end
            end
            
            # Read the raw time values (they're Float64 internally)
            # Use the underlying array if possible
            time_data = var.var[:]  # Access the underlying Float64 array
            
            # Write to output
            v_out[:] = time_data
            
            println("  Time variable copied successfully")
            continue  # Skip the rest of the loop for time variable
        end
        
        # For non-time variables, proceed as before
        # Get the concrete data type
        T = eltype(var)
        if T isa Union
            # Extract the non-missing type
            non_missing_types = [t for t in Base.uniontypes(T) if t != Missing]
            if !isempty(non_missing_types)
                T = non_missing_types[1]
            else
                T = Float32
            end
        end
        
        # Define output variable with same type
        v_out = defVar(ds_out, varname, T, dims)

        # Copy variable attributes with proper type handling
        for (att, v) in var.attrib
            try
                # Special handling for numeric attributes that need type conversion
                if att in ("_FillValue", "missing_value", "valid_min", "valid_max", "valid_range")
                    # Convert to the variable's type
                    if v isa Number
                        v_out.attrib[att] = convert(T, v)
                    else
                        v_out.attrib[att] = v
                    end
                else
                    # Try to keep original type
                    v_out.attrib[att] = v
                end
            catch e
                println("Warning: Could not copy attribute '$att' for variable '$varname': $e")
                # Try string conversion as last resort
                try
                    v_out.attrib[att] = string(v)
                catch
                    println("Skipping attribute '$att'")
                end
            end
        end

        # Build slicing indices
        inds = Any[Colon() for _ in 1:length(dims)]
        for (i, d) in pairs(dims)
            if d == "X"
                inds[i] = ix
            elseif d == "Y"
                inds[i] = iy
            # else: keep as Colon()
            end
        end

        # Read and subset the variable
        data = var[inds...]
        
        # Handle missing values if present
        if eltype(data) <: Union{Missing, T}
            # Get fill value (use original if exists, otherwise default)
            fill_val = haskey(var.attrib, "_FillValue") ? var.attrib["_FillValue"] : 
                      T <: AbstractFloat ? T(NaN) : 
                      T <: Integer ? T(typemin(T)) : 
                      T(-9999)
            
            # Convert to appropriate type
            fill_val_converted = convert(T, fill_val)
            
            # Replace missing values with fill value
            data_converted = coalesce.(data, fill_val_converted)
            v_out[:] = data_converted
        else
            # No missing values, write directly
            v_out[:] = data
        end
    end

    # Copy global attributes
    for (att, v) in ds.attrib
        try
            ds_out.attrib[att] = v
        catch e
            println("Warning: Could not copy global attribute '$att': $e")
            try
                ds_out.attrib[att] = string(v)
            catch
                println("Skipping global attribute '$att'")
            end
        end
    end

    close(ds)
    close(ds_out)

    println("Subset written to $outfile")
end

# ------------------------
# Main execution
# ------------------------

# Define spatial subregion
    xlim = (340000, 436000) # Outer Oslofjord 
    ylim = (85000, 160000)  # Outer Oslofjord
 #   xlim = (390000, 408000) # Drammensfjorden from Outer Oslofjord
 #   ylim = (139000, 160000)  # Drammensfjorden from Outer Oslofjord
#Spcify date range

function run_loop()
    start_date = Date(2020, 01, 01)
    end_date =   Date(2020, 12, 31)

#----------------------------
    current_date = start_date

while current_date <= end_date
    
    println(current_date)
    #for date_of_file in 2020070100:100:2020073100
    #date_of_file = 2020013100
    #    folder = "/mnt/c/HOME/PROG/julia/oslofjord251203/"
    #    filename = "NorKyst-800m_ZDEPTHS_his.an.$date_of_file.nc"
        date_of_file = "$(Dates.format(current_date, "yyyymmdd"))"
        filename = "NKz_800m_$(date_of_file)00.nc"
    #    infile  = joinpath(folder, filename)
    #    infile  = joinpath(homedir(),"FjordSim_data", "oslofjord", "NorKyst", filename)    
    infile  = joinpath(homedir(),"src","gcp-chem-sim-private","data","input","NorKyst", filename)
        println("Input file: $infile")
    #    outfile = joinpath(homedir(),"FjordSim_data", "oslofjord", "NorKyst3", "NKz_$date_of_file.nc")
    outfile = joinpath(homedir(),"src","gcp-chem-sim-private","data","input","NorKyst3", "NKz_$date_of_file.nc")

    println("Starting subset operation...")
    subset_netcdf(infile, outfile; xlim=xlim, ylim=ylim)
    println("file: $filename")
    
    current_date += Day(1)
end
end
run_loop()
println("Done!")