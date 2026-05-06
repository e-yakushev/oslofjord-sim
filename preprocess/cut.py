#!/usr/bin/env python3
"""
Script to cut the last n records from a NetCDF file along the time dimension.
"""
import argparse
import os
from pathlib import Path

import xarray as xr

BASEDIR = Path(os.environ["PROJECT_ROOT"])


def cut_last_n_records(input_file, output_dir, n_records):
    """
    Cut the last n records from all variables along the time dimension.

    Parameters:
    -----------
    input_file : str or Path
        Path to input NetCDF file
    output_dir : str or Path
        Directory to save the output file
    n_records : int
        Number of records to keep from the end
    """
    # Open the dataset
    print(f"Opening {input_file}...")
    ds = xr.open_dataset(input_file)

    # Get the time dimension name (usually 'time' or 'ocean_time')
    time_dim = None
    for dim in ds.dims:
        if "time" in dim.lower():
            time_dim = dim
            break

    if time_dim is None:
        raise ValueError("Could not find time dimension in the dataset")

    total_records = ds.dims[time_dim]
    print(f"Total records in {time_dim}: {total_records}")
    print(f"Cutting last {n_records} records...")

    # Select the last n records
    ds_cut = ds.isel({time_dim: slice(-n_records, None)})

    print(f"Selected records: {ds_cut.dims[time_dim]}")

    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Create output filename
    input_path = Path(input_file)
    output_file = output_path / input_path.name

    # Save to NetCDF
    print(f"Saving to {output_file}...")
    ds_cut.to_netcdf(output_file)

    # Close datasets
    ds.close()
    ds_cut.close()

    print(f"Done! Output saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Cut the last n records from a NetCDF file along the time dimension"
    )
    parser.add_argument(
        "-i",
        "--input",
        default=str(BASEDIR / "data" / "output" / "snapshots_ocean_1.nc"),
        help="Input NetCDF file path",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=str(BASEDIR / "data" / "input"),
        help="Output directory path",
    )
    parser.add_argument(
        "-n",
        "--n-records",
        type=int,
        default=10,
        help="Number of records to keep from the end",
    )

    args = parser.parse_args()

    cut_last_n_records(args.input, args.output, args.n_records)
