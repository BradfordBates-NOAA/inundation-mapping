import os
import argparse
import xarray as xr
import fsspec
import numpy as np
import dask

URL = 's3://noaa-nwm-retro-v2-zarr-pds'


def compute_max_flows(start_datetime, end_datetime, output_flow_file):

    if not os.access(os.path.dirname(output_flow_file), os.W_OK):
        print(f"You do not have permissions to write to the directory: " + os.path.dirname(output_flow_file))
        return

    print("Connecting to NWM Retro S3...")
    ds = xr.open_zarr(fsspec.get_mapper(URL, anon=True), consolidated=True)

    # Get max streamflow
    print("Extracting max flows...")
    period_streamflow = ds['streamflow'].sel(time=slice(start_datetime, end_datetime))
    print(f'Size of streamflow variable for provided time range: {period_streamflow.nbytes/1e9:.1f} GB')

    var_max = period_streamflow.max(dim='time').compute()

    # Assuming var_max is a pandas Series after the .compute()
    df = var_max.to_pandas().to_frame(name='discharge')

    print("Writing discharge file: " + output_flow_file)
    df.to_csv(output_flow_file)

    print("Max flows generation complete.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Produce CSV flow files from NWM 2.1 Retrospective.')
    parser.add_argument(
        '-s',
        '--start-datetime',
        help='The start date and time. Format: YYYY-MM-DD HH:MM',
        required=True,
    )
    parser.add_argument(
        '-e',
        '--end-datetime',
        help='The end date and time. Format: YYYY-MM-DD HH:MM',
        required=True,
    )
    parser.add_argument(
        '-o',
        '--output-flow-file',
        help='Full path to output flow file. A CSV with feature_id and discharge (in CMS) is written.',
        required=True,
    )
    args = vars(parser.parse_args())
    compute_max_flows(**args)