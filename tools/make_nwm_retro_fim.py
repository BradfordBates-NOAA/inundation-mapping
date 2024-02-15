import os
import argparse
from timeit import default_timer as timer

from inundate_mosaic_wrapper import produce_mosaicked_inundation
from get_nwm_max_flows_retro import compute_max_flows


def make_nwm_retro_fim(
    start_datetime,
    end_datetime,
    hydrofabric_dir,
    hucs,
    inundation_raster=None,
    inundation_polygon=None,
    depths_raster=None,
    map_filename=None,
    mask=None,
    unit_attribute_name="huc8",
    num_workers=1,
    remove_intermediate=True,
    verbose=False,
    is_mosaic_for_branches=False,
):
    # TODO allow for a dataframe to be returned instead of writing to CSV and reading back in
    flow_file = os.path.join(os.path.dirname(inundation_polygon), ((start_datetime + end_datetime).replace('-','_').replace(':','_').replace(' ','')) + '.csv')
    compute_max_flows(start_datetime, end_datetime, flow_file)

    print("Generating inundation map with max flow data...")
    produce_mosaicked_inundation(
        hydrofabric_dir,
        hucs,
        flow_file,
        inundation_raster,
        inundation_polygon,
        depths_raster,
        map_filename,
        mask,
        unit_attribute_name,
        num_workers,
        remove_intermediate,
        verbose,
        is_mosaic_for_branches
        )

    print("Inundation mapping for max flows over specified time period is complete.")
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Produce CSV flow files from NWM 2.1 Retrospective.')

    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Helpful utility to produce mosaicked inundation extents (raster and poly) and depths given a datetime and HUC."
    )

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
        "-y",
        "--hydrofabric_dir",
        help="Directory path to FIM hydrofabric by processing unit.",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-u", "--hucs", help="List of HUCS to run", required=True, default="", type=str, nargs="+"
    )
    parser.add_argument(
        "-i", "--inundation-raster", help="Inundation raster output.", required=False, default=None, type=str
    )
    parser.add_argument(
        "-p",
        "--inundation-polygon",
        help="Inundation polygon output. Only writes if designated.",
        required=True,
        default=None,
        type=str,
    )
    parser.add_argument(
        "-d",
        "--depths-raster",
        help="Depths raster output. Only writes if designated. Appends HUC code in batch mode.",
        required=False,
        default=None,
        type=str,
    )
    parser.add_argument(
        "-m",
        "--map-filename",
        help="Path to write output map file CSV (optional). Default is None.",
        required=False,
        default=None,
        type=str,
    )
    parser.add_argument("-k", "--mask", help="Name of mask file.", required=False, default=None, type=str)
    parser.add_argument(
        "-a",
        "--unit_attribute_name",
        help='Name of attribute column in map_file. Default is "huc8".',
        required=False,
        default="huc8",
        type=str,
    )
    parser.add_argument("-w", "--num-workers", help="Number of workers.", required=False, default=1, type=int)
    parser.add_argument(
        "-r",
        "--remove-intermediate",
        help="Remove intermediate products, i.e. individual branch inundation.",
        required=False,
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Verbose printing. Not tested.",
        required=False,
        default=False,
        action="store_true",
    )

    # Extract to dictionary and run
    start = timer()
    make_nwm_retro_fim(**vars(parser.parse_args()))
    print(f"Completed in {round((timer() - start)/60, 2)} minutes.")
