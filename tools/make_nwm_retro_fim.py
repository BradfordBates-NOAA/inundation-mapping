import os
import argparse
from timeit import default_timer as timer
import pandas as pd
import geopandas as gpd
from concurrent.futures import ProcessPoolExecutor, as_completed, wait

from inundate_mosaic_wrapper import produce_mosaicked_inundation
from get_nwm_max_flows_retro import compute_max_flows


def get_bounds_from_hucs(hucs, watershed_boundary_gpkg):
    
    print("Calculating extent from provided HUC8s...")
    # If the "hucs" argument is really one huc, convert it to a list.
    # Otherwise read into list from pandas
    if os.path.exists(hucs[0]):
        df = pd.read_csv(hucs[0])
        huc_list = df['huc8'].tolist()
    else:
        huc_list = hucs

    huc_list_zfill = []
    for huc in huc_list:
        huc_list_zfill.append(str(huc).zfill(8))

    if not os.path.exists(watershed_boundary_gpkg):
        print(watershed_boundary_gpkg + " doesn't exist.")
    
    # Load the GeoPackage
    gdf = gpd.read_file(watershed_boundary_gpkg, layer='WBDHU8')
    #print(gdf.crs)
    #print(gdf.columns)
    # Filter for features where 'huc8' is in 'huc_list'
    filtered_gdf = gdf[gdf['HUC8'].isin(huc_list_zfill)].to_crs("EPSG:4326")
    bounds = filtered_gdf.total_bounds

    return bounds


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
    branch_workers=1,
    huc_workers=1,
    remove_intermediate=True,
    verbose=False,
    is_mosaic_for_branches=False,
    watershed_boundary_gpkg=r'/data/inputs/wbd/WBD_National.gpkg'
):

    # Check job numbers and raise error if necessary
    total_cpus_requested = branch_workers * huc_workers
    total_cpus_available = os.cpu_count() - 1
    if total_cpus_requested > total_cpus_available:
        raise ValueError(
            'The HUC job number, {}, multiplied by the inundate job number, {}, '
            'exceeds your machine\'s available CPU count minus one: {}. '
            'Please lower the job_number_huc or job_number_inundate '
            'values accordingly.'.format(huc_workers, branch_workers, total_cpus_available)
        )

    #bounds = get_bounds_from_hucs(hucs, watershed_boundary_gpkg)
    bounds = []
    # Build list of HUCs to process
    if os.path.exists(hucs[0]):
        df = pd.read_csv(hucs[0])
        huc_list = df['huc8'].tolist()
    else:
        huc_list = hucs
    huc_list_zfill = []
    for huc in huc_list:
        huc_list_zfill.append(str(huc).zfill(8))

    # TODO allow for a dataframe to be returned instead of writing to CSV and reading back in
    flow_file = os.path.join(os.path.dirname(inundation_polygon), ((start_datetime + end_datetime).replace('-','_').replace(':','_').replace(' ','')) + '.csv')
    compute_max_flows(start_datetime, end_datetime, flow_file, bounds)

    inundation_polygons_to_merge = []
    inundation_tifs_to_merge = []
    depths_tifs_to_merge = []

    inundation_polygon_huc, inundation_raster_huc, depths_raster_huc = None, None, None

    
    print("Generating inundation map with max flow data...")
    with ProcessPoolExecutor(max_workers=huc_workers) as executor:
        for huc in huc_list_zfill:
            if '.gpkg' in inundation_polygon:
                inundation_polygon_huc = inundation_polygon.replace('.gpkg', f'_{huc}_mosaic.gpkg')
                inundation_polygons_to_merge.append(inundation_polygon_huc)
            
            if '.shp' in inundation_polygon:
                inundation_polygon_huc = inundation_polygon.replace('.shp', f'_{huc}_mosiac.shp')
                inundation_polygons_to_merge.append(inundation_polygon_huc)

            if inundation_raster != None:
                inundation_raster_huc = inundation_raster.replace('.tif', f'_{huc}_mosaic.tif')
                inundation_tifs_to_merge.append(inundation_raster_huc)
            
            if depths_raster != None:
                depths_raster_huc = depths_raster.replace('.tif', f'_{huc}_mosaic.tif')
                depths_tifs_to_merge.append(depths_raster_huc)

    #         produce_mosaicked_inundation(
    #                 hydrofabric_dir,
    #                 [huc],
    #                 flow_file,
    #                 inundation_raster_huc,
    #                 inundation_polygon_huc,
    #                 depths_raster_huc,
    #                 map_filename,
    #                 mask,
    #                 unit_attribute_name,
    #                 branch_workers,
    #                 remove_intermediate,
    #                 verbose,
    #                 is_mosaic_for_branches
    # )
                    
            executor.submit(
                produce_mosaicked_inundation,
                hydrofabric_dir,
                [huc],
                flow_file,
                inundation_raster_huc,
                inundation_polygon_huc,
                depths_raster_huc,
                map_filename,
                mask,
                unit_attribute_name,
                branch_workers,
                remove_intermediate,
                verbose,
                is_mosaic_for_branches
            )
    
    # Final step: merge all polygons, tifs, and depth grids into respective files.
    # Read each GeoPackage into a list of GeoDataFrames
    print("Merging all mosaicked polygons...")
    gdfs = [gpd.read_file(polygon) for polygon in inundation_polygons_to_merge]
    merged_gdf = pd.concat(gdfs, ignore_index=True)
    merged_gdf.to_file(inundation_polygon, driver='GPKG')

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
        "-u", "--hucs", help="List of HUCS to run, or a CSV with huc8 as the header.", required=True, default="", type=str, nargs="+"
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
    parser.add_argument("-bw", "--branch-workers", help="Number of concurrent branch jobs.", required=False, default=1, type=int)
    parser.add_argument("-hw", "--huc-workers", help="Number of concurrent HUCs to process.", required=False, default=1, type=int)
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
    parser.add_argument(
        "-wb",
        "--watershed-boundary-gpkg",
        help="WBD geopackage.",
        required=False,
        default=r'/data/inputs/wbd/WBD_National.gpkg',
    )

    # Extract to dictionary and run
    start = timer()
    make_nwm_retro_fim(**vars(parser.parse_args()))
    print(f"Completed in {round((timer() - start)/60, 2)} minutes.")
