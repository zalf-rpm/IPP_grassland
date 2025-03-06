#This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

from datetime import datetime
import capnp
from collections import defaultdict
import copy
import csv
from datetime import date, timedelta
import json
import math
import numpy as np
import os
from pathlib import Path
from pyproj import CRS, Transformer
import sqlite3
import sqlite3 as cas_sq3
import sys
import time
import zmq
import geopandas as gpd
import rasterio
from rasterio.transform import from_origin
from rasterio import features
import monica_io3
import soil_io3
import monica_run_lib as Mrunlib

PATHS = {
    # adjust the local path to your environment
    "re-local-remote": {
        # "include-file-base-path": "/home/berg/GitHub/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "/beegfs/common/data/climate/",
        # mounted path to archive or hard drive with climate data
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "./debug-out/",
    },
    # adjust the local path to your environmen
    "ow-local-remote": {
        # "include-file-base-path": "/home/berg/GitHub/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "/beegfs/common/data/climate/",
        # mounted path to archive or hard drive with climate data
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "./debug-out/",
    },
    "mbm-local-remote": {
        # "include-file-base-path": "/home/berg/GitHub/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive or hard drive with climate data
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "./debug-out/",
    },
    "mbm-local-local": {
        # "include-file-base-path": "/home/berg/GitHub/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive or hard drive with climate data
        "monica-path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
    },
    "hpc-local-remote": {
        # "path-to-climate-dir": "/beegfs/common/data/soil/global_soil_dataset_for_earth_system_modeling/",
        # mounted path to archive or hard drive with climate data
        "path-to-climate-dir": "/beegfs/common/data/climate/",  # mounted path to archive or hard drive with climate data
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "/out/debug-out/",
    }
}

DATA_SOIL_DB = "germany/buek200.sqlite"
DATA_GRID_HEIGHT = "germany/dem_1000_25832_etrs89-utm32n.asc"
DATA_GRID_SLOPE = "germany/slope_1000_25832_etrs89-utm32n.asc"
DATA_GRID_LAND_USE = "germany/landuse_1000_31469_gk5.asc"
DATA_GRID_SOIL = "germany/buek200_1000_25832_etrs89-utm32n.asc"
#DATA_GRID_SOIL_OW = "germany/buek200_1000_25832_etrs89-utm32n_OW.asc"
#DATA_GRID_CROPS = "germany/permanent-grass-mask-BB_1000_25832_etrs89-utm32n.asc"  #nodata was 0 in the file assigned below no data realigned to be -9999
DATA_GRID_CROPS = "germany/permanent-grass-mask-BB_1000_25832_etrs89-utm32n-realigned.asc" #nodata = -9999
#DATA_GRID_CROPS = "germany/OWgermany-crop-ww_1000_25832_etrs89-utm32n.asc"  # Added as a cropmap for winter wheat OW
# ORIGINAL DATA_GRID_SOIL = "germany/buek200_1000_25832_etrs89-utm32n.asc"
# DATA_GRID_CROPS = "germany/crops-all2017-2019_1000_25832_etrs89-utm32n.asc"
# DATA_GRID_CROPS = "germany/dwd-stations-pheno_1000_25832_etrs89-utm32n.asc"
# DATA_GRID_CROPS = "germany/germany-complete_1000_25832_etrs89-utm32n.asc"
# DATA_GRID_IRRIGATION = "germany/irrigation_1000_25832_etrs89-utm32n_wc_18.asc"
DATA_GRID_GW_MIN = "germany/gwl-min_1000_25832_etrs89-utm32.asc"  # min groundwater level map
DATA_GRID_GW_MAX = "germany/gwl-max_1000_25832_etrs89-utm32.asc"  # max groundwater level map
DATA_GRID_GW_MEAN = "germany/gwl-mean_1000_25832_etrs89-utm32.asc"  # mean groundwater level map
TEMPLATE_PATH_LATLON = "{path_to_climate_dir}/latlon-to-rowcol.json"
# TEMPLATE_PATH_LATLON = "data/latlon-to-rowcol.json"
TEMPLATE_PATH_CLIMATE_CSV = "{gcm}/{rcm}/{scenario}/{ensmem}/{version}/row-{crow}/col-{ccol}.csv"

# Additional data for masking the regions
NUTS1_REGIONS = "data/germany/NUTS250_N1.shp"

TEMPLATE_PATH_HARVEST = "{path_to_data_dir}/projects/monica-germany/ILR_SEED_HARVEST_doys_{crop_id}.csv"

gdf = gpd.read_file(NUTS1_REGIONS)

DEBUG_DONOT_SEND = False
DEBUG_WRITE = False
DEBUG_ROWS = 10
DEBUG_WRITE_FOLDER = "./debug_out"
DEBUG_WRITE_CLIMATE = False

def update_config(config, argv, print_config=False, allow_new_keys=False):
    if len(argv) > 1:
        for arg in argv[1:]:
            kv = arg.split("=", maxsplit=1)
            if len(kv) < 2:
                continue
            k, v = kv
            if len(k) > 1 and k[:2] == "--":
                k = k[2:]
            if allow_new_keys or k in config:
                config[k] = v.lower() == "true" if v.lower() in ["true", "false"] else v
        if print_config:
            print(config)

def run_producer(server=None, port=None, channel_server=None, channel_port=None):
    context = zmq.Context()
    socket = context.socket(zmq.PUSH)  # pylint: disable=no-member
    # config_and_no_data_socket = context.socket(zmq.PUSH)

    config = {
        "mode": "mbm-local-remote",
        "port": port if port else "6666",
        "server": server if server else "login01.cluster.zalf.de",
        "channel-port": channel_port if channel_port else "9998",
        "channel-server": channel_server if channel_server else "localhost",#"login01.cluster.zalf.de",
        "start-row": "0",
        "end-row": "-1",
        "row_cols": "[[216,507]]",
        "path_to_dem_grid": "",
        "sim.json": "sim_calibration.json",
        "crop.json": "crop.json",
        "site.json": "site.json",
        "setups-file": "sim_setups_calibration_VK.csv",
        "run-setups": "[1]",
        "reader_sr": None,
        "path_to_out": "out/",
        "only_nuts3_region_ids": "[]",  # "[10]",
    }

    update_config(config, sys.argv, print_config=True, allow_new_keys=False)

    path_to_out_file = config["path_to_out"] + "/producer.out"
    if not os.path.exists(config["path_to_out"]):
        try:
            os.makedirs(config["path_to_out"])
        except OSError:
            print("run-calibration-producer.py: Couldn't create dir:", config["path_to_out"], "!")
    with open(path_to_out_file, "a") as _:
        _.write(f"config: {config}\n")

    with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
        _.write(f"{datetime.now()} start producer in producer\n") 

    nuts3_region_ids = json.loads(config["only_nuts3_region_ids"])

    # select paths 
    paths = PATHS[config["mode"]]
    # open soil db connection
    soil_db_con = sqlite3.connect(paths["path-to-data-dir"] + DATA_SOIL_DB)
    # soil_db_con = cas_sq3.connect(paths["path-to-data-dir"] + DATA_SOIL_DB) #CAS.
    # connect to monica proxy (if local, it will try to connect to a locally started monica)
    socket.connect("tcp://" + config["server"] + ":" + str(config["port"]))

    # read setup from csv file
    setups = Mrunlib.read_sim_setups(config["setups-file"])
    rs_ranges = config["run-setups"][1:-1].split(",")
    run_setups = []
    for rsr in rs_ranges:
        rs_r = rsr.split("-")
        if 1 < len(rs_r) <= 2:
            run_setups.extend(range(int(rs_r[0]), int(rs_r[1])+1))
        elif len(rs_r) == 1:
            run_setups.append(int(rs_r[0]))
    #run_setups = json.loads(config["run-setups"])
    print("read sim setups: ", config["setups-file"])

    with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
        _.write(f"{datetime.now()} setup read\n") 

    # transforms geospatial coordinates from one coordinate reference system to another
    # transform wgs84 into gk5
    soil_crs_to_x_transformers = {}
    wgs84_crs = CRS.from_epsg(4326)
    utm32_crs = CRS.from_epsg(25832)
    # transformers[wgs84] = Transformer.from_crs(wgs84_crs, gk5_crs, always_xy=True)

    ilr_seed_harvest_data = defaultdict(
        lambda: {"interpolate": None, "data": defaultdict(dict), "is-winter-crop": None})

    # Load grids
    ## note numpy is able to load from a compressed file, ending with .gz or .bz2

    # soil data
    path_to_soil_grid = paths["path-to-data-dir"] + DATA_GRID_SOIL
    soil_epsg_code = int(path_to_soil_grid.split("/")[-1].split("_")[2])
    soil_crs = CRS.from_epsg(soil_epsg_code)
    if wgs84_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[wgs84_crs] = Transformer.from_crs(soil_crs, wgs84_crs)
    soil_metadata, _ = Mrunlib.read_header(path_to_soil_grid)
    soil_grid = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)
    print("read: ", path_to_soil_grid)

    # height data for germany
    path_to_dem_grid = paths["path-to-data-dir"] + DATA_GRID_HEIGHT
    dem_epsg_code = int(path_to_dem_grid.split("/")[-1].split("_")[2])
    dem_crs = CRS.from_epsg(dem_epsg_code)
    if dem_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[dem_crs] = Transformer.from_crs(soil_crs, dem_crs)
    dem_metadata, _ = Mrunlib.read_header(path_to_dem_grid)
    dem_grid = np.loadtxt(path_to_dem_grid, dtype=float, skiprows=6)
    dem_interpolate = Mrunlib.create_ascii_grid_interpolator(dem_grid, dem_metadata)
    print("read: ", path_to_dem_grid)

    # slope data
    path_to_slope_grid = paths["path-to-data-dir"] + DATA_GRID_SLOPE
    slope_epsg_code = int(path_to_slope_grid.split("/")[-1].split("_")[2])
    slope_crs = CRS.from_epsg(slope_epsg_code)
    if slope_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[slope_crs] = Transformer.from_crs(soil_crs, slope_crs)
    slope_metadata, _ = Mrunlib.read_header(path_to_slope_grid)
    slope_grid = np.loadtxt(path_to_slope_grid, dtype=float, skiprows=6)
    slope_interpolate = Mrunlib.create_ascii_grid_interpolator(slope_grid, slope_metadata)
    print("read: ", path_to_slope_grid)

    # land use data
    path_to_landuse_grid = paths["path-to-data-dir"] + DATA_GRID_LAND_USE
    landuse_epsg_code = int(path_to_landuse_grid.split("/")[-1].split("_")[2])
    landuse_crs = CRS.from_epsg(landuse_epsg_code)
    if landuse_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[landuse_crs] = Transformer.from_crs(soil_crs, landuse_crs)
    landuse_meta, _ = Mrunlib.read_header(path_to_landuse_grid)
    landuse_grid = np.loadtxt(path_to_landuse_grid, dtype=int, skiprows=6)
    landuse_interpolate = Mrunlib.create_ascii_grid_interpolator(landuse_grid, landuse_meta)
    print("read: ", path_to_landuse_grid)

    # crop mask data
    path_to_crop_grid = paths["path-to-data-dir"] + DATA_GRID_CROPS
    crop_epsg_code = int(path_to_crop_grid.split("/")[-1].split("_")[2])
    crop_crs = CRS.from_epsg(crop_epsg_code)
    if crop_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[crop_crs] = Transformer.from_crs(soil_crs, crop_crs)
    crop_meta, _ = Mrunlib.read_header(path_to_crop_grid)
    crop_grid = np.loadtxt(path_to_crop_grid, dtype=int, skiprows=6)
    crop_interpolate = Mrunlib.create_ascii_grid_interpolator(crop_grid, crop_meta)
    print("read: ", path_to_crop_grid)

    # # irrigation data
    # path_to_irrigation_grid = paths["path-to-data-dir"] + DATA_GRID_IRRIGATION
    # irrigation_epsg_code = int(path_to_irrigation_grid.split("/")[-1].split("_")[2])
    # irrigation_crs = CRS.from_epsg(irrigation_epsg_code)
    # if irrigation_crs not in soil_crs_to_x_transformers:
    #     soil_crs_to_x_transformers[irrigation_crs] = Transformer.from_crs(soil_crs, irrigation_crs)
    # irrigation_metadata, _ = Mrunlib.read_header(path_to_irrigation_grid)
    # irrigation_grid = np.loadtxt(path_to_irrigation_grid, dtype=int, skiprows=6)
    # irrigation_interpolate = Mrunlib.create_ascii_grid_interpolator(irrigation_grid, irrigation_metadata, False)
    # print("read: ", path_to_irrigation_grid)

    # min groundwater level data
    path_to_gw_min_grid = paths["path-to-data-dir"] + DATA_GRID_GW_MIN
    gw_min_epsg_code = int(path_to_gw_min_grid.split("/")[-1].split("_")[2])
    gw_min_crs = CRS.from_epsg(gw_min_epsg_code)
    if gw_min_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[gw_min_crs] = Transformer.from_crs(soil_crs, gw_min_crs)
    gw_min_metadata, _ = Mrunlib.read_header(path_to_gw_min_grid)
    gw_min_grid = np.loadtxt(path_to_gw_min_grid, dtype=float, skiprows=6)
    gw_min_interpolate = Mrunlib.create_ascii_grid_interpolator(gw_min_grid, gw_min_metadata)
    print("read: ", path_to_gw_min_grid)

    # max groundwater level data
    path_to_gw_max_grid = paths["path-to-data-dir"] + DATA_GRID_GW_MAX
    gw_max_epsg_code = int(path_to_gw_max_grid.split("/")[-1].split("_")[2])
    gw_max_crs = CRS.from_epsg(gw_max_epsg_code)
    if gw_max_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[gw_max_crs] = Transformer.from_crs(soil_crs, gw_max_crs)
    gw_max_metadata, _ = Mrunlib.read_header(path_to_gw_max_grid)
    gw_max_grid = np.loadtxt(path_to_gw_max_grid, dtype=float, skiprows=6)
    gw_max_interpolate = Mrunlib.create_ascii_grid_interpolator(gw_max_grid, gw_max_metadata)
    print("read: ", path_to_gw_max_grid)

    # mean groundwater level data
    path_to_gw_mean_grid = paths["path-to-data-dir"] + DATA_GRID_GW_MEAN
    gw_mean_epsg_code = int(path_to_gw_mean_grid.split("/")[-1].split("_")[2])
    gw_mean_crs = CRS.from_epsg(gw_mean_epsg_code)
    if gw_mean_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[gw_mean_crs] = Transformer.from_crs(soil_crs, gw_mean_crs)
    gw_mean_metadata, _ = Mrunlib.read_header(path_to_gw_mean_grid)
    gw_mean_grid = np.loadtxt(path_to_gw_mean_grid, dtype=float, skiprows=6)
    gw_mean_interpolate = Mrunlib.create_ascii_grid_interpolator(gw_mean_grid, gw_mean_metadata)
    print("read: ", path_to_gw_mean_grid)

    with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
        _.write(f"{datetime.now()} grids load producer\n\n")

        # Create the function for the mask. This function will later use the additional column in a setup file!

    def create_mask_from_shapefile(NUTS1_REGIONS, region_name, path_to_soil_grid):
        regions_df = gpd.read_file(NUTS1_REGIONS)
        region = regions_df[regions_df["NUTS_NAME"] == region_name]

        # This is needed to read the transformation data correctly from the file. With the original opening it does not work
        with rasterio.open(path_to_soil_grid) as dataset:
            soil_grid = dataset.read(1)
            transform = dataset.transform

        rows, cols = soil_grid.shape
        mask = rasterio.features.geometry_mask([region.geometry.values[0]], out_shape=(rows, cols), transform=transform,
                                               invert=True)
        return mask

    listOfClimateFiles = set()
    start_time = time.perf_counter()

    if len(run_setups) > 1 and run_setups[0] not in setups:
        return
    else:
        setup_id = run_setups[0]

    channel = context.socket(zmq.PULL)
    channel.connect("tcp://" + config["channel-server"] + ":" + config["channel-port"])
    #conman = common.ConnectionManager()
    #reader = conman.try_connect(config["reader_sr"], cast_as=fbp_capnp.Channel.Reader, retry_secs=1)
    sent_env_count = 0
    while True:
        params = channel.recv_json() # keys: MaxAssimilationRate, AssimilateReallocation, RootPenetrationRate
        print("received params: ", params)
        # check for end of data from in port
        if params == "done":
            break

        with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
            _.write(f"{datetime.now()} connected\n")

        env_template = None
        start_setup_time = None
        if "only_nuts3_region_ids" in params:
            nuts3_region_ids = params["only_nuts3_region_ids"]
            del params["only_nuts3_region_ids"]

        start_setup_time = time.perf_counter()

        setup = setups[setup_id]
        gcm = setup["gcm"]
        rcm = setup["rcm"]
        scenario = setup["scenario"]
        ensmem = setup["ensmem"]
        version = setup["version"]
        crop_id = setup["crop-id"]
        region_name = setup["region_name"]

        ## extract crop_id from crop-id name that has possible an extenstion
        crop_id_short = crop_id.split('_')[0]

        #with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
        #    _.write(f"{datetime.now()} setup started producer\n")

        if region_name and len(region_name) > 0:
            # Create the soil mask for the specific region
            path_to_soil_grid_ow = paths["path-to-data-dir"] + DATA_GRID_SOIL
            mask = create_mask_from_shapefile(NUTS1_REGIONS, region_name, path_to_soil_grid_ow)

            # Apply the soil mask to the soil grid
            soil_grid_copy = soil_grid.copy()
            soil_grid[mask == False] = -8888
            soil_grid[soil_grid_copy == -9999] = -9999

        #with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
        #    _.write(f"{datetime.now()} crop added producer\n")

        cdict = {}
        # path to latlon-to-rowcol.json
        # path = TEMPLATE_PATH_LATLON.format(path_to_climate_dir=paths["path-to-climate-dir"] + setup["climate_path_to_latlon_file"] + "/")
        path = TEMPLATE_PATH_LATLON.format(
            path_to_climate_dir=paths["path-to-climate-dir"] + setup["climate_path_to_latlon_file"] + "/")
        climate_data_interpolator = Mrunlib.create_climate_geoGrid_interpolator_from_json_file(path, wgs84_crs,
                                                                                                      soil_crs, cdict)
        print("created climate_data to gk5 interpolator: ", path)

        #with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
        #    _.write(f"{datetime.now()} climate data read producer\n")

        # read template sim.json
        with open(setup.get("sim.json", config["sim.json"])) as _:
            sim_json = json.load(_)
        # change start and end date according to setup
        if setup["start_date"]:
            sim_json["climate.csv-options"]["start-date"] = str(setup["start_date"])
        if setup["end_date"]:
            sim_json["climate.csv-options"]["end-date"] = str(setup["end_date"])
            # sim_json["include-file-base-path"] = paths["include-file-base-path"]

        # read template site.json
        with open(setup.get("site.json", config["site.json"])) as _:
            site_json = json.load(_)

        #with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
        #    _.write(f"{datetime.now()} read site and sim json producer\n\n")

        if len(scenario) > 0 and scenario[:3].lower() == "rcp":
            site_json["EnvironmentParameters"]["rcp"] = scenario

        # read template crop.json
        with open(setup.get("crop.json", config["crop.json"])) as _:
            crop_json = json.load(_)
            ps = crop_json["crops"]["RYE"]["cropParams"]
            for pname, pval in params.items():
                if pname == "SpecificLeafArea":
                    ps["cultivar"][pname][0] *= pval
                    ps["cultivar"][pname][1] *= pval
                    ps["cultivar"][pname][2] *= pval
                    ps["cultivar"][pname][3] *= pval
                    ps["cultivar"][pname][4] *= pval
                    ps["cultivar"][pname][5] *= pval
                elif pname == "DroughtStressThreshold":
                    ps["cultivar"][pname][0] *= pval
                    ps["cultivar"][pname][1] *= pval
                    ps["cultivar"][pname][2] *= pval
                    ps["cultivar"][pname][3] *= pval
                    ps["cultivar"][pname][4] *= pval
                    ps["cultivar"][pname][5] *= pval
                elif pname == "StageKcFactor":
                    ps["cultivar"][pname][0] *= pval
                    ps["cultivar"][pname][1] *= pval
                    ps["cultivar"][pname][2] *= pval
                    ps["cultivar"][pname][3] *= pval
                    ps["cultivar"][pname][4] *= pval
                    ps["cultivar"][pname][5] *= pval
                elif pname == "CropSpecificMaxRootingDepth":
                    ps["cultivar"][pname] = pval

        crop_json["CropParameters"]["__enable_vernalisation_factor_fix__"] = setup[
            "use_vernalisation_fix"] if "use_vernalisation_fix" in setup else False

        # create environment template from json templates
        env_template = monica_io3.create_env_json_from_json_config({
            "crop": crop_json,
            "site": site_json,
            "sim": sim_json,
            "climate": ""
        })

        scols = int(soil_metadata["ncols"])
        srows = int(soil_metadata["nrows"])
        scellsize = int(soil_metadata["cellsize"])
        xllcorner = int(soil_metadata["xllcorner"])
        yllcorner = int(soil_metadata["yllcorner"])
        nodata_value = int(soil_metadata["nodata_value"])

        # unknown_soil_ids = set()
        soil_id_cache = {}
        #print("All Rows x Cols: " + str(srows) + "x" + str(scols))
        # cs__ = open("coord_mapping_etrs89-utm32n_to_wgs84-latlon.csv", "w")
        # cs__.write("row,col,center_25832_etrs89-utm32n_r,center_25832_etrs89-utm32n_h,center_lat,center_lon\n")

        row_cols = json.loads(config["row_cols"])
        for srow, scol in row_cols:
        #for srow in range(0, srows):
            #print(srow, end=", ")
            print("srow:", srow, "scol:", scol)

            #if srow < int(config["start-row"]):
            #    continue
            #elif int(config["end-row"]) > 0 and srow > int(config["end-row"]):
            #    break

            #for scol in range(0, scols):
            soil_id = int(soil_grid[srow, scol])
            if soil_id == nodata_value:
                continue

            # get coordinate of clostest climate element of real soil-cell
            sh = yllcorner + (scellsize / 2) + (srows - srow - 1) * scellsize
            sr = xllcorner + (scellsize / 2) + scol * scellsize
            # inter = crow/ccol encoded into integer
            crow, ccol = climate_data_interpolator(sr, sh)

            #crop_grid_id = int(crop_grid[srow, scol])
            # print(crop_grid_id)
            #if crop_grid_id != 1 or soil_id == -8888:
            #    continue

            tcoords = {}

            if soil_id in soil_id_cache:
                soil_profile = soil_id_cache[soil_id]
            else:
                #soil_profile = soil_io3.soil_parameters(soil_db_con, soil_id)
                soil_profile_group = soil_io3.get_soil_profile_group(soil_db_con, soil_id)
                soil_profile = None
                if len(soil_profile_group) > 0 and len(soil_profile_group[0]) > 0:
                    most_layers = {"layers": None, "no": 0}
                    for p in soil_profile_group[0][1]:
                        if p["id"] == 1:
                            soil_profile = p["layers"]
                            break
                        else:
                            if len(p["layers"]) > most_layers["no"]:
                                most_layers["layers"] = p["layers"]
                                most_layers["no"] = len(p["layers"])
                    if not soil_profile and most_layers["layers"]:
                        soil_profile = most_layers["layers"]
                    else:
                        continue
                else:
                    continue
                soil_id_cache[soil_id] = soil_profile
            if not soil_profile or len(soil_profile) == 0:
                continue

            # check if current grid cell is used for agriculture
            if setup["landcover"]:
                if landuse_crs not in tcoords:
                    tcoords[landuse_crs] = soil_crs_to_x_transformers[landuse_crs].transform(sr, sh)
                lur, luh = tcoords[landuse_crs]
                landuse_id = landuse_interpolate(lur, luh)
                if landuse_id not in [2, 3, 4]:
                    continue

            if dem_crs not in tcoords:
                tcoords[dem_crs] = soil_crs_to_x_transformers[dem_crs].transform(sr, sh)
            demr, demh = tcoords[dem_crs]
            height_nn = dem_interpolate(demr, demh)

            if slope_crs not in tcoords:
                tcoords[slope_crs] = soil_crs_to_x_transformers[slope_crs].transform(sr, sh)
            slr, slh = tcoords[slope_crs]
            slope = slope_interpolate(slr, slh)

            if gw_min_crs not in tcoords:
                tcoords[gw_min_crs] = soil_crs_to_x_transformers[gw_min_crs].transform(sr, sh)
            gw_min_r, gw_min_h = tcoords[gw_min_crs]
            min_groundwater_depth = gw_min_interpolate(gw_min_r, gw_min_h)

            if gw_max_crs not in tcoords:
                tcoords[gw_max_crs] = soil_crs_to_x_transformers[gw_max_crs].transform(sr, sh)
            gw_max_r, gw_max_h = tcoords[gw_max_crs]
            max_groundwater_depth = gw_max_interpolate(gw_max_r, gw_max_h)

            if gw_mean_crs not in tcoords:
                tcoords[gw_mean_crs] = soil_crs_to_x_transformers[gw_mean_crs].transform(sr, sh)
            gw_mean_r, gw_mean_h = tcoords[gw_mean_crs]
            mean_groundwater_depth = gw_mean_interpolate(gw_mean_r, gw_mean_h)

            # if irrigation_crs not in tcoords:
            #     tcoords[irrigation_crs] = soil_crs_to_x_transformers[irrigation_crs].transform(sr, sh)
            # irr_r, irr_h = tcoords[irrigation_crs]
            # irrigation = int(irrigation_interpolate(irr_r, irr_h))

            env_template["params"]["userCropParameters"]["__enable_T_response_leaf_expansion__"] = setup[
                "LeafExtensionModifier"]

            # print("soil:", soil_profile)
            env_template["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

            # setting groundwater level
            # if setup["groundwater-level"]:
            #     groundwaterlevel = 20
            #     layer_depth = 0
            #     for layer in soil_profile:
            #         if layer.get("is_in_groundwater", False):
            #             groundwaterlevel = layer_depth
            #             # print("setting groundwaterlevel of soil_id:", str(soil_id), "to", groundwaterlevel, "m")
            #             break
            #         layer_depth += Mrunlib.get_value(layer["Thickness"])
            #     env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
            #     env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [
            #         max(0, groundwaterlevel - 0.2), "m"]
            #     env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [
            #         groundwaterlevel + 0.2, "m"]

            if setup["groundwater-level"] == "MINMAX":
                # Assign min and max groundwater depths to the environment template
                env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
                env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [
                    min_groundwater_depth / 100, "m"]
                env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [
                    max_groundwater_depth / 100, "m"]
            # elif setup["groundwater-level"] == "MEAN":
            #     # Assign mean groundwater depth to the environment template
            #     env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
            #     env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [
            #         mean_groundwater_depth / 100, "m"]
            #     env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [
            #         mean_groundwater_depth / 100, "m"]
            elif setup["groundwater-level"] == "MIN":
                # Assign min and max groundwater depths to the environment template
                env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
                env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [
                    min_groundwater_depth / 100, "m"]
                env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [
                    min_groundwater_depth / 100, "m"]
            elif setup["groundwater-level"] == "MAX":
                # Assign min and max groundwater depths to the environment template
                env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
                env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [
                    max_groundwater_depth / 100, "m"]
                env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [
                    max_groundwater_depth / 100, "m"]

            # setting impenetrable layer
            if setup["impenetrable-layer"]:
                impenetrable_layer_depth = Mrunlib.get_value(
                    env_template["params"]["userEnvironmentParameters"]["LeachingDepth"])
                layer_depth = 0
                for layer in soil_profile:
                    if layer.get("is_impenetrable", False):
                        impenetrable_layer_depth = layer_depth
                        # print("setting leaching depth of soil_id:", str(soil_id), "to", impenetrable_layer_depth, "m")
                        break
                    layer_depth += Mrunlib.get_value(layer["Thickness"])
                env_template["params"]["userEnvironmentParameters"]["LeachingDepth"] = [impenetrable_layer_depth,
                                                                                        "m"]
                env_template["params"]["siteParameters"]["ImpenetrableLayerDepth"] = [impenetrable_layer_depth, "m"]

            if setup["elevation"]:
                env_template["params"]["siteParameters"]["heightNN"] = float(height_nn)

            if setup["slope"]:
                env_template["params"]["siteParameters"]["slope"] = slope / 100.0

            if setup["latitude"]:
                clat, _ = cdict[(crow, ccol)]
                env_template["params"]["siteParameters"]["Latitude"] = clat

            if setup["CO2"]:
                env_template["params"]["userEnvironmentParameters"]["AtmosphericCO2"] = float(setup["CO2"])

            if setup["O3"]:
                env_template["params"]["userEnvironmentParameters"]["AtmosphericO3"] = float(setup["O3"])

            if setup["FieldConditionModifier"]:
                env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"]["species"][
                    "FieldConditionModifier"] = float(setup["FieldConditionModifier"])

            if setup["StageTemperatureSum"]:
                stage_ts = setup["StageTemperatureSum"].split('_')
                stage_ts = [int(temp_sum) for temp_sum in stage_ts]
                orig_stage_ts = env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"]["cultivar"]["="][
                    "StageTemperatureSum"][0]
                if len(stage_ts) != len(orig_stage_ts):
                    stage_ts = orig_stage_ts
                    print('The provided StageTemperatureSum array is not '
                          'sufficiently long. Falling back to original StageTemperatureSum')

                env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"]["cultivar"]["="][
                    "StageTemperatureSum"][0] = stage_ts

            env_template["params"]["simulationParameters"]["UseNMinMineralFertilisingMethod"] = setup[
                "fertilization"]
            env_template["params"]["simulationParameters"]["UseAutomaticIrrigation"] = setup["irrigation"]

            env_template["params"]["simulationParameters"]["NitrogenResponseOn"] = setup["NitrogenResponseOn"]
            env_template["params"]["simulationParameters"]["WaterDeficitResponseOn"] = setup[
                "WaterDeficitResponseOn"]
            env_template["params"]["simulationParameters"]["EmergenceMoistureControlOn"] = setup[
                "EmergenceMoistureControlOn"]
            env_template["params"]["simulationParameters"]["EmergenceFloodingControlOn"] = setup[
                "EmergenceFloodingControlOn"]

            env_template["csvViaHeaderOptions"] = sim_json["climate.csv-options"]

            subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(gcm=gcm, rcm=rcm, scenario=scenario, ensmem=ensmem,
                                                              version=version, crow=str(crow), ccol=str(ccol))
            for _ in range(4):
                subpath_to_csv = subpath_to_csv.replace("//", "/")
            env_template["pathToClimateCSV"] = [
                paths["monica-path-to-climate-dir"] + setup["climate_path_to_csvs"] + "/" + subpath_to_csv]
            if setup["incl_hist"]:
                hist_subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(gcm=gcm, rcm=rcm, scenario="historical",
                                                                       ensmem=ensmem, version=version,
                                                                       crow=str(crow), ccol=str(ccol))
                for _ in range(4):
                    hist_subpath_to_csv = hist_subpath_to_csv.replace("//", "/")
                env_template["pathToClimateCSV"].insert(0, paths["monica-path-to-climate-dir"] + setup[
                    "climate_path_to_csvs"] + "/" + hist_subpath_to_csv)
            print("pathToClimateCSV:", env_template["pathToClimateCSV"])
            if DEBUG_WRITE_CLIMATE:
                listOfClimateFiles.add(subpath_to_csv)

            env_template["customId"] = {
                "setup_id": setup_id,
                "srow": srow, "scol": scol,
                "crow": int(crow), "ccol": int(ccol),
                "soil_id": soil_id,
                "env_id": sent_env_count,
                "nodata": False
            }

            sent_env_count += 1

            #with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
            #    _.write(f"{datetime.now()} Sending jobs out (producer)\n")

            if not DEBUG_DONOT_SEND:
                socket.send_json(env_template)
                print("sent env ", sent_env_count, " customId: ", env_template["customId"])

            #with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
            #    _.write(f"{datetime.now()} Ended jobs (producer)\n")

        # write debug output, as json file
        if DEBUG_WRITE:
            debug_write_folder = paths["path-debug-write-folder"]
            if not os.path.exists(debug_write_folder):
                os.makedirs(debug_write_folder)
            if sent_env_count < DEBUG_ROWS:

                path_to_debug_file = debug_write_folder + "/row_" + str(sent_env_count - 1) + "_" + str(
                    setup_id) + ".json"

                if not os.path.isfile(path_to_debug_file):
                    with open(path_to_debug_file, "w") as _:
                        _.write(json.dumps(env_template))
                else:
                    print("WARNING: Row ", (sent_env_count - 1), " already exists")
        # print("unknown_soil_ids:", unknown_soil_ids)

        with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
            _.write(f"{datetime.now()} Sending final last message (producer) \n")
        # send a last message will be just forwarded by monica to signify last
        if env_template:
            env_template["pathToClimateCSV"] = ""
            env_template["customId"] = {
                "no_of_sent_envs": sent_env_count,
                "nodata": True
            }
            socket.send_json(env_template)

            # print("crows/cols:", crows_cols)
        # cs__.close()
        stop_setup_time = time.perf_counter()
        print("\nSetup ", setup_id, ":", sent_env_count, " envs took ", (stop_setup_time - start_setup_time),
              " seconds")
        sent_env_count = 0

    stop_time = time.perf_counter()

    # write summary of used json files
    if DEBUG_WRITE_CLIMATE:
        debug_write_folder = paths["path-debug-write-folder"]
        if not os.path.exists(debug_write_folder):
            os.makedirs(debug_write_folder)

        path_to_climate_summary = debug_write_folder + "/climate_file_list" + ".csv"
        with open(path_to_climate_summary, "w") as _:
            _.write('\n'.join(listOfClimateFiles))

    try:
        print("sending ", (sent_env_count - 1), " envs took ", (stop_time - start_time), " seconds")
        # print("ran from ", start, "/", row_cols[start], " to ", end, "/", row_cols[end]
        print("exiting run_producer()")
    except Exception:
        raise


if __name__ == "__main__":
    run_producer()
