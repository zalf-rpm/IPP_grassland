#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

from collections import defaultdict
import copy
import csv
from datetime import date, timedelta
import json
import math
import numpy as np
import os
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
from irrigation_manager import IrrigationManager

PATHS = {
    # adjust the local path to your environment
    "cj-local-remote": {
        # "include-file-base-path": "/home/berg/GitHub/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "D:/projects/KlimErtrag/",  # local path
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
        "path-debug-write-folder": "./debug-out/",
    },

    "remoteProducer-remoteMonica": {
        # "include-file-base-path": "/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "/data/",  # mounted path to archive or hard drive with climate data
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
DATA_GRID_SOIL_OW = "germany/buek200_1000_25832_etrs89-utm32n_OW.asc"
DATA_GRID_CROPS = "germany/OWgermany-crop-ww_1000_25832_etrs89-utm32n.asc"  # Added as a cropmap for winter wheat OW
# ORIGINAL DATA_GRID_SOIL = "germany/buek200_1000_25832_etrs89-utm32n.asc"
# DATA_GRID_CROPS = "germany/crops-all2017-2019_1000_25832_etrs89-utm32n.asc"
# DATA_GRID_CROPS = "germany/dwd-stations-pheno_1000_25832_etrs89-utm32n.asc"
# DATA_GRID_CROPS = "germany/germany-complete_1000_25832_etrs89-utm32n.asc"
DATA_GRID_IRRIGATION = "germany/irrigation_1000_25832_etrs89-utm32n_wc_18.asc"
TEMPLATE_PATH_LATLON = "{path_to_climate_dir}/latlon-to-rowcol.json"
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


## Add an argument in the run_producer function and make a loop with changing of the value of the additional parameter (sensitivity analysis)
## Make a list of the parameter values first

# commandline parameters e.g "server=localhost port=6666 shared_id=2"
def run_producer(server={"server": None, "port": None}, shared_id=None):
    context = zmq.Context()
    socket = context.socket(zmq.PUSH)  # pylint: disable=no-member
    # config_and_no_data_socket = context.socket(zmq.PUSH)

    config = {
        "mode": "mbm-local-remote",
        "server-port": server["port"] if server["port"] else "6666",
        "server": server["server"] if server["server"] else "login01.cluster.zalf.de",
        "start-row": "0",
        "end-row": "-1",
        "path_to_dem_grid": "",
        "sim.json": "sim.json",
        "crop.json": "crop.json",
        "site.json": "site.json",
        "setups-file": "sim_setups_VK.csv",
        "run-setups": "[1]",
        "shared_id": shared_id
    }

    # read commandline args only if script is invoked directly from commandline
    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=")
            if k in config:
                config[k] = v

    print("config:", config)

    # select paths 
    paths = PATHS[config["mode"]]
    # open soil db connection
    soil_db_con = sqlite3.connect(paths["path-to-data-dir"] + DATA_SOIL_DB)
    # soil_db_con = cas_sq3.connect(paths["path-to-data-dir"] + DATA_SOIL_DB) #CAS.
    # connect to monica proxy (if local, it will try to connect to a locally started monica)
    socket.connect("tcp://" + config["server"] + ":" + str(config["server-port"]))

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

    # irrigation data
    path_to_irrigation_grid = paths["path-to-data-dir"] + DATA_GRID_IRRIGATION
    irrigation_epsg_code = int(path_to_irrigation_grid.split("/")[-1].split("_")[2])
    irrigation_crs = CRS.from_epsg(irrigation_epsg_code)
    if irrigation_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[irrigation_crs] = Transformer.from_crs(soil_crs, irrigation_crs)
    irrigation_metadata, _ = Mrunlib.read_header(path_to_irrigation_grid)
    irrigation_grid = np.loadtxt(path_to_irrigation_grid, dtype=int, skiprows=6)
    irrigation_interpolate = Mrunlib.create_ascii_grid_interpolator(irrigation_grid, irrigation_metadata, False)
    print("read: ", path_to_irrigation_grid)

    # initialize irrigation manager
    irrigation_manager = IrrigationManager("irrigated_crops.json")

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

    sent_env_count = 0
    start_time = time.perf_counter()

    listOfClimateFiles = set()

    # run calculations for each setup
    for _, setup_id in enumerate(run_setups):

        if setup_id not in setups:
            continue
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

        if region_name and len(region_name) > 0:
            # Create the soil mask for the specific region
            path_to_soil_grid_ow = paths["path-to-data-dir"] + DATA_GRID_SOIL_OW
            mask = create_mask_from_shapefile(NUTS1_REGIONS, region_name, path_to_soil_grid_ow)

            # Apply the soil mask to the soil grid
            soil_grid_copy = soil_grid.copy()
            soil_grid[mask == False] = -8888
            soil_grid[soil_grid_copy == -9999] = -9999

        # add crop id from setup file
        try:
            # read seed/harvest dates for each crop_id
            path_harvest = TEMPLATE_PATH_HARVEST.format(path_to_data_dir=paths["path-to-data-dir"],
                                                        crop_id=crop_id_short)
            print("created seed harvest gk5 interpolator and read data: ", path_harvest)
            Mrunlib.create_seed_harvest_geoGrid_interpolator_and_read_data(path_harvest, wgs84_crs, utm32_crs,
                                                                           ilr_seed_harvest_data)
        except IOError:
            path_harvest = TEMPLATE_PATH_HARVEST.format(path_to_data_dir=paths["path-to-data-dir"],
                                                        crop_id=crop_id_short)
            print("Couldn't read file:", path_harvest)
            continue

        cdict = {}
        # path to latlon-to-rowcol.json
        # path = TEMPLATE_PATH_LATLON.format(path_to_climate_dir=paths["path-to-climate-dir"] + setup["climate_path_to_latlon_file"] + "/")
        path = TEMPLATE_PATH_LATLON.format(
            path_to_climate_dir=paths["path-to-climate-dir"] + setup["climate_path_to_latlon_file"] + "/")
        climate_data_interpolator = Mrunlib.create_climate_geoGrid_interpolator_from_json_file(path, wgs84_crs,
                                                                                               soil_crs, cdict)
        print("created climate_data to gk5 interpolator: ", path)

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

        if len(scenario) > 0 and scenario[:3].lower() == "rcp":
            site_json["EnvironmentParameters"]["rcp"] = scenario

        # read template crop.json
        with open(setup.get("crop.json", config["crop.json"])) as _:
            crop_json = json.load(_)

        crop_json["CropParameters"]["__enable_vernalisation_factor_fix__"] = setup[
            "use_vernalisation_fix"] if "use_vernalisation_fix" in setup else False

        # set the current crop used for this run id
        crop_json["cropRotation"][2] = crop_id

        # create environment template from json templates
        env_template = monica_io3.create_env_json_from_json_config({
            "crop": crop_json,
            "site": site_json,
            "sim": sim_json,
            "climate": ""
        })

        # set shared id in template
        if config["shared_id"]:
            env_template["sharedId"] = config["shared_id"]

        scols = int(soil_metadata["ncols"])
        srows = int(soil_metadata["nrows"])
        scellsize = int(soil_metadata["cellsize"])
        xllcorner = int(soil_metadata["xllcorner"])
        yllcorner = int(soil_metadata["yllcorner"])
        nodata_value = int(soil_metadata["nodata_value"])

        # unknown_soil_ids = set()
        soil_id_cache = {}
        print("All Rows x Cols: " + str(srows) + "x" + str(scols))
        # cs__ = open("coord_mapping_etrs89-utm32n_to_wgs84-latlon.csv", "w")
        # cs__.write("row,col,center_25832_etrs89-utm32n_r,center_25832_etrs89-utm32n_h,center_lat,center_lon\n")

        for srow in range(0, srows):
            print(srow, end=", ")

            if srow < int(config["start-row"]):
                continue
            elif int(config["end-row"]) > 0 and srow > int(config["end-row"]):
                break

            for scol in range(0, scols):
                soil_id = int(soil_grid[srow, scol])
                if soil_id == nodata_value:
                    continue

                # get coordinate of clostest climate element of real soil-cell
                sh = yllcorner + (scellsize / 2) + (srows - srow - 1) * scellsize
                sr = xllcorner + (scellsize / 2) + scol * scellsize
                # inter = crow/ccol encoded into integer
                crow, ccol = climate_data_interpolator(sr, sh)

                crop_grid_id = int(crop_grid[srow, scol])
                # print(crop_grid_id)
                if crop_grid_id != 1 or soil_id == -8888:
                    # print("row/col:", srow, "/", scol, "is not a crop pixel.")
                    env_template["customId"] = {
                        "setup_id": setup_id,
                        "srow": srow, "scol": scol,
                        "crow": int(crow), "ccol": int(ccol),
                        "soil_id": soil_id,
                        "env_id": sent_env_count,
                        "nodata": True,
                    }
                    if not DEBUG_DONOT_SEND:
                        socket.send_json(env_template)
                        # print("sent nodata env ", sent_env_count, " customId: ", env_template["customId"])
                        sent_env_count += 1
                    continue

                tcoords = {}

                """
                lon, lat = soil_crs_to_x_transformers[wgs84_crs].transform(sr, sh)
                try:
                    int(lon)
                    int(lat)
                except Exception as e:
                    lon, lat = wgs84_ip(sr, sh)

                cs__.write(str(srow) + "," + str(scol) + "," + str(sr) + "," + str(sh) + "," + str(lat) + "," + str(lon) + "\n")
                continue
                """

                if soil_id in soil_id_cache:
                    soil_profile = soil_id_cache[soil_id]
                else:
                    soil_profile = soil_io3.soil_parameters(soil_db_con, soil_id)
                    soil_id_cache[soil_id] = soil_profile

                worksteps = env_template["cropRotation"][0]["worksteps"]
                sowing_ws = next(filter(lambda ws: ws["type"][-6:] == "Sowing", worksteps))
                harvest_ws = next(filter(lambda ws: ws["type"][-7:] == "Harvest", worksteps))

                ilr_interpolate = ilr_seed_harvest_data[crop_id_short]["interpolate"]
                seed_harvest_cs = ilr_interpolate(sr, sh) if ilr_interpolate else None

                # set external seed/harvest dates
                if seed_harvest_cs:
                    seed_harvest_data = ilr_seed_harvest_data[crop_id_short]["data"][seed_harvest_cs]
                    if seed_harvest_data:
                        is_winter_crop = ilr_seed_harvest_data[crop_id_short]["is-winter-crop"]

                        if setup[
                            "sowing-date"] == "fixed":  # fixed indicates that regionally fixed sowing dates will be used
                            sowing_date = seed_harvest_data["sowing-date"]
                        elif setup[
                            "sowing-date"] == "auto":  # auto indicates that automatic sowng dates will be used that vary between regions
                            sowing_date = seed_harvest_data["latest-sowing-date"]
                        elif setup[
                            "sowing-date"] == "fixed1":  # fixed1 indicates that a fixed sowing date will be used that is the same for entire germany
                            sowing_date = sowing_ws["date"]

                        sds = [int(x) for x in sowing_date.split("-")]
                        sd = date(2001, sds[1], sds[2])
                        sdoy = sd.timetuple().tm_yday

                        if setup[
                            "harvest-date"] == "fixed":  # fixed indicates that regionally fixed harvest dates will be used
                            harvest_date = seed_harvest_data["harvest-date"]
                        elif setup[
                            "harvest-date"] == "auto":  # auto indicates that automatic harvest dates will be used that vary between regions
                            harvest_date = seed_harvest_data["latest-harvest-date"]
                        elif setup[
                            "harvest-date"] == "auto1":  # fixed1 indicates that a fixed harvest date will be used that is the same for entire germany
                            harvest_date = harvest_ws["latest-date"]

                        # print("sowing_date:", sowing_date, "harvest_date:", harvest_date)
                        # print("sowing_date:", sowing_ws["date"], "harvest_date:", sowing_ws["date"])

                        hds = [int(x) for x in harvest_date.split("-")]
                        hd = date(2001, hds[1], hds[2])
                        hdoy = hd.timetuple().tm_yday

                        esds = [int(x) for x in seed_harvest_data["earliest-sowing-date"].split("-")]
                        esd = date(2001, esds[1], esds[2])

                        # sowing after harvest should probably never occur in both fixed setup!
                        if setup["sowing-date"] == "fixed" and setup["harvest-date"] == "fixed":
                            # calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy-1))
                            if is_winter_crop:
                                calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy - 1))
                            else:
                                calc_harvest_date = date(2000, 12, 31) + timedelta(days=hdoy)
                            sowing_ws["date"] = seed_harvest_data["sowing-date"]
                            harvest_ws["date"] = "{:04d}-{:02d}-{:02d}".format(hds[0], calc_harvest_date.month,
                                                                               calc_harvest_date.day)
                            print("dates: ", int(seed_harvest_cs), ":", sowing_ws["date"])
                            print("dates: ", int(seed_harvest_cs), ":", harvest_ws["date"])

                        elif setup["sowing-date"] == "fixed" and setup["harvest-date"] == "auto":
                            if is_winter_crop:
                                calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy - 1))
                            else:
                                calc_harvest_date = date(2000, 12, 31) + timedelta(days=hdoy)
                            sowing_ws["date"] = seed_harvest_data["sowing-date"]
                            harvest_ws["latest-date"] = "{:04d}-{:02d}-{:02d}".format(hds[0], calc_harvest_date.month,
                                                                                      calc_harvest_date.day)
                            print("dates: ", int(seed_harvest_cs), ":", sowing_ws["date"])
                            print("dates: ", int(seed_harvest_cs), ":", harvest_ws["latest-date"])

                        elif setup["sowing-date"] == "fixed" and setup["harvest-date"] == "auto1":
                            if is_winter_crop:
                                calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy - 1))
                            else:
                                calc_harvest_date = date(2000, 12, 31) + timedelta(days=hdoy)
                            sowing_ws["date"] = seed_harvest_data["sowing-date"]
                            harvest_ws["latest-date"] = "{:04d}-{:02d}-{:02d}".format(hds[0], hds[1], hds[2])
                            print("dates: ", int(seed_harvest_cs), ":", sowing_ws["date"])
                            print("dates: ", int(seed_harvest_cs), ":", harvest_ws["latest-date"])

                        elif setup["sowing-date"] == "auto" and setup["harvest-date"] == "fixed":
                            sowing_ws["earliest-date"] = seed_harvest_data["earliest-sowing-date"] if esd > date(
                                esd.year, 6, 20) else "{:04d}-{:02d}-{:02d}".format(sds[0], 6, 20)
                            calc_sowing_date = date(2000, 12, 31) + timedelta(days=max(hdoy + 1, sdoy))
                            sowing_ws["latest-date"] = "{:04d}-{:02d}-{:02d}".format(sds[0], calc_sowing_date.month,
                                                                                     calc_sowing_date.day)
                            harvest_ws["date"] = seed_harvest_data["harvest-date"]
                            print("dates: ", int(seed_harvest_cs), ":", sowing_ws["earliest-date"], "<",
                                  sowing_ws["latest-date"])
                            print("dates: ", int(seed_harvest_cs), ":", harvest_ws["date"])

                        elif setup["sowing-date"] == "auto" and setup["harvest-date"] == "auto":
                            sowing_ws["earliest-date"] = seed_harvest_data["earliest-sowing-date"] if esd > date(
                                esd.year, 6, 20) else "{:04d}-{:02d}-{:02d}".format(sds[0], 6, 20)
                            if is_winter_crop:
                                calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy - 1))
                            else:
                                calc_harvest_date = date(2000, 12, 31) + timedelta(days=hdoy)
                            sowing_ws["latest-date"] = seed_harvest_data["latest-sowing-date"]
                            harvest_ws["latest-date"] = "{:04d}-{:02d}-{:02d}".format(hds[0], calc_harvest_date.month,
                                                                                      calc_harvest_date.day)
                            print("dates: ", int(seed_harvest_cs), ":", sowing_ws["earliest-date"], "<",
                                  sowing_ws["latest-date"])
                            print("dates: ", int(seed_harvest_cs), ":", harvest_ws["latest-date"])

                        elif setup["sowing-date"] == "fixed1" and setup["harvest-date"] == "fixed":
                            # calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy-1))
                            if is_winter_crop:
                                calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy - 1))
                            else:
                                calc_harvest_date = date(2000, 12, 31) + timedelta(days=hdoy)
                            sowing_ws["date"] = sowing_date
                            # print(seed_harvest_data["sowing-date"])
                            harvest_ws["date"] = "{:04d}-{:02d}-{:02d}".format(hds[0], calc_harvest_date.month,
                                                                               calc_harvest_date.day)
                            print("dates: ", int(seed_harvest_cs), ":", sowing_ws["date"])
                            print("dates: ", int(seed_harvest_cs), ":", harvest_ws["date"])

                    # print("dates: ", int(seed_harvest_cs), ":", sowing_ws["earliest-date"], "<", sowing_ws["latest-date"] )
                    # print("dates: ", int(seed_harvest_cs), ":", harvest_ws["latest-date"], "<", sowing_ws["earliest-date"], "<", sowing_ws["latest-date"] )

                    # print("dates: ", int(seed_harvest_cs), ":", sowing_ws["date"])
                    # print("dates: ", int(seed_harvest_cs), ":", harvest_ws["date"])

                if len(soil_profile) == 0:
                    # print("row/col:", srow, "/", scol, "has unknown soil_id:", soil_id)
                    # unknown_soil_ids.add(soil_id)

                    env_template["customId"] = {
                        "setup_id": setup_id,
                        "srow": srow, "scol": scol,
                        "crow": int(crow), "ccol": int(ccol),
                        "soil_id": soil_id,
                        "env_id": sent_env_count,
                        "nodata": True,
                    }
                    if and not DEBUG_DONOT_SEND:
                        socket.send_json(env_template)
                        # print("sent nodata env ", sent_env_count, " customId: ", env_template["customId"])
                        sent_env_count += 1
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

                if irrigation_crs not in tcoords:
                    tcoords[irrigation_crs] = soil_crs_to_x_transformers[irrigation_crs].transform(sr, sh)
                irr_r, irr_h = tcoords[irrigation_crs]
                irrigation = int(irrigation_interpolate(irr_r, irr_h))

                env_template["params"]["userCropParameters"]["__enable_T_response_leaf_expansion__"] = setup[
                    "LeafExtensionModifier"]

                # print("soil:", soil_profile)
                env_template["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

                # setting groundwater level
                if setup["groundwater-level"]:
                    groundwaterlevel = 20
                    layer_depth = 0
                    for layer in soil_profile:
                        if layer.get("is_in_groundwater", False):
                            groundwaterlevel = layer_depth
                            # print("setting groundwaterlevel of soil_id:", str(soil_id), "to", groundwaterlevel, "m")
                            break
                        layer_depth += Mrunlib.get_value(layer["Thickness"])
                    env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
                    env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [
                        max(0, groundwaterlevel - 0.2), "m"]
                    env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [
                        groundwaterlevel + 0.2, "m"]

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
                    orig_stage_ts = env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"]["cultivar"][
                        "StageTemperatureSum"][0]
                    if len(stage_ts) != len(orig_stage_ts):
                        stage_ts = orig_stage_ts
                        print('The provided StageTemperatureSum array is not '
                              'sufficiently long. Falling back to original StageTemperatureSum')

                    env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"]["cultivar"][
                        "StageTemperatureSum"][0] = stage_ts

                env_template["params"]["simulationParameters"]["UseNMinMineralFertilisingMethod"] = setup[
                    "fertilization"]

                # set UseAutomaticIrrigation to True if irrigation setup is True and irrigation is 1
                if setup["irrigation"] and irrigation == 1:
                    # check if the crop type is in the irrigated crops map
                    if irrigation_manager.should_be_irrigated_by_crop_id()(setup["crop-id"]):
                        env_template["params"]["simulationParameters"]["UseAutomaticIrrigation"] = True
                        # add default values for irrigation amount and threshold
                        env_template["params"]["simulationParameters"]["AutoIrrigationParams"]["amount"] = [10, "mm"]
                        env_template["params"]["simulationParameters"]["AutoIrrigationParams"]["threshold"] = 0.3
                    else:
                        env_template["params"]["simulationParameters"]["UseAutomaticIrrigation"] = False
                        # reset irrigation amount and threshold
                        env_template["params"]["simulationParameters"]["AutoIrrigationParams"]["amount"] = [0, "mm"]
                        env_template["params"]["simulationParameters"]["AutoIrrigationParams"]["threshold"] = 0.9

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

                print("Harvest type:", setup["harvest-date"])
                print("Srow: ", env_template["customId"]["srow"], "Scol:", env_template["customId"]["scol"])
                harvest_ws = next(
                    filter(lambda ws: ws["type"][-7:] == "Harvest", env_template["cropRotation"][0]["worksteps"]))
                if setup["harvest-date"] == "fixed":
                    print("Harvest-date:", harvest_ws["date"])
                elif setup["harvest-date"] == "auto":
                    print("Harvest-date:", harvest_ws["latest-date"])

                if not DEBUG_DONOT_SEND:
                    socket.send_json(env_template)
                    print("sent env ", sent_env_count, " customId: ", env_template["customId"])

                sent_env_count += 1

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

        if env_template and is_sensitivity_analysis:
            env_template["pathToClimateCSV"] = ""
            env_template["customId"] = {
                "setup_id": setup_id,
                "no_of_sent_envs": sent_env_count,
            }
            socket.send_json(env_template)

            # print("crows/cols:", crows_cols)
        # cs__.close()
        stop_setup_time = time.perf_counter()
        print("\nSetup ", sent_env_count, " envs took ", (stop_setup_time - start_setup_time), " seconds")
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
