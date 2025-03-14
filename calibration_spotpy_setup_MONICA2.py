#!/usr/bin/python
# -*- coding: UTF-8
import uuid

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
from datetime import datetime
import geopandas as gpd
import json
import logging
import numpy as np
import os
from pathlib import Path
from pyproj import CRS, Transformer
import rasterio
import spotpy
import sqlite3
import time
import re
import zmq
from rasterio import features
import monica_io3
import soil_io3
import monica_run_lib

PATHS = {
    "mbm-local-remote": {
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "./debug-out/",
        "path-to-output-dir": "/out/out/",
        "path-to-csv-output-dir": "/out/csv-out/"
    },
    "mbm-local-local": {
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        "monica-path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "/out/debug-out/",
        "path-to-output-dir": "/out/out/",
        "path-to-csv-output-dir": "/out/csv-out/"
    },
    "hpc-local-remote": {
        "path-to-climate-dir": "/beegfs/common/data/climate/",  # mounted path to archive or hard drive with climate data
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "/out/debug-out/",
        "path-to-output-dir": "/out/out/",
        "path-to-csv-output-dir": "/out/csv-out/"
    }
}

DATA_SOIL_DB = "germany/buek200.sqlite"
DATA_GRID_HEIGHT = "germany/dem_1000_25832_etrs89-utm32n.asc"
DATA_GRID_SLOPE = "germany/slope_1000_25832_etrs89-utm32n.asc"
DATA_GRID_LAND_USE = "germany/landuse_1000_31469_gk5.asc"
DATA_GRID_SOIL = "germany/buek200_1000_25832_etrs89-utm32n.asc"
#DATA_GRID_CROPS = "germany/permanent-grass-mask-BB_1000_25832_etrs89-utm32n.asc"  #nodata was 0 in the file assigned below no data realigned to be -9999
DATA_GRID_CROPS = "germany/permanent-grass-mask-BB_1000_25832_etrs89-utm32n-realigned.asc" #nodata = -9999
DATA_GRID_GW_MIN = "germany/gwl-min_1000_25832_etrs89-utm32.asc"  # min groundwater level map
DATA_GRID_GW_MAX = "germany/gwl-max_1000_25832_etrs89-utm32.asc"  # max groundwater level map
DATA_GRID_GW_MEAN = "germany/gwl-mean_1000_25832_etrs89-utm32.asc"  # mean groundwater level map
TEMPLATE_PATH_LATLON = "{path_to_climate_dir}/latlon-to-rowcol.json"
# TEMPLATE_PATH_LATLON = "data/latlon-to-rowcol.json"
TEMPLATE_PATH_CLIMATE_CSV = "{gcm}/{rcm}/{scenario}/{ensmem}/{version}/row-{crow}/col-{ccol}.csv"

# Additional data for masking the regions
NUTS1_REGIONS = "data/germany/NUTS250_N1.shp"

DEBUG_DONOT_SEND = False
DEBUG_WRITE = False
DEBUG_ROWS = 10
DEBUG_WRITE_FOLDER = "./debug_out"
DEBUG_WRITE_CLIMATE = False


def flatten_biomasses_dict(biomasses: dict):
    flattend_dict = []
    for year in sorted(biomasses.keys()):
        flattend_dict.extend(biomasses[year])
    return flattend_dict

class spot_setup(object):
    def __init__(self, user_params, observations, monicas_host, monicas_in_port, monicas_out_port,
                 calib_row_cols, setup_id, setup, path_to_out, mode):
        self.env_template = None
        self.orig_crop_params = None
        self.user_params = user_params
        self.params = []
        self.observations = observations
        self.obs_flat_list = flatten_biomasses_dict(observations)
        self.path_to_out = path_to_out
        self.path_to_prod_out_file = path_to_out + "/spot_setup.out"
        self.context = zmq.Context()
        self.prod_socket = self.context.socket(zmq.PUSH)
        self.prod_socket.connect(f"tcp://{monicas_host}:{monicas_in_port}")
        self.cons_socket = self.context.socket(zmq.DEALER)
        self.shared_id = str(uuid.uuid4())
        self.cons_socket.setsockopt_string(zmq.ROUTING_ID, self.shared_id)
        self.cons_socket.RCVTIMEO = 60000
        self.cons_socket.connect(f"tcp://{monicas_host}:{monicas_out_port}")
        self.calib_row_cols = calib_row_cols
        self.path_to_prod_out_file = f"{self.path_to_out}/producer.out"
        self.path_to_cons_out_file = f"{self.path_to_out}/consumer.out"
        self.gdf = gpd.read_file(NUTS1_REGIONS)
        self.path = None
        self.setup_id = setup_id
        self.setup = setup
        self.paths = PATHS[mode]

        self.init_producer()

        if not os.path.exists(path_to_out):
            try:
                os.makedirs(path_to_out)
            except OSError:
                print("spot_setup.__init__: Couldn't create dir:", path_to_out, "!")

        with open(self.path_to_prod_out_file, "a") as _:
            _.write(f"observations: {self.observations}\n")
            _.write(f"obs_flat_list: {self.obs_flat_list}\n")

        for par in user_params:
            par_name = par["name"]
            if "array" in par:
                par["name"] = f"{par_name}_{par['array']}"  # spotpy does not allow two parameters to have the same name
                del par["array"]
            if "derive_function" not in par:  # spotpy does not care about derived params
                self.params.append(spotpy.parameter.Uniform(**par))

    def create_mask_from_shapefile(self, NUTS1_REGIONS, region_name, path_to_soil_grid):
        regions_df = gpd.read_file(NUTS1_REGIONS)
        region = regions_df[regions_df["NUTS_NAME"] == region_name]

        # This is needed to read the transformation data correctly from the file. With the original opening it does not work
        with rasterio.open(path_to_soil_grid) as dataset:
            soil_grid = dataset.read(1)
            transform = dataset.transform

        rows, cols = soil_grid.shape
        mask = rasterio.features.geometry_mask([region.geometry.values[0]], out_shape=(rows, cols),
                                               transform=transform,
                                               invert=True)
        return mask

    def init_producer(self):
        if not os.path.exists(self.path_to_out):
            try:
                os.makedirs(self.path_to_out)
            except OSError:
                print("run-calibration-producer.py: Couldn't create dir:", self.path_to_out, "!")
        with open(self.path_to_prod_out_file, "a") as _:
            pass
            #_.write(f"config: {config}\n")

        with open(f"{self.path_to_out}/spot_setup.out", "a") as _:
            _.write(f"{datetime.now()} start producer in producer\n")

        # open soil db connection
        self.soil_db_con = sqlite3.connect(self.paths["path-to-data-dir"] + DATA_SOIL_DB)

        # transforms geospatial coordinates from one coordinate reference system to another
        # transform wgs84 into gk5
        self.soil_crs_to_x_transformers = {}
        self.wgs84_crs = CRS.from_epsg(4326)
        utm32_crs = CRS.from_epsg(25832)
        # transformers[wgs84] = Transformer.from_crs(wgs84_crs, gk5_crs, always_xy=True)

        # soil data
        path_to_soil_grid = self.paths["path-to-data-dir"] + DATA_GRID_SOIL
        soil_epsg_code = int(path_to_soil_grid.split("/")[-1].split("_")[2])
        self.soil_crs = CRS.from_epsg(soil_epsg_code)
        if self.wgs84_crs not in self.soil_crs_to_x_transformers:
            self.soil_crs_to_x_transformers[self.wgs84_crs] = Transformer.from_crs(self.soil_crs, self.wgs84_crs)
        soil_metadata, _ = monica_run_lib.read_header(path_to_soil_grid)
        self.soil_grid = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)
        print("read: ", path_to_soil_grid)
        self.scols = int(soil_metadata["ncols"])
        self.srows = int(soil_metadata["nrows"])
        self.scellsize = int(soil_metadata["cellsize"])
        self.xllcorner = int(soil_metadata["xllcorner"])
        self.yllcorner = int(soil_metadata["yllcorner"])
        self.nodata_value = int(soil_metadata["nodata_value"])


        # height data for germany
        path_to_dem_grid = self.paths["path-to-data-dir"] + DATA_GRID_HEIGHT
        dem_epsg_code = int(path_to_dem_grid.split("/")[-1].split("_")[2])
        self.dem_crs = CRS.from_epsg(dem_epsg_code)
        if self.dem_crs not in self.soil_crs_to_x_transformers:
            self.soil_crs_to_x_transformers[self.dem_crs] = Transformer.from_crs(self.soil_crs, self.dem_crs)
        dem_metadata, _ = monica_run_lib.read_header(path_to_dem_grid)
        dem_grid = np.loadtxt(path_to_dem_grid, dtype=float, skiprows=6)
        self.dem_interpolate = monica_run_lib.create_ascii_grid_interpolator(dem_grid, dem_metadata)
        print("read: ", path_to_dem_grid)

        # slope data
        path_to_slope_grid = self.paths["path-to-data-dir"] + DATA_GRID_SLOPE
        slope_epsg_code = int(path_to_slope_grid.split("/")[-1].split("_")[2])
        self.slope_crs = CRS.from_epsg(slope_epsg_code)
        if self.slope_crs not in self.soil_crs_to_x_transformers:
            self.soil_crs_to_x_transformers[self.slope_crs] = Transformer.from_crs(self.soil_crs, self.slope_crs)
        slope_metadata, _ = monica_run_lib.read_header(path_to_slope_grid)
        slope_grid = np.loadtxt(path_to_slope_grid, dtype=float, skiprows=6)
        self.slope_interpolate = monica_run_lib.create_ascii_grid_interpolator(slope_grid, slope_metadata)
        print("read: ", path_to_slope_grid)

        # land use data
        path_to_landuse_grid = self.paths["path-to-data-dir"] + DATA_GRID_LAND_USE
        landuse_epsg_code = int(path_to_landuse_grid.split("/")[-1].split("_")[2])
        self.landuse_crs = CRS.from_epsg(landuse_epsg_code)
        if self.landuse_crs not in self.soil_crs_to_x_transformers:
            self.soil_crs_to_x_transformers[self.landuse_crs] = Transformer.from_crs(self.soil_crs, self.landuse_crs)
        landuse_meta, _ = monica_run_lib.read_header(path_to_landuse_grid)
        landuse_grid = np.loadtxt(path_to_landuse_grid, dtype=int, skiprows=6)
        self.landuse_interpolate = monica_run_lib.create_ascii_grid_interpolator(landuse_grid, landuse_meta)
        print("read: ", path_to_landuse_grid)

        # crop mask data
        path_to_crop_grid = self.paths["path-to-data-dir"] + DATA_GRID_CROPS
        crop_epsg_code = int(path_to_crop_grid.split("/")[-1].split("_")[2])
        self.crop_crs = CRS.from_epsg(crop_epsg_code)
        if self.crop_crs not in self.soil_crs_to_x_transformers:
            self.soil_crs_to_x_transformers[self.crop_crs] = Transformer.from_crs(self.soil_crs, self.crop_crs)
        crop_meta, _ = monica_run_lib.read_header(path_to_crop_grid)
        crop_grid = np.loadtxt(path_to_crop_grid, dtype=int, skiprows=6)
        self.crop_interpolate = monica_run_lib.create_ascii_grid_interpolator(crop_grid, crop_meta)
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
        path_to_gw_min_grid = self.paths["path-to-data-dir"] + DATA_GRID_GW_MIN
        gw_min_epsg_code = int(path_to_gw_min_grid.split("/")[-1].split("_")[2])
        self.gw_min_crs = CRS.from_epsg(gw_min_epsg_code)
        if self.gw_min_crs not in self.soil_crs_to_x_transformers:
            self.soil_crs_to_x_transformers[self.gw_min_crs] = Transformer.from_crs(self.soil_crs, self.gw_min_crs)
        gw_min_metadata, _ = monica_run_lib.read_header(path_to_gw_min_grid)
        gw_min_grid = np.loadtxt(path_to_gw_min_grid, dtype=float, skiprows=6)
        self.gw_min_interpolate = monica_run_lib.create_ascii_grid_interpolator(gw_min_grid, gw_min_metadata)
        print("read: ", path_to_gw_min_grid)

        # max groundwater level data
        path_to_gw_max_grid = self.paths["path-to-data-dir"] + DATA_GRID_GW_MAX
        gw_max_epsg_code = int(path_to_gw_max_grid.split("/")[-1].split("_")[2])
        self.gw_max_crs = CRS.from_epsg(gw_max_epsg_code)
        if self.gw_max_crs not in self.soil_crs_to_x_transformers:
            self.soil_crs_to_x_transformers[self.gw_max_crs] = Transformer.from_crs(self.soil_crs, self.gw_max_crs)
        gw_max_metadata, _ = monica_run_lib.read_header(path_to_gw_max_grid)
        gw_max_grid = np.loadtxt(path_to_gw_max_grid, dtype=float, skiprows=6)
        self.gw_max_interpolate = monica_run_lib.create_ascii_grid_interpolator(gw_max_grid, gw_max_metadata)
        print("read: ", path_to_gw_max_grid)

        # mean groundwater level data
        path_to_gw_mean_grid = self.paths["path-to-data-dir"] + DATA_GRID_GW_MEAN
        gw_mean_epsg_code = int(path_to_gw_mean_grid.split("/")[-1].split("_")[2])
        self.gw_mean_crs = CRS.from_epsg(gw_mean_epsg_code)
        if self.gw_mean_crs not in self.soil_crs_to_x_transformers:
            self.soil_crs_to_x_transformers[self.gw_mean_crs] = Transformer.from_crs(self.soil_crs, self.gw_mean_crs)
        gw_mean_metadata, _ = monica_run_lib.read_header(path_to_gw_mean_grid)
        gw_mean_grid = np.loadtxt(path_to_gw_mean_grid, dtype=float, skiprows=6)
        self.gw_mean_interpolate = monica_run_lib.create_ascii_grid_interpolator(gw_mean_grid, gw_mean_metadata)
        print("read: ", path_to_gw_mean_grid)

        self.cdict = {}
        path = TEMPLATE_PATH_LATLON.format(
            path_to_climate_dir=self.paths["path-to-climate-dir"] + self.setup["climate_path_to_latlon_file"] + "/")
        self.climate_data_interpolator = monica_run_lib.create_climate_geoGrid_interpolator_from_json_file(path,
                                                                                                      self.wgs84_crs,
                                                                                                      self.soil_crs,
                                                                                                      self.cdict)
        print("created climate_data to gk5 interpolator: ", path)

        with open(f"{self.path_to_out}/spot_setup.out", "a") as _:
            _.write(f"{datetime.now()} grids load producer\n\n")

        with open(self.setup["sim.json"]) as _:
            sim_json = json.load(_)
        # change start and end date according to setup
        if self.setup["start_date"]:
            sim_json["climate.csv-options"]["start-date"] = str(self.setup["start_date"])
        if self.setup["end_date"]:
            sim_json["climate.csv-options"]["end-date"] = str(self.setup["end_date"])
            # sim_json["include-file-base-path"] = paths["include-file-base-path"]

        # read template site.json
        with open(self.setup["site.json"]) as _:
            site_json = json.load(_)

        scenario = self.setup["scenario"]
        if len(scenario) > 0 and scenario[:3].lower() == "rcp":
            site_json["EnvironmentParameters"]["rcp"] = scenario

        # read template crop.json
        with open(self.setup["crop.json"]) as _:
            crop_json = json.load(_)

        crop_json["CropParameters"]["__enable_vernalisation_factor_fix__"] = self.setup[
            "use_vernalisation_fix"] if "use_vernalisation_fix" in self.setup else False

        # create environment template from json templates
        self.env_template = monica_io3.create_env_json_from_json_config({
            "crop": crop_json,
            "site": site_json,
            "sim": sim_json,
            "climate": ""
        })
        if self.shared_id:
            self.env_template["sharedId"] = self.shared_id
        self.env_template["csvViaHeaderOptions"] = sim_json["climate.csv-options"]
        self.env_template["params"]["userCropParameters"]["__enable_T_response_leaf_expansion__"] = self.setup[
            "LeafExtensionModifier"]
        self.env_template["params"]["simulationParameters"]["UseNMinMineralFertilisingMethod"] = self.setup[
            "fertilization"]
        self.env_template["params"]["simulationParameters"]["UseAutomaticIrrigation"] = self.setup["irrigation"]
        self.env_template["params"]["simulationParameters"]["NitrogenResponseOn"] = self.setup["NitrogenResponseOn"]
        self.env_template["params"]["simulationParameters"]["WaterDeficitResponseOn"] = self.setup[
            "WaterDeficitResponseOn"]
        self.env_template["params"]["simulationParameters"]["EmergenceMoistureControlOn"] = self.setup[
            "EmergenceMoistureControlOn"]
        self.env_template["params"]["simulationParameters"]["EmergenceFloodingControlOn"] = self.setup[
            "EmergenceFloodingControlOn"]

        self.orig_crop_params = copy.deepcopy(self.env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"])


    def run_producer(self, params):
        sent_env_count = 0

        start_setup_time = time.perf_counter()

        gcm = self.setup["gcm"]
        rcm = self.setup["rcm"]
        scenario = self.setup["scenario"]
        ensmem = self.setup["ensmem"]
        version = self.setup["version"]
        crop_id = self.setup["crop-id"]
        region_name = self.setup["region_name"]

        if region_name and len(region_name) > 0:
            # Create the soil mask for the specific region
            path_to_soil_grid_ow = self.paths["path-to-data-dir"] + DATA_GRID_SOIL
            mask = self.create_mask_from_shapefile(NUTS1_REGIONS, region_name, path_to_soil_grid_ow)

            # Apply the soil mask to the soil grid
            soil_grid_copy = self.soil_grid.copy()
            self.soil_grid[mask == False] = -8888
            self.soil_grid[soil_grid_copy == -9999] = -9999

        # set current calibration parameter values
        ps = copy.deepcopy(self.orig_crop_params)
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
        self.env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"] = ps


        # unknown_soil_ids = set()
        soil_id_cache = {}
        for srow, scol in self.calib_row_cols:
            #print("srow:", srow, "scol:", scol)

            soil_id = int(self.soil_grid[srow, scol])
            if soil_id == self.nodata_value:
                continue

            # get coordinate of clostest climate element of real soil-cell
            sh = self.yllcorner + (self.scellsize / 2) + (self.srows - srow - 1) * self.scellsize
            sr = self.xllcorner + (self.scellsize / 2) + scol * self.scellsize
            crow, ccol = self.climate_data_interpolator(sr, sh)

            tcoords = {}
            def get_coords(crs, sr, sh):
                if crs not in tcoords:
                    tcoords[crs] = self.soil_crs_to_x_transformers[crs].transform(sr, sh)
                return tcoords[crs]

            if soil_id in soil_id_cache:
                soil_profile = soil_id_cache[soil_id]
            else:
                soil_profile_group = soil_io3.get_soil_profile_group(self.soil_db_con, soil_id)
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
                        with open(self.path_to_prod_out_file, "a") as _:
                            _.write(f"no most_layers for soil_profile with id {soil_id}\n")
                        continue
                else:
                    with open(self.path_to_prod_out_file, "a") as _:
                        _.write(f"soil_profile for id {soil_id} has no valid layers\n")
                    continue
                soil_id_cache[soil_id] = soil_profile
            if not soil_profile or len(soil_profile) == 0:
                with open(self.path_to_prod_out_file, "a") as _:
                    _.write(f"soil_profile for id {soil_id} not valid\n")
                continue
            else:
                # print("soil:", soil_profile)
                self.env_template["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

            # check if current grid cell is used for agriculture
            if self.setup["landcover"]:
                lur, luh = get_coords(self.landuse_crs, sr, sh)
                landuse_id = self.landuse_interpolate(lur, luh)
                if landuse_id not in [2, 3, 4]:
                    continue

            demr, demh = get_coords(self.dem_crs, sr, sh)
            height_nn = self.dem_interpolate(demr, demh)

            slr, slh = get_coords(self.slope_crs, sr, sh)
            slope = self.slope_interpolate(slr, slh)

            gw_min_r, gw_min_h = get_coords(self.gw_min_crs, sr, sh)
            min_groundwater_depth = self.gw_min_interpolate(gw_min_r, gw_min_h)

            gw_max_r, gw_max_h = get_coords(self.gw_max_crs, sr, sh)
            max_groundwater_depth = self.gw_max_interpolate(gw_max_r, gw_max_h)

            gw_mean_r, gw_mean_h = get_coords(self.gw_mean_crs, sr, sh)
            mean_groundwater_depth = self.gw_mean_interpolate(gw_mean_r, gw_mean_h)

            # irr_r, irr_h = get_coords(self.irrigation_crs, sr, sh)
            # irrigation = int(self.irrigation_interpolate(irr_r, irr_h))

            # setting groundwater level
            # if self.setup["groundwater-level"]:
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

            if self.setup["groundwater-level"] == "MINMAX":
                # Assign min and max groundwater depths to the environment template
                self.env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
                self.env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [
                    min_groundwater_depth / 100, "m"]
                self.env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [
                    max_groundwater_depth / 100, "m"]
            # elif setup["groundwater-level"] == "MEAN":
            #     # Assign mean groundwater depth to the environment template
            #     env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
            #     env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [
            #         mean_groundwater_depth / 100, "m"]
            #     env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [
            #         mean_groundwater_depth / 100, "m"]
            elif self.setup["groundwater-level"] == "MIN":
                # Assign min and max groundwater depths to the environment template
                self.env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
                self.env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [
                    min_groundwater_depth / 100, "m"]
                self.env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [
                    min_groundwater_depth / 100, "m"]
            elif self.setup["groundwater-level"] == "MAX":
                # Assign min and max groundwater depths to the environment template
                self.env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
                self.env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [
                    max_groundwater_depth / 100, "m"]
                self.env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [
                    max_groundwater_depth / 100, "m"]

            # setting impenetrable layer
            if self.setup["impenetrable-layer"]:
                impenetrable_layer_depth = monica_run_lib.get_value(
                    self.env_template["params"]["userEnvironmentParameters"]["LeachingDepth"])
                layer_depth = 0
                for layer in soil_profile:
                    if layer.get("is_impenetrable", False):
                        impenetrable_layer_depth = layer_depth
                        # print("setting leaching depth of soil_id:", str(soil_id), "to", impenetrable_layer_depth, "m")
                        break
                    layer_depth += monica_run_lib.get_value(layer["Thickness"])
                self.env_template["params"]["userEnvironmentParameters"]["LeachingDepth"] = [impenetrable_layer_depth,
                                                                                        "m"]
                self.env_template["params"]["siteParameters"]["ImpenetrableLayerDepth"] = [impenetrable_layer_depth, "m"]

            if self.setup["elevation"]:
                self.env_template["params"]["siteParameters"]["heightNN"] = float(height_nn)

            if self.setup["slope"]:
                self.env_template["params"]["siteParameters"]["slope"] = slope / 100.0

            if self.setup["latitude"]:
                clat, _ = self.cdict[(crow, ccol)]
                self.env_template["params"]["siteParameters"]["Latitude"] = clat

            if self.setup["CO2"]:
                self.env_template["params"]["userEnvironmentParameters"]["AtmosphericCO2"] = float(self.setup["CO2"])

            if self.setup["O3"]:
                self.env_template["params"]["userEnvironmentParameters"]["AtmosphericO3"] = float(self.setup["O3"])

            if self.setup["FieldConditionModifier"]:
                self.env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"]["species"][
                    "FieldConditionModifier"] = float(self.setup["FieldConditionModifier"])

            if self.setup["StageTemperatureSum"]:
                stage_ts = self.setup["StageTemperatureSum"].split('_')
                stage_ts = [int(temp_sum) for temp_sum in stage_ts]
                orig_stage_ts = \
                    self.env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"]["cultivar"]["="][
                    "StageTemperatureSum"][0]
                if len(stage_ts) != len(orig_stage_ts):
                    stage_ts = orig_stage_ts
                    print('The provided StageTemperatureSum array is not '
                          'sufficiently long. Falling back to original StageTemperatureSum')

                self.env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"]["cultivar"]["="][
                    "StageTemperatureSum"][0] = stage_ts

            subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(gcm=gcm, rcm=rcm, scenario=scenario, ensmem=ensmem,
                                                              version=version, crow=str(crow), ccol=str(ccol))
            for _ in range(4):
                subpath_to_csv = subpath_to_csv.replace("//", "/")
            self.env_template["pathToClimateCSV"] = [
                self.paths["monica-path-to-climate-dir"] + self.setup["climate_path_to_csvs"] + "/" + subpath_to_csv]
            if self.setup["incl_hist"]:
                hist_subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(gcm=gcm, rcm=rcm, scenario="historical",
                                                                       ensmem=ensmem, version=version,
                                                                       crow=str(crow), ccol=str(ccol))
                for _ in range(4):
                    hist_subpath_to_csv = hist_subpath_to_csv.replace("//", "/")
                self.env_template["pathToClimateCSV"].insert(0, self.paths["monica-path-to-climate-dir"] + self.setup[
                    "climate_path_to_csvs"] + "/" + hist_subpath_to_csv)
            #print("pathToClimateCSV:", self.env_template["pathToClimateCSV"])

            self.env_template["customId"] = {
                "setup_id": self.setup_id,
                "srow": srow, "scol": scol,
                "crow": int(crow), "ccol": int(ccol),
                "soil_id": soil_id,
                "env_id": sent_env_count,
                "nodata": False,
                "shared-id": self.shared_id,
            }

            sent_env_count += 1

            self.prod_socket.send_json(self.env_template)
            print("sent env ", sent_env_count, " customId: ", self.env_template["customId"])

        with open(f"{self.path_to_out}/spot_setup.out", "a") as _:
            _.write(f"{datetime.now()} Sending final last message (producer) \n")
        # send a last message will be just forwarded by monica to signify last
        if self.env_template:
            self.env_template["pathToClimateCSV"] = ""
            self.env_template["customId"] = {
                "no_of_sent_envs": sent_env_count,
                "nodata": True
            }
            self.prod_socket.send_json(self.env_template)

            # print("crows/cols:", crows_cols)
        # cs__.close()
        stop_setup_time = time.perf_counter()
        #print("\nSetup ", self.setup_id, ":", sent_env_count, " envs took ", (stop_setup_time - start_setup_time), " seconds")

    def run_consumer(self):
        year_to_biomasses = defaultdict(list)
        envs_received = 0
        no_of_envs_expected = None
        while True:
            try:
                msg: dict = self.cons_socket.recv_json()  # encoding="latin-1"

                custom_id = msg["customId"]
                if "no_of_sent_envs" in custom_id:
                    no_of_envs_expected = custom_id["no_of_sent_envs"]
                else:
                    envs_received += 1
                    for data in msg.get("data", []):
                        results = data.get("results", [])
                        for vals in results:
                            if "Year" in vals and "exportedCutBiomass" in vals:
                                year_to_biomasses[int(vals["Year"])].append(vals["exportedCutBiomass"])

                if no_of_envs_expected == envs_received:  # and writer:
                    with open(self.path_to_cons_out_file, "a") as _:
                        _.write(f"{datetime.now()} last expected env received\n")
                    # print("last expected env received")
                    return year_to_biomasses

            except zmq.error.Again as _e:
                with open(self.path_to_cons_out_file, "a") as _:
                    _.write(f"no response from the server (with {self.cons_socket.RCVTIMEO} ms timeout)\n")
                print(f"no response from the server (with timeout={self.cons_socket.RCVTIMEO} ms)")
                break
        return None


    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # vector = MaxAssimilationRate, AssimilateReallocation, RootPenetrationRate
        params = dict(zip(vector.name, vector))
        self.run_producer(params)
        year_to_biomasses = self.run_consumer()
        sim_list = flatten_biomasses_dict(year_to_biomasses)
        print("len(sim_list):", len(sim_list), "== len(self.obs_list):", len(self.obs_flat_list), flush=True)
        with open(self.path_to_prod_out_file, "a") as _:
            _.write(f"{datetime.now()}  len(sim_list): {len(sim_list)} == len(self.obs_list): {len(self.obs_flat_list)}\n")
        # besides the order the length of observation results and simulation results should be the same
        assert len(sim_list) == len(self.obs_flat_list)
        return sim_list if len(sim_list) > 0 else None

    def evaluation(self):
        return self.obs_flat_list

    def objectivefunction(self, simulation, evaluation):
        #return unbiased_rmse_RB(evaluation, simulation)
        return spotpy.objectivefunctions.rmse(evaluation, simulation)
        #return calculate_percentage_difference_new(evaluation, simulation)
        #return calculate_weighted_rmse(evaluation, simulation, self.weight_per_region)

def calculate_weighted_rmse(evaluation, simulation, weight_per_region):
    """
    Calculate the weighted RMSE (Root Mean Squared Error).

    .. math::

        RMSE_weighted = \\sqrt{\\frac{\\sum_{i=1}^N w_i * (sim_i - obs_i)^2}{\\sum_{i=1}^N w_i}}

        w_i = \\frac{pixels_i}{\\sum pixels}

    :param evaluation: Observed data to compare with simulation data.
    :type evaluation: list or numpy array of numeric values

    :param simulation: Simulated data to compare with evaluation data.
    :type simulation: list or numpy array of numeric values

    :param pixels: Number of pixels in each region.
    :type pixels: list or numpy array of numeric values

    :return: Weighted RMSE
    :rtype: float
    """


    if len(evaluation) == len(simulation):
        obs = np.array(evaluation)
        sim = np.array(simulation)
        #pixels_array = np.array(pixels) It's used for weights calculation and they're already stored in the csv file
        
        # Calculate weights
        #weights = pixels_array / np.sum(pixels_array)
        
        # Weighted squared differences
        weighted_squared_diff = weight_per_region * (sim - obs) ** 2
        
        # Weighted RMSE calculation
        weighted_rmse = np.sqrt(np.nansum(weighted_squared_diff) / 1)
        
        return weighted_rmse
    else:
        logging.warning("evaluation, simulation, and pixels must have the same non-zero length.")
        return np.nan



def calculate_percentage_difference_new(evaluation: list, simulation: list) -> float:
    """
    Calculate the mean absolute percentage difference between observed (evaluation) and simulated values.

    :param evaluation: Observed data to compare with simulation data.
    :type evaluation: list of numeric values

    :param simulation: Simulation data to compare with evaluation data.
    :type simulation: list of numeric values

    :return: Mean absolute percentage difference, or np.nan if division by zero occurs or
             if the lengths do not match.
    :rtype: float
    """
    if len(evaluation) == len(simulation) > 0:
        percentage_differences = []
        for eval_val, sim_val in zip(evaluation, simulation):
            if eval_val != 0:
                percentage_difference = abs((eval_val - sim_val) / eval_val * 100)
                percentage_differences.append(percentage_difference)
            else:
                percentage_differences.append(np.nan)  # Handle division by zero
        
        # Convert list to numpy array to handle np.nan and compute the mean
        percentage_differences = np.array(percentage_differences)
        
        # Return the mean, ignoring nan values
        return np.nanmean(percentage_differences)
    else:
        logging.warning("evaluation and simulation lists do not have the same length or are empty.")
        return np.nan



def rmse(evaluation, simulation):
    """
    Root Mean Squared Error

        .. math::

         RMSE=\\sqrt{\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})^2}

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Root Mean Squared Error
    :rtype: float
    """
    if len(evaluation) == len(simulation) > 0:
        return np.sqrt(mse(evaluation, simulation))
    else:
        logging.warning("evaluation and simulation lists do not have the same length.")
        return np.nan



def mse(evaluation, simulation):
    """
    Mean Squared Error

        .. math::

         MSE=\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})^2

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Mean Squared Error
    :rtype: float
    """

    if len(evaluation) == len(simulation):
        obs, sim = np.array(evaluation), np.array(simulation)
        mse = np.nanmean((obs - sim) ** 2)
        return mse
    else:
        #logging.warning(
        #    "evaluation and simulation lists does not have the same length."
        #)
        return np.nan


def calculate_mbe(evaluation, simulation):
    """
    Calculate the Mean Bias Error.

    Args:
    y_true (array-like): True values.
    y_pred (array-like): Predicted values.

    Returns:
    float: Mean Bias Error.
    """
    if len(evaluation) == len(simulation):
        obs, sim = np.array(evaluation), np.array(simulation)
        mbe = np.nanmean(sim-obs)#swapped obs-sim to sim-obs
    return mbe


def rmse(evaluation, simulation):
    """
    Root Mean Squared Error

        .. math::

         RMSE=\\sqrt{\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})^2}

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Root Mean Squared Error
    :rtype: float
    """
    if len(evaluation) == len(simulation) > 0:
        return np.sqrt(mse(evaluation, simulation))
    else:
        logging.warning("evaluation and simulation lists do not have the same length.")
        return np.nan

#RB addition

def unbiased_rmse_RB(evaluation, simulation):
    """
    Calculate RMSE with prior bias correction.

    Args:
    y_true (array-like): True values.
    y_pred (array-like): Predicted values.

    Returns:
    float: RMSE after bias correction.
    """

    mbe = calculate_mbe(evaluation, simulation)
    print("      mean bias estimate        =",mbe)
    #subtract mean bias estimate from simulation. 
    #This will give bias corrected simulation.
    y_pred_adjusted = simulation - mbe
    return rmse(evaluation, y_pred_adjusted)



def mse(evaluation, simulation):
    """
    Mean Squared Error

        .. math::

         MSE=\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})^2

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Mean Squared Error
    :rtype: float
    """

    if len(evaluation) == len(simulation):
        obs, sim = np.array(evaluation), np.array(simulation)
        mse = np.nanmean((obs - sim) ** 2)
        return mse
    else:
        #logging.warning(
        #    "evaluation and simulation lists does not have the same length."
        #)
        return np.nan







