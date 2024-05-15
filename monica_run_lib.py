#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Susanne Schulz <susanne.schulz@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

import csv
import json
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from pyproj import Transformer
from datetime import date, timedelta


def read_csv(path_to_setups_csv, key="run-id"):
    "read sim setup from csv file"
    with open(path_to_setups_csv) as _:
        key_to_data = {}
        # determine seperator char
        dialect = csv.Sniffer().sniff(_.read(), delimiters=';,\t')
        _.seek(0)
        # read csv with seperator char
        reader = csv.reader(_, dialect)
        header_cols = next(reader)
        for row in reader:
            data = {}
            for i, header_col in enumerate(header_cols):
                value = row[i]
                if value.lower() in ["true", "false"]:
                    value = value.lower() == "true"
                if header_col == key:
                    value = int(value)
                data[header_col] = value
            key_to_data[int(data[key])] = data
        return key_to_data


def read_sim_setups(path_to_setups_csv):
    "read sim setup from csv file"
    with open(path_to_setups_csv) as setup_file:
        setups = {}
        # determine seperator char
        dialect = csv.Sniffer().sniff(setup_file.read(), delimiters=';,\t')
        setup_file.seek(0)
        # read csv with seperator char
        reader = csv.reader(setup_file, dialect)
        header_cols = next(reader)
        for row in reader:
            data = {}
            for i, header_col in enumerate(header_cols):
                value = row[i]
                if value.lower() in ["true", "false"]:
                    value = value.lower() == "true"
                if i == 0:
                    value = int(value)
                data[header_col] = value
            setups[int(data["run-id"])] = data
        return setups


def read_header(path_to_ascii_grid_file):
    "read metadata from esri ascii grid file"
    metadata = {}
    header_str = ""
    with open(path_to_ascii_grid_file) as _:
        for i in range(0, 6):
            line = _.readline()
            header_str += line
            sline = [x for x in line.split() if len(x) > 0]
            if len(sline) > 1:
                metadata[sline[0].strip().lower()] = float(sline[1].strip())
    return metadata, header_str


def create_ascii_grid_interpolator(grid, meta_data, ignore_nodata=True):
    "read an ascii grid into a map, without the no-data values"
    "grid - 2D array of values"

    rows, cols = grid.shape

    cellsize = int(meta_data["cellsize"])
    xll = int(meta_data["xllcorner"])
    yll = int(meta_data["yllcorner"])
    nodata_value = meta_data["nodata_value"]

    xll_center = xll + cellsize // 2
    yll_center = yll + cellsize // 2
    yul_center = yll_center + (rows - 1)*cellsize

    points = []
    values = []

    for row in range(rows):
        for col in range(cols):
            value = grid[row, col]
            if ignore_nodata and value == nodata_value:
                continue
            r = xll_center + col * cellsize
            h = yul_center - row * cellsize
            points.append([r, h])
            values.append(value)

    return NearestNDInterpolator(np.array(points), np.array(values))


def get_value(list_or_value):
   return list_or_value[0] if isinstance(list_or_value, list) else list_or_value


def create_seed_harvest_geoGrid_interpolator_and_read_data(path_to_csv_file, worldGeodeticSys84, geoTargetGrid, ilr_seed_harvest_data):
    "read seed/harvest dates and apoint climate stations"

    wintercrop = {
        "WW": True,
        "SW": False,
        "WR": True,
        "WRa": True,
        "WB": True,
        "SM": False,
        "GM": False,
        "SBee": False,
        "SU": False,
        "SB": False,
        "SWR": True, 
        "CLALF": False,
        "PO": False,
        "RYE": True
    }

    with open(path_to_csv_file) as _:
        reader = csv.reader(_)

        #print "reading:", path_to_csv_file

        # skip header line
        next(reader)

        points = [] # climate station position (lat, long transformed to a geoTargetGrid, e.g gk5)
        values = [] # climate station ids

        transformer = Transformer.from_crs(worldGeodeticSys84, geoTargetGrid, always_xy=True) 

        prev_cs = None
        prev_lat_lon = [None, None]
        #data_at_cs = defaultdict()
        for row in reader:
            
            # first column, climate station
            cs = int(row[0])

            # if new climate station, store the data of the old climate station
            if prev_cs is not None and cs != prev_cs:

                llat, llon = prev_lat_lon
                #r_geoTargetGrid, h_geoTargetGrid = transform(worldGeodeticSys84, geoTargetGrid, llon, llat)
                r_geoTargetGrid, h_geoTargetGrid = transformer.transform(llon, llat)
                    
                points.append([r_geoTargetGrid, h_geoTargetGrid])
                values.append(prev_cs)

            crop_id = row[3]
            is_wintercrop = wintercrop[crop_id]
            ilr_seed_harvest_data[crop_id]["is-winter-crop"] = is_wintercrop

            base_date = date(2001, 1, 1)

            sdoy = int(float(row[4]))
            ilr_seed_harvest_data[crop_id]["data"][cs]["sowing-doy"] = sdoy
            sd = base_date + timedelta(days = sdoy - 1)
            ilr_seed_harvest_data[crop_id]["data"][cs]["sowing-date"] = "0000-{:02d}-{:02d}".format(sd.month, sd.day)

            esdoy = int(float(row[8]))
            ilr_seed_harvest_data[crop_id]["data"][cs]["earliest-sowing-doy"] = esdoy
            esd = base_date + timedelta(days = esdoy - 1)
            ilr_seed_harvest_data[crop_id]["data"][cs]["earliest-sowing-date"] = "0000-{:02d}-{:02d}".format(esd.month, esd.day)

            lsdoy = int(float(row[9]))
            ilr_seed_harvest_data[crop_id]["data"][cs]["latest-sowing-doy"] = lsdoy
            lsd = base_date + timedelta(days = lsdoy - 1)
            ilr_seed_harvest_data[crop_id]["data"][cs]["latest-sowing-date"] = "0000-{:02d}-{:02d}".format(lsd.month, lsd.day)

            digit = 1 if is_wintercrop else 0
            if crop_id == 'CLALF': digit = 2
            
            hdoy = int(float(row[6]))
            ilr_seed_harvest_data[crop_id]["data"][cs]["harvest-doy"] = hdoy
            hd = base_date + timedelta(days = hdoy - 1)
            ilr_seed_harvest_data[crop_id]["data"][cs]["harvest-date"] = "000{}-{:02d}-{:02d}".format(digit, hd.month, hd.day)

            ehdoy = int(float(row[10]))
            ilr_seed_harvest_data[crop_id]["data"][cs]["earliest-harvest-doy"] = ehdoy
            ehd = base_date + timedelta(days = ehdoy - 1)
            ilr_seed_harvest_data[crop_id]["data"][cs]["earliest-harvest-date"] = "000{}-{:02d}-{:02d}".format(digit, ehd.month, ehd.day)

            lhdoy = int(float(row[11]))
            ilr_seed_harvest_data[crop_id]["data"][cs]["latest-harvest-doy"] = lhdoy
            lhd = base_date + timedelta(days = lhdoy - 1)
            ilr_seed_harvest_data[crop_id]["data"][cs]["latest-harvest-date"] = "000{}-{:02d}-{:02d}".format(digit, lhd.month, lhd.day)

            lat = float(row[1])
            lon = float(row[2])
            prev_lat_lon = (lat, lon)      
            prev_cs = cs

        ilr_seed_harvest_data[crop_id]["interpolate"] = NearestNDInterpolator(np.array(points), np.array(values))


def create_climate_geoGrid_interpolator_from_json_file(path_to_latlon_to_rowcol_file, worldGeodeticSys84, geoTargetGrid, cdict):
    "create interpolator from json list of lat/lon to row/col mappings"
    with open(path_to_latlon_to_rowcol_file) as _:
        points = []
        values = []

        transformer = Transformer.from_crs(worldGeodeticSys84, geoTargetGrid, always_xy=True) 

        for latlon, rowcol in json.load(_):
            row, col = rowcol
            clat, clon = latlon
            try:
                cr_geoTargetGrid, ch_geoTargetGrid = transformer.transform(clon, clat)
                cdict[(row, col)] = (round(clat, 4), round(clon, 4))
                points.append([cr_geoTargetGrid, ch_geoTargetGrid])
                values.append((row, col))
                #print "row:", row, "col:", col, "clat:", clat, "clon:", clon, "h:", h, "r:", r, "val:", values[i]
            except Exception as e:
                print("row/col:", (row,col), "clat/clon:", (clat, clon), "cr/ch:", (cr_geoTargetGrid, ch_geoTargetGrid), "Exception:", e)
                continue

        return NearestNDInterpolator(np.array(points), np.array(values))

