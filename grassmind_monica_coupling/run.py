#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg-mohnicke@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

import asyncio
import capnp
from datetime import date
import fileinput
import io
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess as sp
import sys
import time
import uuid
from zalfmas_common import common
import zalfmas_capnp_schemas
from zalfmas_services.crop import monica_crop_service
from zalfmas_common.model import monica_io
from zalfmas_common.soil import soil_io
from zalfmas_fbp.run import channels as chans, ports as fbp_ports
capnp_path = Path(os.path.dirname(zalfmas_capnp_schemas.__file__))
sys.path.append(str(capnp_path))
import climate_capnp
import common_capnp
import fbp_capnp
import model_capnp
import crop_capnp
sys.path.append(str(capnp_path / "model" / "monica"))
import monica_management_capnp as mgmt_capnp
import monica_params_capnp
import monica_state_capnp

standalone_config_mbm_lin = {
    "task_id_offset": "0",
    "row": "230",#220,c454
    "col": "423", #"403",
    "rcp": "85",
    "start_year": "2021",
    "end_year": "2023",
    "path_to_channel": "/home/berg/GitHub/monica/_cmake_debug/common/channel",
    "path_to_daily_monica_fbp_component": "/home/berg/GitHub/monica/_cmake_debug/daily-monica-fbp-component",
    "path_to_monica_parameters_dir": "/home/berg/GitHub/monica-parameters",
    "path_to_formind_exe": "/home/berg/GitHub/grassmind_zalf/_cmake_debug/formind",
    "path_to_lat_lon_soil_json": "/home/berg/Desktop/valeh/GRASSMIND/latlon_to_rowcol_with_soilmap.json",
    "path_to_full_weather_file": "/home/berg/Desktop/valeh/weatherData/{row:03}/daily_mean_RES1_C{col:03}R{row:03}.csv",
    "path_to_grassmind_weather_file": "/home/berg/Desktop/valeh/GRASSMIND/4Zalf_10102024_rcp{rcp}/formind_parameters/Climate/daily_mean_RES1_C{col:03}R{row:03}.csv_Grassmind.txt",
    "path_to_grassmind_soil_file": "/home/berg/Desktop/valeh/GRASSMIND/4Zalf_10102024_rcp{rcp}/formind_parameters/Soil/soil_R{row:03}C{col:03}.txt",
    "path_to_grassmind_param_file": "/home/berg/Desktop/valeh/GRASSMIND/4Zalf_10102024_rcp{rcp}/formind_parameters/parameter_R{row:03}C{col:03}I41.par",
    "path_to_result_div": "/home/berg/Desktop/valeh/GRASSMIND/4Zalf_10102024_rcp{rcp}/results/parameter_R{row:03}C{col:03}I41.div",
    "path_to_result_bt": "/home/berg/Desktop/valeh/GRASSMIND/4Zalf_10102024_rcp{rcp}/results/parameter_R{row:03}C{col:03}I41.bt",
    "path_to_biomass_output_file": "biomass_outputs/biomass_rcp{rcp}_R{row}C{col}.csv",
}
standalone_config_rpm_hpc = {
    "task_id_offset": "0",
    "row": "220",
    "col": "454", #"403",
    "rcp": "26",
    "start_year": "2021",
    "end_year": "2098",
    "path_to_channel": "/beegfs/common/data/grassmind/monica/_cmake_release/common/channel",
    "path_to_daily_monica_fbp_component": "/beegfs/common/data/grassmind/monica/_cmake_release/daily-monica-fbp-component",
    "path_to_monica_parameters_dir": "/beegfs/common/data/grassmind/monica-parameters",
    "path_to_formind_exe": "/beegfs/common/data/grassmind/grassmind_zalf/src/formind",
    "path_to_lat_lon_soil_json": "/beegfs/common/data/grassmind/data/latlon_to_rowcol_with_soilmap.json",
    "path_to_full_weather_file": "/beegfs/common/data/grassmind/data/weather_data/rcp_{rcp}/{row:03}/daily_mean_RES1_C{col:03}R{row:03}.csv",
    "path_to_grassmind_weather_file": "/beegfs/common/data/grassmind/data/formind_parameters/Climate/daily_mean_RES1_C{col:03}R{row:03}.csv_Grassmind.txt",
    "path_to_grassmind_soil_file": "/beegfs/common/data/grassmind/data/formind_parameters/Soil/soil_R{row:03}C{col:03}.txt",
    "path_to_grassmind_param_file": "/beegfs/common/data/grassmind/data/formind_parameters/parameter_R{row:03}C{col:03}I41.par",
    "path_to_result_div": "/beegfs/common/data/grassmind/data/results/parameter_R{row:03}C{col:03}I41.div",
    "path_to_result_bt": "/beegfs/common/data/grassmind/data/results/parameter_R{row:03}C{col:03}I41.bt",
    "path_to_biomass_output_file": "biomass_outputs/biomass_rcp{rcp}_R{row}C{col}.csv",
}
standalone_config_mbm_win = {
    "task_id_offset": "0",
    "row": "220",
    "col": "403",
    "rcp": "26",
    "start_year": "2021",
    "end_year": "2098",
    "path_to_channel": "C:/Users/berg/development/monica_win64_3.6.36.daily_fbp_component/bin/channel.exe",
    "path_to_daily_monica_fbp_component": "C:/Users/berg/development/monica_win64_3.6.36.daily_fbp_component/bin/daily-monica-fbp-component.exe",
    "path_to_monica_parameters_dir": "C:/Users/berg/development/monica_win64_3.6.36.daily_fbp_component/monica-parameters",
    "path_to_formind_exe": "C:/Users/berg/Desktop/valeh/4Zalf_10102024_rcp{rcp}/formind.exe",
    "path_to_lat_lon_soil_json": "C:/Users/berg/Desktop/valeh/latlon_to_rowcol_with_soilmap.json",
    "path_to_full_weather_file": "C:/Users/berg/Desktop/valeh/weatherData/{row:03}/daily_mean_RES1_C{col:03}R{row:03}.csv",
    "path_to_grassmind_weather_file": "C:/Users/berg/Desktop/valeh/4Zalf_10102024_rcp{rcp}/formind_parameters\Climate/daily_mean_RES1_C{col:03}R{row:03}.csv_Grassmind.txt",
    "path_to_grassmind_soil_file": "C:/Users/berg/Desktop/valeh/4Zalf_10102024_rcp{rcp}/formind_parameters/Soil/soil_R{row:03}C{col:03}.txt",
    "path_to_grassmind_param_file": "C:/Users/berg/Desktop/valeh/4Zalf_10102024_rcp{rcp}/formind_parameters/parameter_R{row:03}C{col:03}I41.par",
    "path_to_result_div": "C:/Users/berg/Desktop/valeh/4Zalf_10102024_rcp{rcp}/results/parameter_R{row:03}C{col:03}I41.div",
    "path_to_result_bt": "C:/Users/berg/Desktop/valeh/4Zalf_10102024_rcp{rcp}/results/parameter_R{row:03}C{col:03}I41.bt",
    "path_to_biomass_output_file": "biomass_outputs/biomass_rcp{rcp}_R{row}C{col}.csv",
}
standalone_config_vk_win = {
    "task_id_offset": "0",
    "row": "220",
    "col": "403",
    "rcp": "26",
    "start_year": "2021",
    "end_year": "2098",
    "path_to_channel": "C:/Users/khaledi/development/monica_win64_3.6.37/bin/channel.exe",
    "path_to_daily_monica_fbp_component": "C:/Users/khaledi/development/monica_win64_3.6.36.daily_fbp_component/bin/daily-monica-fbp-component.exe",
    "path_to_monica_parameters_dir": "C:/Users/khaledi/development/monica_win64_3.6.36.daily_fbp_component/monica-parameters",
    "path_to_formind_exe": "E:/4Zalf_10102024_rcp{rcp}/formind.exe",
    "path_to_lat_lon_soil_json": "E:/4Zalf_10102024_rcp{rcp}/latlon_to_rowcol_with_soilmap.json",
    "path_to_full_weather_file": "E:/4Zalf_10102024_rcp26/weatherData/{row:03}/daily_mean_RES1_C{col:03}R{row:03}.csv",
    "path_to_grassmind_weather_file": "E:/4Zalf_10102024_rcp{rcp}/formind_parameters/Climate/daily_mean_RES1_C{col:03}R{row:03}.csv_Grassmind.txt",
    "path_to_grassmind_soil_file": "E:/4Zalf_10102024_rcp{rcp}/formind_parameters/Soil/soil_R{row:03}C{col:03}.txt",
    "path_to_grassmind_param_file": "E:/4Zalf_10102024_rcp{rcp}/formind_parameters/parameter_R{row:03}C{col:03}I41.par",
    "path_to_result_div": "E:/4Zalf_10102024_rcp{rcp}/results/parameter_R{row:03}C{col:03}I41.div",
    "path_to_result_bt": "E:/4Zalf_10102024_rcp{rcp}/results/parameter_R{row:03}C{col:03}I41.bt",
    "path_to_biomass_output_file": "biomass_outputs/biomass_rcp{rcp}_R{row}C{col}.csv",
}
async def main(config: dict):
    common.update_config(config, sys.argv, print_config=True, allow_new_keys=False)

    con_man = common.ConnectionManager()
    channels = []
    procs = []

    slurm_array_job_id = os.getenv("SLURM_ARRAY_JOB_ID", None)
    slurm_task_id = os.getenv("SLURM_ARRAY_TASK_ID", None)
    task_id_offset = int(config["task_id_offset"])
    if slurm_task_id:
        # iterate the weather file folder
        one_param_file = config["path_to_grassmind_param_file"].format(row=config["row"], col=config["col"], rcp=config["rcp"])
        params_dir = os.path.dirname(one_param_file)
        all_param_files = []
        for entry in os.listdir(params_dir):
            if entry.startswith("parameter_R") and entry.endswith("I41.par"):
                all_param_files.append(os.path.join(params_dir, entry))
        all_param_files.sort()
        selected_file = all_param_files[int(slurm_task_id) + task_id_offset - 1]
        row = int(selected_file[-14:-11])
        col = int(selected_file[-10:-7])
    else:
        row = int(config["row"])
        col = int(config["col"])

    paths = {
        "full_weather": config["path_to_full_weather_file"].format(row=row, col=col, rcp=config["rcp"]),
        "weather": config["path_to_grassmind_weather_file"].format(row=row, col=col, rcp=config["rcp"]),
        "soil": config["path_to_grassmind_soil_file"].format(row=row, col=col, rcp=config["rcp"]),
        "params": config["path_to_grassmind_param_file"].format(row=row, col=col, rcp=config["rcp"]),
        "formind": config["path_to_formind_exe"].format(row=row, col=col, rcp=config["rcp"]),
        "div": config["path_to_result_div"].format(row=row, col=col, rcp=config["rcp"]),
        "bt": config["path_to_result_bt"].format(row=row, col=col, rcp=config["rcp"]),
        "shm": f"/dev/shm/{uuid.uuid4()}/",
        "biomass_out": config["path_to_biomass_output_file"].format(job_id=slurm_array_job_id,
                                                                    rcp=config["rcp"], row=row, col=col),
        "lat_lon_soil": config["path_to_lat_lon_soil_json"],
    }

    # copy grassmind files into ramdisk
    if os.path.exists("/dev/shm/"):
        os.makedirs(os.path.join(paths["shm"], "grassmind"))

        # copy formind executable to speed up execution
        shutil.copy(paths["formind"], os.path.join(paths["shm"], "grassmind"))
        paths["formind"] = os.path.join(paths["shm"], "grassmind", "formind")

        # copy lat/lon to row/col to soil_id file
        shutil.copy(paths["lat_lon_soil"], paths["shm"])
        paths["lat_lon_soil"] = os.path.join(paths["shm"], os.path.basename(paths["lat_lon_soil"]))

        params_dir = os.path.dirname(paths["params"])
        # copy pin file
        shutil.copy(os.path.join(params_dir, "init41.pin"), os.path.join(paths["shm"], "grassmind"))

        # create and copy observation folder
        os.makedirs(os.path.join(paths["shm"], "grassmind", "Observation"))
        shutil.copytree(os.path.join(params_dir, "Observation"), os.path.join(paths["shm"], "grassmind", "Observation"), dirs_exist_ok=True)

        # create and copy management folder
        os.makedirs(os.path.join(paths["shm"], "grassmind", "Management"))
        shutil.copytree(os.path.join(params_dir, "Management"), os.path.join(paths["shm"], "grassmind", "Management"), dirs_exist_ok=True)

        # weather data file will be created from full weather data
        os.makedirs(os.path.join(paths["shm"], "grassmind", "Climate"))
        #shutil.copy(paths["weather"], os.path.join(paths["shm"], "grassmind", "Climate"))
        # update path to weather file to be created
        paths["weather"] = os.path.join(paths["shm"], "grassmind", "Climate", os.path.basename(paths["weather"]))

        # soil data file will be created from MONICA soil
        os.makedirs(os.path.join(paths["shm"], "grassmind", "Soil"))
        #shutil.copy(paths["soil"], os.path.join(paths["shm"], "grassmind", "Soil"))
        # update path to the soil file to be created
        paths["soil"] = os.path.join(paths["shm"], "grassmind", "Soil", os.path.basename(paths["soil"]))

        # copy params file
        shutil.copy(paths["params"], os.path.join(paths["shm"], "grassmind"))
        # update path to params file
        paths["params"] = os.path.join(paths["shm"], "grassmind", os.path.basename(paths["params"]))

        # create results dir
        os.makedirs(os.path.join(paths["shm"], "results"), exist_ok=True)
        # update path to the two grasssmind output files
        paths["div"] = os.path.join(paths["shm"], "results", os.path.basename(paths["div"]))
        paths["bt"] = os.path.join(paths["shm"], "results", os.path.basename(paths["bt"]))

        # create biomass output folder
        os.makedirs(os.path.join(paths["shm"], "biomass_out"), exist_ok=True)
        # update path to biomass file to be created
        paths["biomass_out"] = os.path.join(paths["shm"], "biomass_out", os.path.basename(paths["biomass_out"]))

    # create biomass output folder (in case of linux and shm already happend)
    if not os.path.exists(os.path.dirname(paths["biomass_out"])):
        os.makedirs(os.path.dirname(paths["biomass_out"]), exist_ok=True)

    # update the parameter file to run just for a single day
    with fileinput.input(paths["params"], inplace=True) as f:
        for line_no, line in enumerate(f):
            if line_no == 8:
                print("float	TimeEnd		0.00274")
            else:
                print(line, end="")

    # create soil profile
    soil_db_con = sqlite3.connect("../data/germany/buek200.sqlite")
    soil_profile = None
    soil_id = None
    soil_id_exists = False
    with open(paths["lat_lon_soil"], "r") as f:
        for lat_lon, row_col, soil_id_ in json.load(f):
            if row_col == [row, col]:
                soil_id = soil_id_[0]
                soil_id_exists = True
                #soil_profile = soil_io.soil_parameters(soil_db_con, soil_id)
                soil_profile_group = soil_io.get_soil_profile_group(soil_db_con, soil_id)
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
                        break
                break
    if not soil_profile:
        print("no soil profile found for row/col:", row, "/", col, "soil_id:", soil_id if soil_id_exists else "-")
        exit(0)

    try:
        first_chan, first_reader_sr, first_writer_sr = chans.start_first_channel(config["path_to_channel"])
        channels.append(first_chan)
        first_reader = await con_man.try_connect(first_reader_sr, cast_as=fbp_capnp.Channel.Reader)

        # create the three channels for the three ports
        channels.append(chans.start_channel(config["path_to_channel"],
                                            "env_in", first_writer_sr, name="env"))
        channels.append(chans.start_channel(config["path_to_channel"],
                                            "events_in", first_writer_sr, name="events"))
        channels.append(chans.start_channel(config["path_to_channel"],
                                            "result_out", first_writer_sr, name="result"))
        channels.append(chans.start_channel(config["path_to_channel"],
                                            "serialized_state_in", first_writer_sr, name="serialized_state"))
        channels.append(chans.start_channel(config["path_to_channel"],
                                            "serialized_state_out", first_writer_sr, name="serialized_state"))
        channels.append(chans.start_channel(config["path_to_channel"],
                                            "port_infos", first_writer_sr, name="port_infos"))#,
                                            #port=9991,
                                            #reader_srts="r_in"))

        port_srs = {"in": {}, "out": {}}
        port_infos_reader_sr = None
        port_infos_writer = None
        port_infos_msg = fbp_capnp.PortInfos.new_message()
        in_ports = []
        out_ports = []
        for i in range(len(channels)-1):
            p = (await first_reader.read()).value.as_struct(common_capnp.Pair)
            c_id = p.fst.as_text()
            info = p.snd.as_struct(fbp_capnp.Channel.StartupInfo)
            print("channel:", c_id, "reader_sr:", info.readerSRs[0], "writer_sr:", info.writerSRs[0])
            if c_id[-3:] == "_in":
                port_name = c_id[:-3]
                in_ports.append({"name": port_name, "sr": info.readerSRs[0]})
                port_srs["in"][port_name] = info.writerSRs[0]
            elif c_id[-4:] == "_out":
                port_name = c_id[:-4]
                out_ports.append({"name": port_name, "sr": info.writerSRs[0]})
                port_srs["out"][port_name] = info.readerSRs[0]
            elif c_id == "port_infos":
                port_infos_writer = await con_man.try_connect(info.writerSRs[0], cast_as=fbp_capnp.Channel.Writer)
                port_infos_reader_sr = info.readerSRs[0]
            else:
                print("Error received unknown channel startupInfo with id:", c_id)
                exit(1)
        port_infos_msg.inPorts = in_ports
        port_infos_msg.outPorts = out_ports

        # run monica daily fbp component
        procs.append(sp.Popen([
            config["path_to_daily_monica_fbp_component"],
            #"--verbose",
            port_infos_reader_sr
        ], env={"MONICA_PARAMETERS": config["path_to_monica_parameters_dir"],
                "systemroot": os.getenv("systemroot", "")},
        stdout=sp.DEVNULL))

        # write the config to the config channel
        await port_infos_writer.write(value=port_infos_msg)

        with open("sim.json") as _:
            sim_json = json.load(_)
        with open("site.json") as _:
            site_json = json.load(_)
            site_json["SiteParameters"]["SoilProfileParameters"] = soil_profile
        with open("crop.json") as _:
            crop_json = json.load(_)
        env_template = monica_io.create_env_json_from_json_config({
            "crop": crop_json,
            "site": site_json,
            "sim": sim_json,
            "climate": ""
        })

        env_writer = await con_man.try_connect(port_srs["in"]["env"], cast_as=fbp_capnp.Channel.Writer)
        env = common_capnp.StructuredText.new_message(value=json.dumps(env_template),
                                                      structure={"json": None})
        await env_writer.write(value=fbp_capnp.IP.new_message(content=env))
        print("send env on env channel")

        state_reader = await con_man.try_connect(port_srs["out"]["serialized_state"], cast_as=fbp_capnp.Channel.Reader)
        state_writer = await con_man.try_connect(port_srs["in"]["serialized_state"], cast_as=fbp_capnp.Channel.Writer)
        event_writer = await con_man.try_connect(port_srs["in"]["events"], cast_as=fbp_capnp.Channel.Writer)
        output_reader = await con_man.try_connect(port_srs["out"]["result"], cast_as=fbp_capnp.Channel.Reader)

        start_year = int(config["start_year"])
        end_year = int(config["end_year"])
        iso_dates = []
        grassmind_climate = ["rain[mm]\tTemperature[degC]\tRadiation[mmolm-2s-1]\tDaylength[h]\tPET[mm]\tCO2[ppm]\n"]
        monica_climate = []
        with open(paths["full_weather"]) as _:
            header = _.readline()[:-1].split(",")
            h2i = {h.replace('"',''): i for i, h in enumerate(header)}
            #_.readline() # skip units
            for line in _.readlines():
                data = line.split(",")
                data[0] = data[0].replace('"','')
                #iso_date = data[h2i["iso-date"]]
                iso_date = data[h2i["date"]]

                if int(iso_date[:4]) < start_year:
                    continue
                if int(iso_date[:4]) > end_year:
                    continue

                iso_dates.append(iso_date)
                grassmind_climate.append(f'{float(data[h2i["precip"]])}\t{float(data[h2i["tavg"]])}\t{53.1 * float(data[h2i["globrad"]])}\t{float(data[h2i["DayLength"]])}\t{float(data[h2i["PET"]])}\t{400}')
                monica_climate.append(mgmt_capnp.Params.DailyWeather.new_message(data=[
                    {"key": "tavg", "value": float(data[h2i["tavg"]])},
                    {"key": "tmin", "value": float(data[h2i["tmin"]])},
                    {"key": "tmax", "value": float(data[h2i["tmax"]])},
                    {"key": "wind", "value": float(data[h2i["wind"]])},
                    {"key": "globrad", "value": float(data[h2i["globrad"]])},
                    {"key": "precip", "value": float(data[h2i["precip"]])},
                    {"key": "relhumid", "value": float(data[h2i["relhumid"]])}
                ]))

        ferts = {
            "AN": {
                "id": "AN",
                "name": "Ammonium Nitrate",
                "carbamid": 0,
                "nh4": 0.5,
                "no3": 0.5,
            }
        }

        # create monica crop
        abs_events = {
            "2021-03-01": create_sowing_event(monica_crop_service.Crop(
                {"id": "Grass_Species4", "name": "Grass Species 4"}, "../data/params/crops/species/Grass_Species4.json",
                {"id": "Grass_CLV4", "name": "Grass CLV 4"}, "../data/params/crops/cultivars/Grass_CLV4.json",
                "../data/params/crops/residues/grass-ley.json")),
        }
        rel_events = {
            #"03-01": create_sowing_event(monica_crop_service.Crop(
            #    {"id": "Grass_Species4", "name": "Grass Species 4"}, "../data/params/crops/species/Grass_Species4.json",
            #    {"id": "Grass_CLV4", "name": "Grass CLV 4"}, "../data/params/crops/cultivars/Grass_CLV4.json",
            #    "../data/params/crops/residues/grass-ley.json")),
            "06-15": create_cutting_event([
                {"organ": "leaf", "value": 0.15, "unit": "lai", "cutOrLeft": "left", "exportPercentage": 100.0},
                {"organ": "shoot", "value": 100, "unit": "biomass", "cutOrLeft": "left", "exportPercentage": 100.0}]),
            "06-20": create_n_demand_fert_event(n_demand=20.0, depth=0.3, partition=ferts["AN"]),
            "08-15": create_cutting_event([
                {"organ": "leaf", "value": 0.4, "unit": "lai", "cutOrLeft": "left", "exportPercentage": 100.0},
                {"organ": "shoot", "value": 100, "unit": "biomass", "cutOrLeft": "left", "exportPercentage": 100.0}]),
            "08-20": create_n_demand_fert_event(n_demand=20.0, depth=0.3, partition=ferts["AN"]),
            "10-15": create_cutting_event(
                [{"organ": "leaf", "value": 0.4, "unit": "lai", "cutOrLeft": "left", "exportPercentage": 100.0},
                 {"organ": "shoot", "value": 100, "unit": "biomass", "cutOrLeft": "left",
                  "exportPercentage": 100.0}]),
            "10-20": create_n_demand_fert_event(n_demand=20.0, depth=0.3, partition=ferts["AN"]),
        }

        #print("wrote openBracket on event channel")
        prev_year = None
        for day_index, iso_date in enumerate(iso_dates[:-365]):
            current_date = date.fromisoformat(iso_date)
            if prev_year != current_date.year:
                #print(current_date.year, end=" ", flush=True)
                prev_year = current_date.year

            await event_writer.write(value=fbp_capnp.IP.new_message(type="openBracket"))

            if iso_date in abs_events:
                current_events = abs_events[iso_date] if abs_events[iso_date] is list else [abs_events[iso_date]]
                for cev in current_events:
                    event = cev(current_date)
                    await event_writer.write(value=fbp_capnp.IP.new_message(content=event))

            rel_date = iso_date[5:]
            if rel_date in rel_events:
                current_events = rel_events[rel_date] if rel_events[rel_date] is list else [rel_events[rel_date]]
                for cev in current_events:
                    event = cev(current_date)
                    await event_writer.write(value=fbp_capnp.IP.new_message(content=event))

            # finally send weather to do daily step
            weather_event = create_weather_event(monica_climate[day_index])(current_date)
            await event_writer.write(value=fbp_capnp.IP.new_message(content=weather_event))
            #print("send weather event for day:", iso_date, "on event channel")

            # save state
            save_state_event = create_save_state_event(10, False)(current_date)
            #print("send saveState event for day:", iso_date, "on event channel")
            await event_writer.write(value=fbp_capnp.IP.new_message(content=save_state_event))

            # read state
            state_msg = await state_reader.read()
            grassmind_total_biomass_kg_per_ha = None
            if state_msg.which() == "value":
                #print("received state for day:", iso_date, "on state channel")
                state_ip = state_msg.value.as_struct(fbp_capnp.IP)
                old_state = state_ip.content.as_struct(monica_state_capnp.RuntimeState)
                new_state, grassmind_total_biomass_kg_per_ha = run_grassmind_on_monica_state(old_state, day_index, grassmind_climate, paths)
                await state_writer.write(value=fbp_capnp.IP.new_message(content=new_state))
            else:
                print("received done on state channel")

            await event_writer.write(value=fbp_capnp.IP.new_message(type="closeBracket"))
            #print("send closeBracket on event channel")

            # reading output
            out_msg = await output_reader.read()
            mo_biomass = None
            if out_msg.which() == "value":
               out_ip = out_msg.value.as_struct(fbp_capnp.IP)
               st = out_ip.content.as_struct(common_capnp.StructuredText)
               r = json.loads(st.value)
               mo_res = r["data"][0]["results"][0]
               mo_biomass = mo_res["AbBiom"]
               #print(r)
            else:
               print("received done on output channel")

            #print(iso_date, "biomass gm:", grassmind_total_biomass_kg_per_ha, "mo:", mo_biomass)
            with open(paths["biomass_out"], "a") as f:
                f.write(f"{iso_date},{mo_biomass}\n")

        await event_writer.close()
        await env_writer.close()
        await state_writer.close()
        time.sleep(3)

        for channel in channels:
            channel.terminate()
        print(f"{os.path.basename(__file__)}: all channels terminated")
        for proc in procs:
            proc.terminate()

    except Exception as e:
        print(f"exception terminated {os.path.basename(__file__)} early. Exception:", e)

        #for process in process_id_to_process.values():
        #    process.terminate()
        for channel in channels:
            channel.terminate()
        for proc in procs:
            proc.terminate()

    final_biomass_out_file = config["path_to_biomass_output_file"].format(job_id=slurm_array_job_id,
                                                                          rcp=config["rcp"], row=row, col=col)
    if not os.path.exists(os.path.dirname(final_biomass_out_file)):
        os.makedirs(os.path.dirname(final_biomass_out_file), exist_ok=True)
    shutil.copy(paths["biomass_out"], final_biomass_out_file)
    shutil.rmtree(paths["shm"], ignore_errors=True)


def run_grassmind_on_monica_state(old_state, day_index, grassmind_climate, paths):
    if not old_state.modelState._has("currentCropModule"):
        return old_state, None

    new_state = old_state
    if not hasattr(run_grassmind_on_monica_state, "initial_param_values"):
        run_grassmind_on_monica_state.initial_param_values = {
            "SpecificLeafArea": list(old_state.modelState.currentCropModule.cultivarParams.specificLeafArea),
            "StageKcFactor": list(old_state.modelState.currentCropModule.cultivarParams.stageKcFactor),
            "DroughtStressThreshold": list(old_state.modelState.currentCropModule.cultivarParams.droughtStressThreshold),
            "CropSpecificMaxRootingDepth": float(old_state.modelState.currentCropModule.cultivarParams.cropSpecificMaxRootingDepth),
        }

    with open(paths["weather"], "wt") as f:
        f.write(grassmind_climate[0])
        for line in grassmind_climate[day_index+1:day_index+1 + 365]:
            f.write(line)
            f.write("\n")

    with open(paths["soil"], "wt") as f:
        f.write(create_grassmind_soil_from_state(old_state))

    p = sp.Popen([paths["formind"], paths["params"]], stdout=sp.DEVNULL)
    p.wait()     #grassmind run

    # read .div file to get the current fractions
    rel_species_abundance = None
    with open(paths["div"]) as f:
        lines = f.readlines()
        rel_species_abundance = list(map(float, lines[4].split("\t")[2:6]))
    with open(paths["bt"]) as f:
        lines = f.readlines()
        total_biomass_kg_per_ha = float(lines[4].split("\t")[1])*100000000.0 # t -> kg ha-1

    if rel_species_abundance:
        params = calc_community_level_params(rel_species_abundance)
        #print(old_state.modelState.currentStepDate, "community params:", params)
        new_state = old_state.as_builder()
        cps = new_state.modelState.currentCropModule.cultivarParams
        cps.specificLeafArea = list(map(lambda v: v * params["SpecificLeafArea"],
                                        run_grassmind_on_monica_state.initial_param_values["SpecificLeafArea"]))
        #print("new specificLeafArea:", cps.specificLeafArea)
        cps.stageKcFactor = list(map(lambda v: v * params["StageKcFactor"],
                                     run_grassmind_on_monica_state.initial_param_values["StageKcFactor"]))
        #print("new stageKcFactor:", cps.stageKcFactor)
        cps.droughtStressThreshold = list(map(lambda v: v * params["DroughtStressThreshold"],
                                              run_grassmind_on_monica_state.initial_param_values[
                                                  "DroughtStressThreshold"]))
        #print("new droughtStressThreshold:", cps.droughtStressThreshold)
        cps.cropSpecificMaxRootingDepth = params["CropSpecificMaxRootingDepth"]
        #print("new cropSpecificMaxRootingDepth:", cps.cropSpecificMaxRootingDepth)
    return new_state, total_biomass_kg_per_ha

def create_grassmind_soil_from_state(state):
    sc = state.modelState.soilColumn
    sb = io.StringIO()
    sb.write("Silt\tClay\tSand\n")
    silt = 1 - sc.layers[0].sps.soilSandContent - sc.layers[0].sps.soilClayContent
    sb.write(f"{silt}\t{sc.layers[0].sps.soilClayContent}\t{sc.layers[0].sps.soilSandContent}\n")
    sb.write("\n")
    sb.write("Layer\tRWC[-]\tFC[V%]\tPWP[V%]\tMinN[gm-2]\tPOR[V%]\tKS[mm/d]\n")
    for i, l in enumerate(sc.layers):
        n_kg_per_m3 = l.soilNO3 + l.soilNO2 + l.soilNH4
        n_kg_per_m2 = n_kg_per_m3 * 0.01 # m3 -> m2 (0.1m layer thickness)
        n_g_per_m2 = n_kg_per_m2 * 1000.0 # kg -> g
        ks = l.sps._get("lambda") * 1000.0 # m/d -> mm/d
        sb.write(f"{i}\t{round(l.soilMoistureM3, 5)}\t{round(l.sps.fieldCapacity*100.0, 2)}\t{round(l.sps.permanentWiltingPoint*100.0,2)}\t{round(n_g_per_m2, 5)}\t{round(l.sps.saturation*100.0, 2)}\t{round(ks,2)}\n")
    return sb.getvalue()

def calc_community_level_params(s_i_s):
    p_i_s = {
        "SpecificLeafArea": {"S_PFTs": [0.36, 0.47, 0.03, 0.36], "min": 0.6, "max": 1.4, "values": [0.8, 0.75, 0.7, 0.6]},
        "StageKcFactor": {"S_PFTs": [0.36, 0.16, 0.03, 0.4], "min": 0.6, "max": 1.4, "values": [1.3, 1.1, 0.91, 0.83]},
        "DroughtStressThreshold": {"S_PFTs": [0.0, 0.16, 0.78, 0.58], "min": 0.2, "max": 1.0, "values": [0.35, 0.5, 0.75, 0.9]},
        "CropSpecificMaxRootingDepth": {"S_PFTs": [0.36, 0.31, 0.03, 0.98], "min": 0.1, "max": 0.3, "values": [0.15, 0.2, 0.25, 0.3]},
    }
    res = {}
    for param_name, ps in p_i_s.items():
        sum_sip_x_si_x_pi = 0
        sum_sip_x_si = 0
        for i, s_i_p in enumerate(ps["S_PFTs"]):
            sum_sip_x_si_x_pi += s_i_p * s_i_s[i] * ps["values"][i]
            sum_sip_x_si += s_i_p * s_i_s[i]
        res[param_name] = ps["min"] + (sum_sip_x_si_x_pi/sum_sip_x_si)*(ps["max"] - ps["min"])
    return res

def create_save_state_event(no_of_prev_days_to_serialize=10, serialize_as_json=False):
    return lambda at: mgmt_capnp.Event.new_message(type="saveState",
                                                   at={"date": {"year": at.year, "month": at.month, "day": at.day}},
                                                   info={"id": at.isoformat()},
                                                   params=mgmt_capnp.Params.SaveState.new_message(
                                                       noOfPreviousDaysSerializedClimateData=no_of_prev_days_to_serialize,
                                                       asJson=serialize_as_json))

def create_weather_event(daily_weather):
    return lambda at: mgmt_capnp.Event.new_message(type="weather",
                                                   at={"date": {"year": at.year, "month": at.month, "day": at.day}},
                                                   info={"id": at.isoformat()},
                                                   params=daily_weather)

def create_sowing_event(crop):
    return lambda at: mgmt_capnp.Event.new_message(type="sowing",
                                                   at={"date": {"year": at.year, "month": at.month, "day": at.day}},
                                                   info={"id": at.isoformat()},
                                                   params=mgmt_capnp.Params.Sowing.new_message(
                                                       cultivar="Grass_CLV4",
                                                       crop=crop))

def create_harvest_event():
    return lambda at: mgmt_capnp.Event.new_message(type="harvest",
                                                   at={"date": {"year": at.year, "month": at.month, "day": at.day}},
                                                   info={"id": at.isoformat()},
                                                   params=mgmt_capnp.Params.Harvest.new_message())

def create_n_fert_event(amount: float, partition: dict):
    return lambda at: mgmt_capnp.Event.new_message(type="mineralFertilization",
                                                   at={"date": {"year": at.year, "month": at.month, "day": at.day}},
                                                   info={"id": at.isoformat()},
                                                   params=mgmt_capnp.Params.MineralFertilization.new_message(
                                                       amount = amount,
                                                       partition = partition))

def create_n_demand_fert_event(n_demand: float, depth: float, partition: dict, stage: int = 1):
    return lambda at: mgmt_capnp.Event.new_message(type="nDemandFertilization",
                                                   at={"date": {"year": at.year, "month": at.month, "day": at.day}},
                                                   info={"id": at.isoformat()},
                                                   params=mgmt_capnp.Params.NDemandFertilization.new_message(
                                                       nDemand = n_demand,
                                                       depth = depth,
                                                       stage = stage,
                                                       partition = partition))

def create_cutting_event(cutting_spec: list[dict]):
    return lambda at: mgmt_capnp.Event.new_message(type="cutting",
                                                   at={"date": {"year": at.year, "month": at.month, "day": at.day}},
                                                   info={"id": at.isoformat()},
                                                   params=mgmt_capnp.Params.Cutting.new_message(cuttingSpec = cutting_spec))


if __name__ == '__main__':
    asyncio.run(capnp.run(main(standalone_config_rpm_hpc)))
    #asyncio.run(capnp.run(main(standalone_config_mbm_lin)))
    #asyncio.run(capnp.run(main(standalone_config_vk_win)))
