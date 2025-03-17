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

import json
import os
import subprocess as sp
import sys

local_run = False

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


def run_parallel_calibrations(server=None, prod_port=None, cons_port=None):
    config = {
        "mode": "mbm-local-remote",#"mbm-local-local", #"mbm-local-remote",
        "prod-port": prod_port if prod_port else "6666",  # local: 6667, remote 6666
        "cons-port": cons_port if cons_port else "7777",  # local: 6667, remote 6666
        "server": server if server else "login01.cluster.zalf.de",
        "setups-file": "sim_setups_calibration_VK.csv",
        "path_to_out": "out/",
        "run-setup": "1",
        "no-of-parallel-calibrations": "20",
        "sm": "10",
        "rcp": "26",
        "path_to_grassmind_biomass_files": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/rpm/projects/monica/project/ipp_grassland/rcp{rcp}/{sm}/" \
            if local_run else "/project/rcp{rcp}/{sm}/",
        "path_to_python": "python" if local_run else "/home/rpm/.conda/envs/clim4cast/bin/python",
        "repetitions": "10",
    }

    update_config(config, sys.argv, print_config=True, allow_new_keys=False)

    sm = int(config["sm"])
    rcp = config["rcp"]
    setup_id = int(config["run-setup"])
    max_no_of_par_calibs = int(config["no-of-parallel-calibrations"])

    path_to_biomass_files_dir = config["path_to_grassmind_biomass_files"].format(rcp=config["rcp"], sm=config["sm"])
    # iterate files in directory
    procs = []
    files = os.listdir(path_to_biomass_files_dir)
    for i, file in enumerate(files):
        rc = file.split("_")[1].split("I")[0]
        row = int(rc[1:4])
        col = int(rc[5:8])

        procs.append(sp.Popen([
            config["path_to_python"],
            "run-calibration2.py",
            f"mode={config['mode']}",
            f"server={config['server']}",
            f"prod-port={config['prod-port']}",
            f"cons-port={config['cons-port']}",
            f"setups-file={config['setups-file']}",
            f"run-setups=[{config['run-setup']}]",
            f"row={row}",
            f"col={col}",
            f"path_to_out={config['path_to_out']}/rcp{rcp}/{sm}/calib_setup{setup_id}_r{row}_c{col}/",
            f"repetitions={config['repetitions']}",
            f"path_to_grassmind_biomass_file={path_to_biomass_files_dir}/{file}",
        ]))

        if len(procs) == max_no_of_par_calibs or (i + 1) == len(files):
            for p in procs:
                p.wait()
                p.terminate()
            procs = []

    # output best parameters for
    path_to_best_params_csv = f"{config['path_to_out']}/best_params_setup-{setup_id}_sm-{sm}_rcp-{rcp}.json"
    with open(path_to_best_params_csv, "a") as f:
        f.write(f"row, col, sm, SpecificLeafArea, StageKcFactor, DroughtStressThreshold, CropSpecificMaxRootingDepth\n")

    path = f"{config['path_to_out']}/rcp{rcp}/{sm}/"
    for dir_ in os.listdir(path):
        r, c = dir_.split("_")[2:]
        path2 = f"{path}/{dir_}/"
        best = {"sm": sm, "rcp": rcp, "row": int(r[1:]), "col": int(c[1:])}
        try:
            with open(f"{path2}/best.out", "r") as f:
                lines = f.readlines()
                for line in lines[4:8]:
                    k, v = line.split(":")
                    best[k] = float(v)
        except FileNotFoundError:
            continue

        with open(path_to_best_params_csv, "a") as f:
            f.write(f"{best['row']}, {best['col']}, {best['sm']}, {best['SpecificLeafArea']}, {best['StageKcFactor']}, {best['DroughtStressThreshold']}, {best['CropSpecificMaxRootingDepth']}\n")

    print("run-parallel-calibrations")

if __name__ == "__main__":
    run_parallel_calibrations()
