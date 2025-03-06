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

import subprocess as sp
import sys

local_run = True

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
        "mode": "mbm-local-remote",
        "prod-port": prod_port if prod_port else "6666",  # local: 6667, remote 6666
        "cons-port": cons_port if cons_port else "7777",  # local: 6667, remote 6666
        "server": server if server else "login01.cluster.zalf.de",
        "setups-file": "sim_setups_calibration_VK.csv",
        "path_to_out": "out/",
        "run-setups": "[1-3]",
        "path_to_python": "python" if local_run else "/home/rpm/.conda/envs/clim4cast/bin/python",
        "repetitions": "1000",
    }

    update_config(config, sys.argv, print_config=True, allow_new_keys=False)

    rs_ranges = config["run-setups"][1:-1].split(",")
    run_setups = []
    for rsr in rs_ranges:
        rs_r = rsr.split("-")
        if 1 < len(rs_r) <= 2:
            run_setups.extend(range(int(rs_r[0]), int(rs_r[1]) + 1))
        elif len(rs_r) == 1:
            run_setups.append(int(rs_r[0]))
    # run_setups = json.loads(config["run-setups"])

    procs = []
    for setup_id in run_setups:
        procs.append(sp.Popen([
            config["path_to_python"],
            "run-calibration2.py",
            f"mode={config['mode']}",
            f"server={config['server']}",
            f"prod-port={config['prod-port']}",
            f"cons-port={config['prod-port']}",
            f"setups-file={config['setups-file']}",
            f"run-setups=[{setup_id}]",
            f"path_to_out={config['path_to_out']}/calib_{setup_id}/",
            f"repetitions={config['repetitions']}",
        ]))

    for p in procs:
        p.wait()
        p.terminate()

    print("run-parallel-calibrations")

if __name__ == "__main__":
    run_parallel_calibrations()
