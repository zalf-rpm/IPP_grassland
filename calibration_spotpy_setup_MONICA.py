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

import capnp
from datetime import datetime
import json
import logging
import os
from pathlib import Path
import numpy as np
import spotpy
import re
import zmq

#PATH_TO_REPO = Path(os.path.realpath(__file__)).parent
#PATH_TO_MAS_INFRASTRUCTURE_REPO = PATH_TO_REPO / "../mas-infrastructure"
#PATH_TO_CAPNP_SCHEMAS = (PATH_TO_MAS_INFRASTRUCTURE_REPO / "capnproto_schemas").resolve()
#abs_imports = [str(PATH_TO_CAPNP_SCHEMAS)]
#fbp_capnp = capnp.load(str(PATH_TO_CAPNP_SCHEMAS / "fbp.capnp"), imports=abs_imports)

def flatten_biomasses_dict(biomasses: dict):
    flattend_dict = []
    for year in sorted(biomasses.keys()):
        flattend_dict.extend(biomasses[year])
    return flattend_dict

class spot_setup(object):
    def __init__(self, user_params, observations, #prod_writer, cons_reader,
                 path_to_out):#, only_nuts3_region_ids, weight_per_region):
        self.user_params = user_params
        self.params = []
        self.observations = observations
        self.obs_flat_list = flatten_biomasses_dict(observations)
        #self.prod_writer = prod_writer
        #self.cons_reader = cons_reader
        self.path_to_out_file = path_to_out + "/spot_setup.out"
        #self.only_nuts3_region_ids = only_nuts3_region_ids
        #self.weight_per_region = weight_per_region
        self.context = zmq.Context()
        self.prod_channel = self.context.socket(zmq.PUSH)
        self.prod_channel.bind("tcp://*:9998")
        self.cons_channel = self.context.socket(zmq.PULL)
        self.cons_channel.bind("tcp://*:9999")

        if not os.path.exists(path_to_out):
            try:
                os.makedirs(path_to_out)
            except OSError:
                print("spot_setup.__init__: Couldn't create dir:", path_to_out, "!")

        with open(self.path_to_out_file, "a") as _:
            _.write(f"observations: {self.observations}\n")
            _.write(f"obs_flat_list: {self.obs_flat_list}\n")

        for par in user_params:
            par_name = par["name"]
            if "array" in par:
                par["name"] = f"{par_name}_{par['array']}"  # spotpy does not allow two parameters to have the same name
                del par["array"]
            if "derive_function" not in par:  # spotpy does not care about derived params
                self.params.append(spotpy.parameter.Uniform(**par))

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # vector = MaxAssimilationRate, AssimilateReallocation, RootPenetrationRate
        msg_content = dict(zip(vector.name, vector))
        #msg_content["only_nuts3_region_ids"] = self.only_nuts3_region_ids
        self.prod_channel.send_json(msg_content)
        #out_ip = fbp_capnp.IP.new_message(content=json.dumps(msg_content))
        #self.prod_writer.write(value=out_ip).wait()
        with open(self.path_to_out_file, "a") as _:
            _.write(f"{datetime.now()} sent params to monica setup: {vector}\n")
        print("sent params to monica setup:", vector, flush=True)

        #msg = self.cons_reader.read().wait()
        ## check for end of data from in port
        #if msg.which() == "done":
        #    return

        year_to_biomasses = self.cons_channel.recv_json()
        #in_ip = msg.value.as_struct(fbp_capnp.IP)
        #s: str = in_ip.content.as_text()
        #nuts3_region_id_and_year_to_avg_yield = json.loads(s)

        #with open(self.path_to_out_file, "a") as _:
        #    _.write(f"{datetime.now()} jsons loaded cal-sp-set-M\n")
        # print("received monica results:", country_id_and_year_to_avg_yield, flush=True)

        # remove all simulation results which are not in the observed list
        sim_list = flatten_biomasses_dict(year_to_biomasses)
        #for d in self.observations:
        #    key = f"{d['id']}|{d['year']}"
        #    if key in nuts3_region_id_and_year_to_avg_yield:
        #        #if np.isnan(d["value"]):
        #        #    sim_list.append(np.nan)
        #        #else:
        #        sim_list.append(nuts3_region_id_and_year_to_avg_yield[key])
        #    else:
        #        sim_list.append(np.nan)

        #with open(self.path_to_out_file, "a") as _:
        #    _.write(f"{datetime.now()} simulation and observation matchedcal-sp-set-M\n\n")

        print("len(sim_list):", len(sim_list), "== len(self.obs_list):", len(self.obs_flat_list), flush=True)
        with open(self.path_to_out_file, "a") as _:
            #_.write(f"received monica results: {country_id_and_year_to_avg_yield}\n")
            _.write(f"{datetime.now()}  len(sim_list): {len(sim_list)} == len(self.obs_list): {len(self.obs_flat_list)}\n")
            #_.write(f"sim_list: {sim_list}\n")
            #_.write(f"obs_list: {self.obs_flat_list}\n")
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







