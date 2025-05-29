"""
Author: Yudong Li (Yudong.Li@nrel.gov)
"""
from g2t import EnzymeSystem
from g2t.core.constants import *

from scikits.odes.odeint import odeint

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from g2t.core.species import *
from g2t.core.enzyme_kinetics import *

from sklearn.metrics import mean_squared_error
from bayes_opt import BayesianOptimization, UtilityFunction
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
from bayes_opt.util import load_logs

import traceback
import os
from importlib import reload 
import logging
reload(logging)
logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.INFO, datefmt='%I:%M:%S')

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

class Optimizer:

    def __init__(self, experiment_t:list, experiemnet_glc: list, experiment_lim: list, c0, logfile="logs.log", restart=False):
        self.exp_t = experiment_t
        self.exp_glc = experiemnet_glc
        self.exp_lim = experiment_lim
        self.c0 = c0
        self.enzSys = EnzymeSystem()
        self.datadir = os.path.join(os.getcwd(), "optimization_log")
        mkdir_p(self.datadir)
        self.logger = logging.getLogger()
        log_filename = os.path.join(self.datadir, f"optimization.log")
        file_handler = logging.FileHandler(filename=log_filename, mode='w')
        file_handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s %(levelname)s:%(message)s', '%I:%M:%S')
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
        self.cachefile = os.path.join(self.datadir, logfile+".json")
        self.cachelogger = JSONLogger(path = self.cachefile, reset=restart)
        self.restart = restart
        self.load_cache = False
        self.check_cache()
    
    def check_cache(self):
        if os.path.isfile(self.cachefile):
            if self.restart:
                print(f"Overwriting logg and start fresh: {self.cachefile}")
                self.cachelogger = JSONLogger(path = self.cachefile, reset=self.restart)
                self.load_cache = False
            else:
                print(f"Loading previous progress: {self.cachefile}")
                self.load_cache = True
        else:
            print(f"Starting fresh, logging into file {self.cachefile}")
            self.load_cache = False
        


    def run_simulation(self):
        def rhs(t, y, ydot):
            self.enzSys.get_rate(y, ydot)
        extra_options = {'old_api': False, 'rtol': 1e-6, 'atol': 1e-12, 'max_steps': 500000}
        try:
            solution = odeint(rhs, np.linspace(0, 60*60*24*6, int(60*60*24*6/5)+1), self.c0, method='bdf', **extra_options)
        except Exception as e:
            self.logger.error(traceback.format_exc())

        time = solution.values.t
        c = solution.values.y

        t2 = np.reshape(time, (len(time), -1))
        self.sim_res = np.hstack((t2, c))
        return self.sim_res
    
    def get_predicted_res(self):
        sim_res = self.run_simulation()
        glc = []
        lim = []
        c = sim_res[:, 1:]
        for t in self.exp_t:
            glc.append(c[sim_res[:,0]==t, Glc].item())
            lim.append(c[sim_res[:,0]==t, LIM].item())
        
        return glc, lim

    def optimize(self, pbounds: dict, steps: int = 2000, restart=False):
        self.restart = restart
        self.check_cache()
        
        def minus_rmse_exp_sim(Kcat_dict):
            for key, value in Kcat_dict.items():
                self.enzSys.set_enzyme_Kcat(key, value)
            sim_glc, sim_lim = self.get_predicted_res()
            total_error = mean_squared_error(sim_glc, self.exp_glc) + mean_squared_error(sim_lim, self.exp_lim)
            return -total_error
        
        optimizer = BayesianOptimization(
            f = None,
            pbounds= pbounds,
            verbose= 2,
            random_state= 1,
            allow_duplicate_points = True,
        )
        self.optimizer = optimizer
        if self.load_cache:
            load_logs(self.optimizer, logs=self.cachefile)
        self.optimizer.subscribe(Events.OPTIMIZATION_STEP, self.cachelogger)

        utility = UtilityFunction(kind="ucb", kappa=2.5, xi=0.0)
        next_point_to_probe = optimizer.suggest(utility)
        self.logger.info(f"Optimization step: 0/{steps}")
        self.logger.info(f"Next point to probe is: {next_point_to_probe}")
        target = minus_rmse_exp_sim(next_point_to_probe)
        optimizer.register(
            params = next_point_to_probe,
            target = target,
        )
        for step in range(steps):
            self.logger.info(f"Optimization step: {step+1}/{steps}")
            next_point_to_probe = optimizer.suggest(utility)
            self.probe_point = next_point_to_probe
            self.logger.info(f"Next point to probe is: {next_point_to_probe}")
            target = minus_rmse_exp_sim(next_point_to_probe)
            optimizer.register(
                params = next_point_to_probe,
                target = target,
            )
            max = optimizer.max
            self.logger.info(optimizer.max)
        self.logger.info(f"Final optimization result: {optimizer.max}")
        return optimizer.max



