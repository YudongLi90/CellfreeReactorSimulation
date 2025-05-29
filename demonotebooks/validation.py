"""
Author: Yudong Li (Yudong.Li@nrel.gov)
"""
from g2t.core.glucose_reactions_packedbed import get_rxn_rates_pkb 
from g2t.core.constants import *
from scikits.odes.odeint import odeint
from scipy import integrate as igt
from g2t.core.species import *
from g2t import get_rxn_rates, EnzymeSystem

import numpy as np
import matplotlib.pyplot as plt

from g2t.core.enzyme_kinetics import *
import pandas as pd

plt.rc('font', family='sans-serif', serif='Arial', size = 14)
plt.rc('text', usetex=False)
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('axes', labelsize=18)

def plot_validation_curve(save=False, savename="ExpValidation", parameters=None):
    c0 = np.zeros(TOTAL_SPECIES)
    c0[Glc] = 100
    c0[ADP] = 2
    c0[ATP] = 2
    c0[FBP] = 1
    c0[P_I] = 10
    c0[NAD] = 0.5
    c0[NADP] = 1.5
    c0[PPi] = 0.5 # from Thiamine pyrophosphate
    c0[CoA] = 1.5
    c0[BPG13] = 0.25
    enzSys = EnzymeSystem(fitExp=True)
    if parameters is not None:
        for k, v in parameters.items():
            enzSys.set_enzyme_Kcat(k, v)

    def kcat_rhs(t, y, ydot): # this is the experiment
        enzSys.get_rate(y, ydot)

    extra_options = {'old_api': False, 'rtol': 1e-6, 'atol': 1e-12, 'max_steps': 50000}
    solution = odeint(kcat_rhs, np.linspace(0, 60*60*24*6, 1000), c0, method='bdf', **extra_options)

    time = solution.values.t
    c = solution.values.y

    pgi = Pgi(packedBed=False)
    N = len(time)
    carbon = np.zeros(N)
    P = np.zeros(N) 
    rate_pgi = np.zeros(N)
    for i in range(N):
        carbon[i] = total_carbon(c[i, :])
        P[i] = total_P(c[i,:])
        rate_pgi[i] = pgi.get_rate(c[i, G6P], c[i, F6P])

    expdata = {}
    expdata["Experimental Data"] = pd.read_excel("GlcToLim_expData.xlsx", sheet_name=0)

    fig= plt.figure(figsize=(7, 6))
    ax1 = fig.subplots()

    ax1.set_xlabel("Time(days)")
    ax1.set_ylabel("Glucose/Limonene (mM)")

    # experiment
    ax1.plot(time/(3600*24), c[:, Glc],c='k', label="Glucose")
    ax1.plot(time/(3600*24), c[:, LIM],c='b', label="Limonene")

    # plot experimental data fig 4a
    for key,val in expdata.items():
        ax1.plot(val["Time(s)"]/(3600*24),val["Glc(mM)"], 'o', markersize=10, c='k', label="Glu, exp")
        ax1.plot(val["Time(s)"]/(3600*24),val["Lim(mM)"], '^', markersize=10, c='b', label="Lim, exp")

    # for sensitivity param sweep
    tt = [1.75,1.75]
    cc = [0,100]

    ax1.legend(loc = "best")
    fig.tight_layout()
    if save:
        plt.savefig(savename+".png", dpi=300, bbox_inches="tight")


