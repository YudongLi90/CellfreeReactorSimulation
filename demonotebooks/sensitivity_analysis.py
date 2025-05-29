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

colors = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]
linestyles = ["-", "--", "-.", ":", "-", "--", "-."]
markers = ["4", "2", "3", "1", "+", "x", "."]

def get_default_c0():
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
    return c0

enzSys = EnzymeSystem(fitExp=True)
def kcat_rhs(t, y, ydot): # this is the experiment
    enzSys.get_rate(y, ydot)


def run_sa(save=False):
    enzyme_lst = ['Hex', 'Pgi', 'Pfk', 'Fba', 'Tpi', 'Gap', 'mGap', 'Pgk', 'Pgm', 'Eno', 'Pyk', 'Pdh', 'PhaA', 'Hmgs', 'Hmgr', 'Mvk', 'Pmvk', 'Mdc', 'Idi', 'Fpps_S82F', 'LimSyn', 'PPase', 'NoxE']
    list_factors = [1, 0.5, 0.1, 2]
    for ez in enzyme_lst:
        sensitivity_analysis(ez, list_factors, save=save, savename=f"./enzyme_sensitivity/sensitivity_{ez}.png")


def run_cofactor_sa(save=False):
    cofactor_lst = [ADP, ATP, P_I, NAD, NADP]
    prefix_lst = ["ADP", "ATP", "Pi", "NAD", "NADP"]
    list_multiplier = [0, 0.1, 1, 10]


    for cofactor, prefix in zip(cofactor_lst, prefix_lst):
        time, c_glc, c_lim, list_conc = compute_cofactor_sensitivity(cofactor, list_multiplier)
        plot_cofactor_sensitive_analysis(time, c_glc, c_lim, list_multiplier, labe_prefix=prefix, save=save, save_name=f"./cofactor_sensitivity/cofactor_{prefix}.png")


def sensitivity_analysis(enzyme_name, list_U_factor, save=False, savename=None):
    time, c_glc, c_lim, list_U = compute_sensitivity(enzyme_name=enzyme_name, list_U_factor=list_U_factor)
    plot_sensitive_analysis(time=time, c_glc=c_glc, c_lim=c_lim, list_U_mutiplier=list_U_factor, labe_prefix=enzyme_name, save=save, save_name=savename)




def compute_cofactor_sensitivity(cofactor, list_conc_multiplier, c0=None):
    if c0 == None:
        c0 = get_default_c0()

    list_conc = c0[cofactor]*np.array(list_conc_multiplier)

    enzSys = EnzymeSystem(fitExp=True)
    def kcat_rhs(t, y, ydot): # this is the experiment
        enzSys.get_rate(y, ydot)

    c_glc = []
    c_lim = []

    for conc in list_conc:
        c0[cofactor] = conc
        
        extra_options = {'old_api': False, 'rtol': 1e-6, 'atol': 1e-12, 'max_steps': 50000}
        solution = odeint(kcat_rhs, np.linspace(0, 60*60*24*6, 1000), c0, method='bdf', **extra_options)

        time = solution.values.t
        c = solution.values.y
        
        c_glc.append(c[:, Glc])
        c_lim.append(c[:, LIM])
    
    return time, c_glc, c_lim, list_conc



def compute_sensitivity(enzyme_name, list_U_factor, c0 = None):
    if c0 == None:
        c0 = get_default_c0()

    enzSys = EnzymeSystem(fitExp=True)
    def kcat_rhs(t, y, ydot): # this is the experiment
        enzSys.get_rate(y, ydot)

    Kcat = enzSys.enzymes[enzyme_name].Kcat
    list_U = enzSys.enzymes[enzyme_name].U * np.array(list_U_factor)
    list_kcats = list_U_factor

    c_glc = []
    c_lim = []

    for factor in list_U_factor:
        Kcat_hex_new = Kcat*factor
        enzSys.set_enzyme_Kcat(enzyme_name, Kcat_hex_new)
        
        extra_options = {'old_api': False, 'rtol': 1e-6, 'atol': 1e-12, 'max_steps': 50000}
        solution = odeint(kcat_rhs, np.linspace(0, 60*60*24*6, 1000), c0, method='bdf', **extra_options)

        time = solution.values.t
        c = solution.values.y
        
        c_glc.append(c[:, Glc])
        c_lim.append(c[:, LIM])
    
    return time, c_glc, c_lim, list_U

def plot_cofactor_sensitive_analysis(time, c_glc, c_lim, list_conc_multiplier, labe_prefix, color_glc= 'k', color_lim= 'b', save=False, save_name=None):
    fig= plt.figure(figsize=(7, 6))
    ax1 = fig.subplots()
    ax2 = ax1.twinx()

    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Glucose (mM)")
    ax2.set_ylabel("Limonene (mM)")

    styles = [ "dotted", 'dashdot', "solid",  "dashed", (0,(3, 5, 1, 5, 1, 5))]
    marker_list = ['s', 'o', '>', '*', 'D']
    # styles = linestyles
    # marker_list = markers

    for i, multiplier in enumerate(list_conc_multiplier):
        ax1.plot(time/(3600*24), c_glc[i], linestyle= styles[i], marker=marker_list[i], markevery=50+10*i, markersize=8, c=color_glc, label=f"{labe_prefix}:{multiplier:.3g}X")
        ax2.plot(time/(3600*24), c_lim[i], linestyle= styles[i], marker=marker_list[i], markevery=50+10*i, markersize=8, c=color_lim,)

    ax1.set_ylim([0, 105])
    ax2.set_ylim([0, 105])

    ax1.set_xlim([0,6])
    ax2.set_xlim([0,6])

    ax1.spines['left'].set_color(color_glc)
    ax1.tick_params(axis='y', colors=color_glc)
    ax1.yaxis.label.set_color(color_glc)

    ax2.spines['left'].set_color(color_glc)
    ax2.spines['right'].set_color(color_lim)
    ax2.tick_params(axis='y', colors=color_lim)
    ax2.yaxis.label.set_color(color_lim)

    ax1.legend(loc = "best")
    fig.tight_layout()

    if save and save_name:
        plt.savefig(save_name, dpi=300, bbox_inches="tight")


def plot_sensitive_analysis(time, c_glc, c_lim, list_U_mutiplier, labe_prefix, color_glc='k', color_lim='b', save=False, save_name=None):
    fig= plt.figure(figsize=(7, 6))
    ax1 = fig.subplots()
    ax2 = ax1.twinx()

    ax1.set_xlabel("Time (days)")
    ax1.set_ylabel("Glucose (mM)")
    ax2.set_ylabel("Limonene (mM)")


    styles = [ 'solid', "dashdot", "dotted", "dashed", (0,(3, 5, 1, 5, 1, 5))]
    marker_list = ['s', 'o', '>', '*', 'D']

    for i, U_multiplier in enumerate(list_U_mutiplier):
        ax1.plot(time/(3600*24), c_glc[i], linestyle= styles[i], marker=marker_list[i], markevery=50+10*i, markersize=10, c=color_glc, label=f"{labe_prefix}:{U_multiplier:0.3g}X")
        ax2.plot(time/(3600*24), c_lim[i], linestyle= styles[i], marker=marker_list[i], markevery=50+10*i, markersize=10, c=color_lim,)

    ax1.set_ylim([0, 105])
    ax2.set_ylim([0, 105])

    ax1.set_xlim([0,6])
    ax2.set_xlim([0,6])

    ax1.spines['left'].set_color(color_glc)
    ax1.tick_params(axis='y', colors=color_glc)
    ax1.yaxis.label.set_color(color_glc)

    ax2.spines['left'].set_color(color_glc)
    ax2.spines['right'].set_color(color_lim)
    ax2.tick_params(axis='y', colors=color_lim)
    ax2.yaxis.label.set_color(color_lim)

    ax1.legend(loc = "best")
    fig.tight_layout()

    if save and save_name:
        plt.savefig(save_name, dpi=300, bbox_inches="tight")

