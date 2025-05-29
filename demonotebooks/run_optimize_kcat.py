"""
Author: Yudong Li (Yudong.Li@nrel.gov)
"""
from g2t.utils.optimizer import Optimizer
from g2t.core.constants import *

import pandas as pd
import numpy as np

if __name__ == "__main__":

    exp_res = pd.read_excel("demonotebooks/GlcToLim_expData.xlsx",sheet_name="Sheet1")

    # initial limonene production
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
    # glutathione 5mM

    optim = Optimizer(exp_res["Time(s)"].to_numpy(), exp_res["Glc(mM)"].to_numpy(), exp_res["Lim(mM)"].to_numpy(), c0)
    pbounds = {}

    # for k, val_lst in optim.enzSys.enzyme_Kcat_range.items():
    #     if k == 'Idi':
    #         pbounds[k] = (val_lst[0], 30000)
    #     else:
    #         pbounds[k] = (val_lst[0], 10000)

    pbounds = {
                "Hex": [2e4, 8e4],
                "NoxE": [651, 4000],
                "PPase": [1000,6000],
                "Gap": [8e4, 21e4],
                "Tpi": [1e4, 18e4],
                "Pgk": [2e4, 27e4],
                "Pyk": [1e4, 37e4],
                }
    optim.optimize(pbounds,steps=200)

