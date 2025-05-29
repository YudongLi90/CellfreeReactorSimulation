# Chemical species
"""
Author: Yudong Li (Yudong.Li@nrel.gov)
"""
import numpy as np
from .constants import *

Composition = np.zeros((TOTAL_SPECIES, 6))
Elements = ["C", "H", "N", "O", "P", "S"]
MWs = [12, 1, 13, 16, 31, 32]
Composition[Glc] = [6, 12, 0, 6, 0, 0] #C6H12O6
Composition[ATP] = [10, 16, 5, 13, 3, 0] #C10H16N5O13P3
Composition[ADP] = [10, 15, 5, 1, 2, 0] #C10H15N5O10P2
Composition[G6P] = [6, 13, 0, 9, 1, 0] # C6H13O9P
Composition[F6P] = [6, 13, 0, 9, 1, 0] #C6H13O9P
Composition[FBP] = [6, 14, 0, 12, 2, 0] # C6H14O12P2
Composition[DHAP] = [3, 7, 0, 6, 1, 0] # C3H7O6P
Composition[G3P] = [3, 7, 0, 6, 1, 0] # C3H7O6P
Composition[BPG13] = [3, 8, 0, 10, 2, 0] #C3H8O10P2
Composition[NAD] = [21, 27, 7, 14, 2, 0] #C21H27N7O14P2 # Nicotinamide adenine dinucleotide
Composition[NADH] = [21, 27, 7, 14, 2, 0] #C21H27N7O14P2 (https://www.nad.com/nad-vs-nadh, should be one more H- than NAD+)
Composition[P_I] = [0, 0, 0, 4, 1, 0] # PO4^3-
Composition[NADP] = [21, 29, 7, 17, 3, 0] #C21H29N7O17P3
Composition[NADPH] = [21, 29, 7, 17, 3, 0] #C21H29N7O17P3
Composition[PG3] = [3, 7, 0, 7, 1, 0] #C3H7O7P
Composition[PG2] = [3, 7, 0, 7, 1, 0] #C3H7O7P
Composition[PEP] = [3, 7, 0, 7, 1, 0] #C3H7O7P
Composition[PYR] = [3, 4, 0, 3, 0, 0] #C3H4O3
Composition[CoA] = [21, 36, 7, 16, 3, 1] #C21H36N7O16P3S
Composition[AcCoA] = [23, 38, 7, 17, 3, 1] #C23H38N7O17P3S
Composition[CO2] = [1, 0, 0, 2, 0, 0]
Composition[HMGCoA] = [27, 44, 7, 20, 3, 1] #C27H44N7O20P3S
Composition[AcAcCoA] = [25, 40, 7, 18, 3, 1] #C25H40N7O18P3S
Composition[MEV] = [6, 12, 0, 4, 0, 0] #C6H12O4(Mevalonic Acid)
Composition[MEVP] = [6, 13, 0, 7, 1, 0] #C6H13O7P Phosphomevalonic acid
Composition[MEVPP] = [6, 13, 0, 10, 2, 0]#C6H13O10P2
Composition[IPP] = [5, 12, 0, 7, 2, 0] #C5H12O7P2
Composition[DMAPP] = [5, 12, 0, 7, 2, 0] #C5H12O7P2
Composition[GPP] = [10, 20, 0, 7, 2, 0] #C10H20O7P2
Composition[PPi] = [0, 0, 0, 7, 2, 0] #P2O7^4- # Inorganic pyrophosphate
Composition[LIM] = [10, 16, 0, 0, 0, 0] #C10H16




def total_carbon(c):
    total_carbon = 0
    for i, specie in enumerate(c):
        total_carbon += specie*Composition[i][0]

    return total_carbon

def total_P(c):
    total_P = 0
    for i, specie in enumerate(c):
        total_P += specie*Composition[i][4]

    return total_P

def molecular_weight(c):
    assert len(c) == 6
    res = 0
    for n, w in zip (c, MWs):
        res += n*w
    return res

def mass_balance(reactants, products):
    w_r = 0
    w_p = 0
    elements_r = []
    elements_p = []
    for react in reactants:
        
        w_r += molecular_weight(Composition[react])
    
    for prod in products:
        w_p += molecular_weight(Composition[prod])

    print(f"Reactant MW: {w_r}, Product MW: {w_p}, difference (react-prod): {w_r - w_p}.")
    return w_r-w_p




    

