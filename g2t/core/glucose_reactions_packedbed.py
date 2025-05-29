"""
Author: Yudong Li (Yudong.Li@nrel.gov)
"""
import numpy as np
from g2t.core.enzyme_kinetics import *
from .constants import *


def get_rxn_rates_pkb(c, rates):
    hex = Hex()
    pgi = Pgi()
    pfk = Pfk()
    fba = Fba()
    tpi = Tpi()
    gap = Gap()
    mgap = mGap()
    noxe = NoxE()
    pgk = Pgk()
    pgm = Pgm()
    eno = Eno()
    pyk = Pyk()
    pdh = Pdh()
    phaA = PhaA()
    hmgs = Hmgs()
    hmgr = Hmgr()
    mvk = Mvk()
    pmvk = Pmvk()
    mdc = Mdc()
    idi = Idi()
    fpps = Fpps_S82F()
    ppase = PPase()
    limSyn = LimSyn()

    lst_en = [hex, pgi, pfk, fba, tpi, gap, mgap, pgk, pgm, eno, pyk, pdh, phaA, hmgs, hmgr, mvk, pmvk, mdc, idi, fpps, limSyn, ppase, noxe]
    l_mg = [0.001,0.03,0.02,0.04,0.004,0.007,0.03,0.01,0.05,0.05,0.01,0.01,0.02,0.03,0.08,0.02,0.04,0.08,0.07,0.02,0.6,0.2,0.05]
    for en, mg in zip(lst_en, l_mg):
        en.set_enzyme_concentration_with_mg(mg, 800)
        assert(abs(en.Kcat*en.c - en.U) < 1)

    rate_Hex = hex.get_rate(c[Glc], c[ATP], c[G6P], c[ADP])
    rate_Pgi = pgi.get_rate(c[G6P], c[F6P])
    rate_PfkB = pfk.get_rate(c[F6P], c[ATP], c[FBP], c[ADP])
    rate_Fba = fba.get_rate(c[FBP], c[DHAP], c[G3P])
    rate_Tpi = tpi.get_rate(c[DHAP], c[G3P])
    rate_Gap = gap.get_rate(c[G3P], c[NAD], c[P_I], c[BPG13], c[NADH])
    rate_mGap = mgap.get_rate(c[G3P], c[NADP], c[P_I], c[BPG13], c[NADPH])
    rate_NoxE = noxe.get_rate(c[NADH])
    rate_Pgk = pgk.get_rate(c[BPG13], c[ADP], c[PG3], c[ATP])
    rate_dPgm = pgm.get_rate(c[PG3], c[PG2])
    rate_Eno = eno.get_rate(c[PG2], c[PEP])
    rate_PykF = pyk.get_rate(c[PEP], c[ADP], c[PYR], c[ATP])
    rate_PDH = pdh.get_rate(c[CoA], c[NAD], c[PYR], c[AcCoA], c[NADH], c[CO2])
    rate_PhaA = phaA.get_rate(c[AcCoA], c[AcAcCoA], c[CoA])
    rate_Hmgs = hmgs.get_rate(c[AcCoA], c[AcAcCoA], c[HMGCoA], c[CoA])
    rate_Hmgr = hmgr.get_rate(c[HMGCoA], c[NADPH], c[MEV], c[NADP], c[CoA])
    rate_Mvk = mvk.get_rate(c[MEV], c[ATP], c[MEVP], c[ADP])
    rate_Pmvk = pmvk.get_rate(c[MEVP], c[ATP], c[MEVPP], c[ADP])
    rate_Mdc = mdc.get_rate(c[MEVPP], c[ATP], c[IPP], c[ADP], c[CO2], c[P_I])
    rate_Idi = idi.get_rate(c[IPP], c[DMAPP])
    rate_FppsS82F = fpps.get_rate(c[IPP], c[DMAPP], c[GPP], c[PPi])
    rate_PPase = ppase.get_rate(c[PPi], c[P_I])
    rate_LimSyn = limSyn.get_rate(c[GPP], c[LIM], c[PPi])

    rates[Glc] = -rate_Hex
    rates[G6P] = rate_Hex - rate_Pgi
    rates[F6P] = rate_Pgi - rate_PfkB
    rates[FBP] = rate_PfkB - rate_Fba
    rates[G3P] = rate_Fba + rate_Tpi - rate_Gap - rate_mGap
    rates[BPG13] = rate_Gap + rate_mGap - rate_Pgk
    rates[PG3] = rate_Pgk - rate_dPgm
    rates[PG2] = rate_dPgm - rate_Eno
    rates[PEP] = rate_Eno - rate_PykF
    rates[PYR] = rate_PykF - rate_PDH
    # rates[AcCoA] = rate_PDH - rate_PhaA - rate_Hmgs
    # rates[AcAcCoA] = 0.5*rate_PhaA - rate_Hmgs
    rates[AcCoA] = rate_PDH - 2*rate_PhaA - rate_Hmgs
    rates[AcAcCoA] = rate_PhaA - rate_Hmgs
    rates[HMGCoA] = rate_Hmgs - rate_Hmgr
    rates[MEV] = rate_Hmgr - rate_Mvk
    rates[MEVP] = rate_Mvk - rate_Pmvk
    rates[MEVPP] = rate_Pmvk - rate_Mdc
    rates[IPP] = rate_Mdc - rate_Idi - rate_FppsS82F
    rates[DMAPP] = rate_Idi - rate_FppsS82F
    rates[GPP] = rate_FppsS82F - rate_LimSyn
    rates[LIM] = rate_LimSyn 

    rates[ATP] = - rate_Hex - rate_PfkB + rate_Pgk + rate_PykF - rate_Mvk - rate_Pmvk - rate_Mdc
    rates[ADP] = - rates[ATP]
    rates[DHAP] = rate_Fba - rate_Tpi
    # rates[NAD] = - (1/3)*rate_Gap + rate_NoxE - rate_PDH
    rates[NAD] = - rate_Gap + rate_NoxE - rate_PDH
    rates[NADH] = -rates[NAD]
    # rates[NADP] = -(2/3)*rate_mGap + 2*rate_Hmgr
    rates[NADP] = -rate_mGap + 2*rate_Hmgr
    rates[NADPH] = -rates[NADP]
    rates[P_I] = -rate_Gap - rate_mGap + rate_Mdc + 2*rate_PPase
    # rates[CoA] = -rate_PDH + 0.5*rate_PhaA + rate_Hmgs + rate_Hmgr
    rates[CoA] = -rate_PDH + rate_PhaA + rate_Hmgs + rate_Hmgr
    rates[CO2] = rate_PDH + rate_Mdc
    rates[PPi] = rate_FppsS82F - rate_PPase + rate_LimSyn


    # rates[Glc] = -rate_Hex
    # rates[G6P] = rate_Hex - rate_Pgi
    # rates[F6P] = rate_Pgi - rate_PfkB
    # rates[FBP] = rate_PfkB - rate_Fba
    # rates[G3P] = 2*rate_Fba + 2*rate_Tpi - rate_Gap - rate_mGap
    # rates[BPG13] = rate_Gap + rate_mGap - rate_Pgk
    # rates[PG3] = rate_Pgk - rate_dPgm
    # rates[PG2] = rate_dPgm - rate_Eno
    # rates[PEP] = rate_Eno - rate_PykF
    # rates[PYR] = rate_PykF - rate_PDH
    # rates[AcCoA] = rate_PDH - rate_PhaA - rate_Hmgs
    # rates[AcAcCoA] = 0.5*rate_PhaA - rate_Hmgs
    # rates[HMGCoA] = rate_Hmgs - rate_Hmgr
    # rates[MEV] = rate_Hmgr - rate_Mvk
    # rates[MEVP] = rate_Mvk - rate_Pmvk
    # rates[MEVPP] = rate_Pmvk - rate_Mdc
    # rates[IPP] = rate_Mdc - rate_Idi - rate_FppsS82F
    # rates[DMAPP] = rate_Idi - rate_FppsS82F
    # rates[GPP] = rate_FppsS82F - rate_LimSyn
    # rates[LIM] = rate_LimSyn 

    # rates[ATP] = - rate_Hex - rate_PfkB + rate_Pgk + rate_PykF - rate_Mvk - rate_Pmvk - rate_Mdc
    # rates[ADP] = - rates[ATP]
    # rates[DHAP] = rate_Fba - rate_Tpi
    # rates[NAD] = - (1/3)*rate_Gap + rate_NoxE - rate_PDH
    # rates[NADH] = -rates[NAD]
    # rates[NADP] = -(2/3)*rate_mGap + 2*rate_Hmgr
    # rates[NADPH] = -rates[NADP]
    # rates[P_I] = -rate_Gap + rate_Mdc + 2*rate_PPase
    # rates[CoA] = -rate_PDH + 0.5*rate_PhaA + rate_Hmgs + rate_Hmgr
    # rates[CO2] = rate_PDH + rate_Mdc
    # rates[PPi] = rate_FppsS82F - rate_PPase
    



if __name__=="__main__":
    c = np.ones(TOTAL_SPECIES)
    rates = np.zeros(TOTAL_SPECIES)
    get_rxn_rates_pkb(c, rates)
    print(rates)




