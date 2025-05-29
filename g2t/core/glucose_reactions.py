"""
Author: Yudong Li (Yudong.Li@nrel.gov)
"""
from g2t.core.enzyme_kinetics import *
from .constants import *

hexokinase = [0.00091, 483.3] #2.7.1.1: hexokinase https://www.brenda-enzymes.org/enzyme.php?ecno=2.7.1.1#TURNOVER%20NUMBER%20[1/s]
pgi = [0.04, 2765] #5.3.1.9: glucose-6-phosphate isomerase  https://www.brenda-enzymes.org/enzyme.php?ecno=5.3.1.9#TURNOVER%20NUMBER%20[1/s]
pfk = [0.016, 4328] #2.7.1.11: 6-phosphofructokinase  https://www.brenda-enzymes.org/enzyme.php?ecno=2.7.1.11#TURNOVER%20NUMBER%20[1/s]
fba = [0.0002, 64.5] #4.1.2.13: fructose-bisphosphate aldolase https://www.brenda-enzymes.org/enzyme.php?ecno=4.1.2.13#TURNOVER%20NUMBER%20[1/s]
tpi = [0.06, 800] #5.3.1.1: triose-phosphate isomerase https://www.brenda-enzymes.org/enzyme.php?ecno=5.3.1.1#TURNOVER%20NUMBER%20[1/s]
gap = [0.002, 442] #1.2.1.12: glyceraldehyde-3-phosphate dehydrogenase (phosphorylating) https://www.brenda-enzymes.org/all_enzymes.php?ecno=1.2.1.12&table=Turnover_Number#TAB
mgap = [0.002, 442] #???? No info, use the same as gap 1.2.1.12: glyceraldehyde-3-phosphate dehydrogenase (phosphorylating)
pgk = [0.78, 2633] #EC 2.7.2.3, phosphoglycerate kinase https://www.brenda-enzymes.org/enzyme.php?ecno=2.7.2.3#TURNOVER%20NUMBER%20[1/s]
dpgm = [1.51, 3226] #EC 5.4.2.12, phosphoglycerate mutase https://www.brenda-enzymes.org/enzyme.php?ecno=5.4.2.12#TURNOVER%20NUMBER%20[1/s]
eno = [0.018, 230] #EC 4.2.1.11, phosphopyruvate hydratase https://www.brenda-enzymes.org/enzyme.php?ecno=4.2.1.11#TURNOVER%20NUMBER%20[1/s]
pykf = [0.38, 3204] #EC 2.7.1.40, pyruvate kinase https://www.brenda-enzymes.org/enzyme.php?ecno=2.7.1.40#TURNOVER%20NUMBER%20[1/s]
pdh = [0.077, 69.9] #EC 1.2.1.104, pyruvate dehydrogenase system https://www.brenda-enzymes.org/enzyme.php?ecno=1.2.1.104#TURNOVER%20NUMBER%20[1/s]
phaa = [2.1, 7] #EC 2.3.1.9, acetyl-CoA C-acetyltransferase https://www.brenda-enzymes.org/enzyme.php?ecno=2.3.1.9#TURNOVER%20NUMBER%20[1/s]
hmgs = [0.0017, 4.6] #EC 2.3.3.10, hydroxymethylglutaryl-CoA synthase https://www.brenda-enzymes.org/enzyme.php?ecno=2.3.3.10#TURNOVER%20NUMBER%20[1/s]
hmgr = [0.023, 0.023] #EC 1.1.1.34, hydroxymethylglutaryl-CoA reductase (NADPH) https://www.brenda-enzymes.org/enzyme.php?ecno=1.1.1.34#TURNOVER%20NUMBER%20[1/s]
mvk = [3.1, 143.8] # EC 2.7.1.36, mevalonate kinase https://www.brenda-enzymes.org/enzyme.php?ecno=2.7.1.36#TURNOVER%20NUMBER%20[1/s]
pmvk = [3.4, 10.2] #EC 2.7.4.2, phosphomevalonate kinase https://www.brenda-enzymes.org/enzyme.php?ecno=2.7.4.2#TURNOVER%20NUMBER%20[1/s] ATP 3.4 phosphomevalonate 10.2
mdc = [30, 40]  #https://www.sciencedirect.com/science/article/pii/S0141813015003827 APT 30/s mevalonate-pp 40/s
idi = [0.000197, 29900] #EC 5.3.3.2, isopentenyl-diphosphate DELTA-isomerase https://www.brenda-enzymes.org/enzyme.php?ecno=5.3.3.2#TURNOVER%20NUMBER%20[1/s]
fpps = [0.49, 33.5] #EC 2.5.1.10, (2E,6E)-farnesyl diphosphate synthase https://www.brenda-enzymes.org/enzyme.php?ecno=2.5.1.10#TURNOVER%20NUMBER%20[1/s]  DMAPP 5.73   This is probably not correct, but we will use this range for now.
limsyn = [0.012, 0.186] #EC 4.2.3.20, (R)-limonene synthase https://www.brenda-enzymes.org/enzyme.php?ecno=4.2.3.20#TURNOVER%20NUMBER%20[1/s]
noxe = [6.67, 428.3] #EC 1.6.3.4, NADH oxidase (H2O-forming) https://www.brenda-enzymes.org/enzyme.php?ecno=1.6.3.4#TURNOVER%20NUMBER%20[1/s]
ppase = [570, 1000] #https://www.sciencedirect.com/science/article/pii/S0021925820511113#:~:text=coli%20pyrophosphatase%20was%20found%20to,%E2%88%9216.6%20kcal%2Fmol).

def get_rxn_rates(c, rates):
    hex = Hex(packedBed=False)
    pgi = Pgi(packedBed=False)
    pfk = Pfk(packedBed=False)
    fba = Fba(packedBed=False)
    tpi = Tpi(packedBed=False)
    gap = Gap(packedBed=False)
    mgap = mGap(packedBed=False)
    noxe = NoxE(packedBed=False)
    pgk = Pgk(packedBed=False)
    pgm = Pgm(packedBed=False)
    eno = Eno(packedBed=False)
    pyk = Pyk(packedBed=False)
    pdh = Pdh(packedBed=False)
    phaa = PhaA(packedBed=False)
    hmgs = Hmgs(packedBed=False)
    hmgr = Hmgr(packedBed=False)
    mvk = Mvk(packedBed=False)
    pmvk = Pmvk(packedBed=False)
    mdc = Mdc(packedBed=False)
    idi = Idi(packedBed=False)
    fpps = Fpps_S82F(packedBed=False)
    ppase = PPase(packedBed=False)
    limsyn = LimSyn(packedBed=False)
    # Glc + ATP -> G6P + ADP
    # C6H12O6 + C10H16N5O13P3 -> C6H13O9P + C10H15N5O10P2
    rate_Hex = hex.get_rate(c[Glc], c[ATP], c[G6P], c[ADP])
    # G6P -> F6P
    # C6H13O9P -> C6H13O9P
    rate_Pgi = pgi.get_rate(c[G6P], c[F6P])
    # C6H13O9P + C10H16N5O13P3 -> C6H14O12P2 + C10H15N5O10P2
    rate_PfkB = pfk.get_rate(c[F6P], c[ATP], c[FBP], c[ADP])
    # C6H14O12P2 -> C3H7O6P + C3H7O6P
    rate_Fba = fba.get_rate(c[FBP], c[DHAP], c[G3P])
    # C3H7O6P -> C3H7O6P
    rate_Tpi = tpi.get_rate(c[DHAP], c[G3P])
    # 3 C3H7O6P + C21H27N7O14P2 + 3 PO4^{3-} + H+ + 2e- -> C21H27N7O14P2 + 3 C3H8O10P2
    # NAD -> NADH: NAD+ + H+ + 2e- = NADH (https://www.nad.com/nad-vs-nadh)
    rate_Gap = gap.get_rate(c[G3P], c[NAD], c[P_I], c[BPG13], c[NADH])
    rate_mGap = mgap.get_rate(c[G3P], c[NADP], c[P_I], c[BPG13], c[NADPH])
    # rate_Gap = gap.get_rate(c[G3P], c[NAD], c[PPi], c[BPG13], c[NADH])
    # rate_mGap = mgap.get_rate(c[G3P], c[NADP], c[PPi], c[BPG13], c[NADPH])
    rate_NoxE = noxe.get_rate(c[NADH])
    rate_Pgk = pgk.get_rate(c[BPG13], c[ADP], c[PG3], c[ATP])
    rate_dPgm = pgm.get_rate(c[PG3], c[PG2])
    rate_Eno = eno.get_rate(c[PG2], c[PEP])
    rate_PykF = pyk.get_rate(c[PEP], c[ADP], c[PYR], c[ATP])
    rate_PDH = pdh.get_rate(c[CoA], c[NAD], c[PYR], c[AcCoA], c[NADH], c[CO2])
    rate_PhaA = phaa.get_rate(c[AcCoA], c[AcAcCoA], c[CoA])
    rate_Hmgs = hmgs.get_rate(c[AcCoA], c[AcAcCoA], c[HMGCoA], c[CoA])
    rate_Hmgr = hmgr.get_rate(c[HMGCoA], c[NADPH], c[MEV], c[NADP], c[CoA])
    rate_Mvk = mvk.get_rate(c[MEV], c[ATP], c[MEVP], c[ADP])
    rate_Pmvk = pmvk.get_rate(c[MEVP], c[ATP], c[MEVPP], c[ADP])
    rate_Mdc = mdc.get_rate(c[MEVPP], c[ATP], c[IPP], c[ADP], c[CO2], c[P_I])
    rate_Idi = idi.get_rate(c[IPP], c[DMAPP])
    rate_FppsS82F = fpps.get_rate(c[IPP], c[DMAPP], c[GPP], c[PPi])
    rate_PPase = ppase.get_rate(c[PPi], c[P_I])
    rate_LimSyn = limsyn.get_rate(c[GPP], c[LIM], c[PPi])

    rates[Glc] = -rate_Hex if c[Glc] >= 0 else max(-rate_Hex, 0)
    rates[G6P] = rate_Hex - rate_Pgi if c[G6P] >= 0 else max(rate_Hex - rate_Pgi, 0)
    rates[F6P] = rate_Pgi - rate_PfkB if c[F6P] >= 0 else max(rate_Pgi - rate_PfkB, 0)
    rates[FBP] = rate_PfkB - rate_Fba if c[FBP] >= 0 else max(rate_PfkB - rate_Fba, 0)
    rates[G3P] = rate_Fba + rate_Tpi - rate_Gap - rate_mGap if c[G3P] >= 0 else max(rate_Fba + rate_Tpi - rate_Gap - rate_mGap, 0)
    rates[BPG13] = rate_Gap + rate_mGap - rate_Pgk if c[BPG13] >= 0 else max(rate_Gap + rate_mGap - rate_Pgk, 0)
    rates[PG3] = rate_Pgk - rate_dPgm if c[PG3] >= 0 else max(rate_Pgk - rate_dPgm, 0)
    rates[PG2] = rate_dPgm - rate_Eno if c[PG2] >= 0 else max(rate_dPgm - rate_Eno, 0)
    rates[PEP] = rate_Eno - rate_PykF if c[PEP] >= 0 else max(rate_Eno - rate_PykF, 0)
    rates[PYR] = rate_PykF - rate_PDH if c[PYR] >= 0 else max(rate_PykF - rate_PDH, 0)
    rates[AcCoA] = rate_PDH - rate_PhaA - rate_Hmgs if c[AcCoA] >= 0 else max(rate_PDH - rate_PhaA - rate_Hmgs, 0)
    rates[AcAcCoA] = 0.5*rate_PhaA - rate_Hmgs if c[AcAcCoA] >= 0 else max(0.5*rate_PhaA - rate_Hmgs, 0)
    rates[HMGCoA] = rate_Hmgs - rate_Hmgr if c[HMGCoA] >= 0 else max(rate_Hmgs - rate_Hmgr, 0)
    rates[MEV] = rate_Hmgr - rate_Mvk if c[MEV] >= 0 else max(rate_Hmgr - rate_Mvk, 0)
    rates[MEVP] = rate_Mvk - rate_Pmvk if c[MEVP] >= 0 else max(rate_Mvk - rate_Pmvk, 0)
    rates[MEVPP] = rate_Pmvk - rate_Mdc if c[MEVPP] >= 0 else max(rate_Pmvk - rate_Mdc, 0)
    rates[IPP] = rate_Mdc - rate_Idi - rate_FppsS82F if c[IPP] >= 0 else max(rate_Mdc - rate_Idi - rate_FppsS82F, 0)
    rates[DMAPP] = rate_Idi - rate_FppsS82F if c[DMAPP] >= 0 else max(rate_Idi - rate_FppsS82F, 0)
    rates[GPP] = rate_FppsS82F - rate_LimSyn if c[GPP] >= 0 else max(rate_FppsS82F - rate_LimSyn, 0)
    rates[LIM] = rate_LimSyn if c[LIM] >= 0 else max(rate_LimSyn, 0)

    rates[ATP] = - rate_Hex - rate_PfkB + rate_Pgk + rate_PykF - rate_Mvk - rate_Pmvk - rate_Mdc if c[ATP] >= 0 else max(- rate_Hex - rate_PfkB + rate_Pgk + rate_PykF - rate_Mvk - rate_Pmvk - rate_Mdc, 0)
    rates[ADP] = - rates[ATP] if c[ADP] >= 0 else max(- rates[ATP], 0)
    rates[DHAP] = rate_Fba - rate_Tpi if c[DHAP] >= 0 else max(rate_Fba - rate_Tpi, 0)
    rates[NAD] = - (1/3)*rate_Gap + rate_NoxE - rate_PDH if c[NAD] >= 0 else max(- (1/3)*rate_Gap + rate_NoxE - rate_PDH, 0)
    rates[NADH] = -rates[NAD] if c[NADH] >= 0 else max(-rates[NAD], 0)
    rates[NADP] = -(2/3)*rate_mGap + 2*rate_Hmgr if c[NADP] >= 0 else max(-(2/3)*rate_mGap + 2*rate_Hmgr, 0)
    rates[NADPH] = -rates[NADP] if c[NADPH] >= 0 else max(-rates[NADP] , 0)
    rates[P_I] = -rate_Gap - rate_mGap + rate_Mdc + 2*rate_PPase if c[P_I] >= 0 else max(-rate_Gap - rate_mGap + rate_Mdc + 2*rate_PPase, 0)
    rates[CoA] = -rate_PDH + 0.5*rate_PhaA + rate_Hmgs + rate_Hmgr if c[CoA] >= 0 else max(-rate_PDH + 0.5*rate_PhaA + rate_Hmgs + rate_Hmgr, 0)
    rates[CO2] = rate_PDH + rate_Mdc if c[CO2] >= 0 else max(rate_PDH + rate_Mdc, 0)
    rates[PPi] = rate_FppsS82F - rate_PPase + rate_LimSyn if c[PPi] >= 0 else max(rate_FppsS82F - rate_PPase + rate_LimSyn, 0)



class EnzymeSystem:

    def __init__(self, fitExp=True, packedBed=False):
        self.l_enzyme_names = ['Hex', 'Pgi', 'Pfk', 'Fba', 'Tpi', 'Gap', 'mGap', 'Pgk', 'Pgm', 'Eno', 'Pyk', 'Pdh', 'PhaA', 'Hmgs', 'Hmgr', 'Mvk', 'Pmvk', 'Mdc', 'Idi', 'Fpps_S82F', 'LimSyn', 'PPase', 'NoxE']
        self.l_enzyme_cls = [Hex, Pgi, Pfk, Fba, Tpi, Gap, mGap, Pgk, Pgm, Eno, Pyk, Pdh, PhaA, Hmgs, Hmgr, Mvk, Pmvk, Mdc, Idi, Fpps_S82F, LimSyn, PPase, NoxE]
        self.l_mg = [0.001,0.03,0.02,0.04,0.004,0.007,0.03,0.01,0.05,0.05,0.01,0.01,0.02,0.03,0.08,0.02,0.04,0.08,0.07,0.02,0.6,0.2,0.05]
        self.l_U = [0.01,0.6,0.2,0.9,0.4,0.1,0.4,1.1,4.6,1.4,2.9,0.2,1.3,0.04,0.4,0.2,0.6,0.3,0.3,0.1,0.006,0.1,2.1]
        self.l_enzyme_kcat_ranges = [hexokinase, pgi, pfk, fba, tpi, gap, mgap, pgk, dpgm, eno, pykf, pdh, phaa, hmgs, hmgr, mvk, pmvk, mdc, idi, fpps, limsyn, ppase, noxe]
        
        self.kcat_fitted = {
            'Hex': 49781.02365551095,
            'Pgi': 73000.0,
            'Pfk': 106000.0,
            'Fba': 139000.0,
            'Tpi': 172000.0,
            'Gap': 205000.0,
            'mGap': 238000.0,
            'Pgk': 271000.0,
            'Pgm': 304000.0,
            'Eno': 337000.0,
            'Pyk': 370000.0,
            'Pdh': 403000.0,
            'PhaA': 436000.0,
            'Hmgs': 469000.0,
            'Hmgr': 502000.0,
            'Mvk': 535000.0,
            'Pmvk': 568000.0,
            'Mdc': 601000.0,
            'Idi': 634000.0,
            'Fpps_S82F': 667000.0,
            'LimSyn': 700000.0,
            'PPase': 2206.3524049873267,
            'NoxE': 1136.3479611885175,
            }
        # self.kcats_fitted = {
        #     'Hex': 49615, # very sensitive, Li params not okay
        #     'Pgi':73000,  # not sensitive to 1 order of magnitude, changes for 2 OoM less, Li okay
        #     'Pfk':106000, # not sensitive to 1 order of magnitude, changes for 2 OoM less, Li okay
        #     'Fba':139000, # no change with 2 OoM, Li okay
        #     'Tpi':172000, # no change 1 OoM, changes for 2 OoM less, Li okay
        #     'Gap':205000, # changes with 1 OoM less, no change with increase up to 2 OoM, Li okay
        #     'mGap':238000, # changes for 2 OoM less, no change with increase up to 2 OoM, Li okay
        #     'Pgk':271000, # changes for 2 OoM less, no change with increase up to 2 OoM, Li okay
        #     'Pgm':304000, # no change with 2 OoM, Li okay
        #     'Eno':337000, # changes for 2 OoM less, no change with increase up to 2 OoM, Li okay
        #     'Pyk':370000, # changes for 2 OoM less, no change with increase up to 2 OoM, Li okay
        #     'Pdh':403000, # no change with 2 OoM, Li okay
        #     'PhaA':436000, # small change for 2 OoM less, no change with increase up to 2 OoM, Li okay
        #     'Hmgs':469000, # no change with 2 OoM, Li okay
        #     'Hmgr':502000, # no change with 2 OoM, Li okay
        #     'Mvk':535000, # no change with 2 OoM, Li okay
        #     'Pmvk':568000, # no change with 2 OoM, Li okay
        #     'Mdc':601000, # no change with 2 OoM, Li okay
        #     'Idi':634000, # no change with 2 OoM, Li okay
        #     'Fpps_S82F':667000, # no change with 2 OoM, Li okay
        #     'LimSyn':700000, # no change with 2 OoM, Li okay
        #     'NoxE':1136, # changes with 1 OoM decrease, no change with increase up to 2 OoM, Li okay
        #     'PPase':2206, # changes with 1 OoM decrease, no change with increase up to 2 OoM, Li okay
        # }
        self.enzymes = {}
        self.enzyme_Kcat_range = {}
        
        for name, e, mg, U, Kcat_range in zip(self.l_enzyme_names, self.l_enzyme_cls, self.l_mg, self.l_U, self.l_enzyme_kcat_ranges):
            en = e(packedBed=packedBed)
            self.enzyme_Kcat_range[name] = Kcat_range
            if packedBed:
                en.set_enzyme_concentration_with_mg(mg, 800)
            if fitExp:
                en.set_enzyme_conc_with_U(U, 800)
            self.enzymes[name] = en
            # self.enzyme.append(en)
            # print(f"Enzyme {en.name} has U computed from Kcat[E]={en.c*en.Kcat:0.2f} compared with table data: {en.U} ")
        if packedBed:
            for name, value in self.kcat_fitted.items():
                self.set_enzyme_Kcat(name, value)
    
    def set_enzyme_Kcat(self, name, Kcat):
        """
        sets Kcat for enzyme name
        """
        self.enzymes[name].set_Kcat(Kcat)

    def get_rate(self, c, rates):
        hex,pgi,pfk,fba,tpi,gap,mgap,pgk,pgm,eno,pyk,pdh,phaa,hmgs,hmgr,mvk,pmvk,mdc,idi,fpps,limsyn,ppase,noxe = [self.enzymes[name] for name in self.l_enzyme_names]
        # Glc + ATP -> G6P + ADP
        # C6H12O6 + C10H16N5O13P3 -> C6H13O9P + C10H15N5O10P2
        rate_Hex = hex.get_rate(c[Glc], c[ATP], c[G6P], c[ADP])
        # G6P -> F6P
        # C6H13O9P -> C6H13O9P
        rate_Pgi = pgi.get_rate(c[G6P], c[F6P])
        # C6H13O9P + C10H16N5O13P3 -> C6H14O12P2 + C10H15N5O10P2
        rate_PfkB = pfk.get_rate(c[F6P], c[ATP], c[FBP], c[ADP])
        # C6H14O12P2 -> C3H7O6P + C3H7O6P
        rate_Fba = fba.get_rate(c[FBP], c[DHAP], c[G3P])
        # C3H7O6P -> C3H7O6P
        rate_Tpi = tpi.get_rate(c[DHAP], c[G3P])
        # 3 C3H7O6P + C21H27N7O14P2 + 3 PO4^{3-} + H+ + 2e- -> C21H27N7O14P2 + 3 C3H8O10P2
        # NAD -> NADH: NAD+ + H+ + 2e- = NADH (https://www.nad.com/nad-vs-nadh)
        rate_Gap = gap.get_rate(c[G3P], c[NAD], c[P_I], c[BPG13], c[NADH])
        rate_mGap = mgap.get_rate(c[G3P], c[NADP], c[P_I], c[BPG13], c[NADPH])
        # rate_Gap = gap.get_rate(c[G3P], c[NAD], c[PPi], c[BPG13], c[NADH])
        # rate_mGap = mgap.get_rate(c[G3P], c[NADP], c[PPi], c[BPG13], c[NADPH])
        rate_NoxE = noxe.get_rate(c[NADH])
        rate_Pgk = pgk.get_rate(c[BPG13], c[ADP], c[PG3], c[ATP])
        rate_dPgm = pgm.get_rate(c[PG3], c[PG2])
        rate_Eno = eno.get_rate(c[PG2], c[PEP])
        rate_PykF = pyk.get_rate(c[PEP], c[ADP], c[PYR], c[ATP])
        rate_PDH = pdh.get_rate(c[CoA], c[NAD], c[PYR], c[AcCoA], c[NADH], c[CO2])
        rate_PhaA = phaa.get_rate(c[AcCoA], c[AcAcCoA], c[CoA])
        rate_Hmgs = hmgs.get_rate(c[AcCoA], c[AcAcCoA], c[HMGCoA], c[CoA])
        rate_Hmgr = hmgr.get_rate(c[HMGCoA], c[NADPH], c[MEV], c[NADP], c[CoA])
        rate_Mvk = mvk.get_rate(c[MEV], c[ATP], c[MEVP], c[ADP])
        rate_Pmvk = pmvk.get_rate(c[MEVP], c[ATP], c[MEVPP], c[ADP])
        rate_Mdc = mdc.get_rate(c[MEVPP], c[ATP], c[IPP], c[ADP], c[CO2], c[P_I])
        rate_Idi = idi.get_rate(c[IPP], c[DMAPP])
        rate_FppsS82F = fpps.get_rate(c[IPP], c[DMAPP], c[GPP], c[PPi])
        rate_PPase = ppase.get_rate(c[PPi], c[P_I])
        rate_LimSyn = limsyn.get_rate(c[GPP], c[LIM], c[PPi])

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

    def get_packedBed_rate(self, c, enzyme_concentration:dict, rates):
        hex,pgi,pfk,fba,tpi,gap,mgap,pgk,pgm,eno,pyk,pdh,phaa,hmgs,hmgr,mvk,pmvk,mdc,idi,fpps,limsyn,ppase,noxe = [self.enzymes[name] for name in self.l_enzyme_names]
        # Glc + ATP -> G6P + ADP
        # C6H12O6 + C10H16N5O13P3 -> C6H13O9P + C10H15N5O10P2
        rate_Hex = hex.get_rate(c[Glc], c[ATP], c[G6P], c[ADP])*enzyme_concentration["Hex"]
        # G6P -> F6P
        # C6H13O9P -> C6H13O9P
        rate_Pgi = pgi.get_rate(c[G6P], c[F6P])*enzyme_concentration["Pgi"]
        # C6H13O9P + C10H16N5O13P3 -> C6H14O12P2 + C10H15N5O10P2
        rate_PfkB = pfk.get_rate(c[F6P], c[ATP], c[FBP], c[ADP])*enzyme_concentration["Pfk"]
        # C6H14O12P2 -> C3H7O6P + C3H7O6P
        rate_Fba = fba.get_rate(c[FBP], c[DHAP], c[G3P])*enzyme_concentration["Fba"]
        # C3H7O6P -> C3H7O6P
        rate_Tpi = tpi.get_rate(c[DHAP], c[G3P])*enzyme_concentration["Tpi"]
        # 3 C3H7O6P + C21H27N7O14P2 + 3 PO4^{3-} + H+ + 2e- -> C21H27N7O14P2 + 3 C3H8O10P2
        # NAD -> NADH: NAD+ + H+ + 2e- = NADH (https://www.nad.com/nad-vs-nadh)
        rate_Gap = gap.get_rate(c[G3P], c[NAD], c[P_I], c[BPG13], c[NADH])*enzyme_concentration["Gap"]
        rate_mGap = mgap.get_rate(c[G3P], c[NADP], c[P_I], c[BPG13], c[NADPH])*enzyme_concentration["mGap"]
        # rate_Gap = gap.get_rate(c[G3P], c[NAD], c[PPi], c[BPG13], c[NADH])
        # rate_mGap = mgap.get_rate(c[G3P], c[NADP], c[PPi], c[BPG13], c[NADPH])
        rate_NoxE = noxe.get_rate(c[NADH])*enzyme_concentration["NoxE"]
        rate_Pgk = pgk.get_rate(c[BPG13], c[ADP], c[PG3], c[ATP])*enzyme_concentration["Pgk"]
        rate_dPgm = pgm.get_rate(c[PG3], c[PG2])*enzyme_concentration["Pgm"]
        rate_Eno = eno.get_rate(c[PG2], c[PEP])*enzyme_concentration["Eno"]
        rate_PykF = pyk.get_rate(c[PEP], c[ADP], c[PYR], c[ATP])*enzyme_concentration["Pyk"]
        rate_PDH = pdh.get_rate(c[CoA], c[NAD], c[PYR], c[AcCoA], c[NADH], c[CO2])*enzyme_concentration["Pdh"]
        rate_PhaA = phaa.get_rate(c[AcCoA], c[AcAcCoA], c[CoA])*enzyme_concentration["PhaA"]
        rate_Hmgs = hmgs.get_rate(c[AcCoA], c[AcAcCoA], c[HMGCoA], c[CoA])*enzyme_concentration["Hmgs"]
        rate_Hmgr = hmgr.get_rate(c[HMGCoA], c[NADPH], c[MEV], c[NADP], c[CoA])*enzyme_concentration["Hmgr"]
        rate_Mvk = mvk.get_rate(c[MEV], c[ATP], c[MEVP], c[ADP])*enzyme_concentration["Mvk"]
        rate_Pmvk = pmvk.get_rate(c[MEVP], c[ATP], c[MEVPP], c[ADP])*enzyme_concentration["Pmvk"]
        rate_Mdc = mdc.get_rate(c[MEVPP], c[ATP], c[IPP], c[ADP], c[CO2], c[P_I])*enzyme_concentration["Mdc"]
        rate_Idi = idi.get_rate(c[IPP], c[DMAPP])*enzyme_concentration["Idi"]
        rate_FppsS82F = fpps.get_rate(c[IPP], c[DMAPP], c[GPP], c[PPi])*enzyme_concentration["Fpps_S82F"]
        rate_PPase = ppase.get_rate(c[PPi], c[P_I])*enzyme_concentration["PPase"]
        rate_LimSyn = limsyn.get_rate(c[GPP], c[LIM], c[PPi])*enzyme_concentration["LimSyn"]

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


if __name__=="__main__":
    import numpy as np
    c = np.ones(TOTAL_SPECIES)
    rates = np.zeros(TOTAL_SPECIES)
    get_rxn_rates(c, rates)
    print(rates)

    c = np.ones(TOTAL_SPECIES)
    rates = np.zeros(TOTAL_SPECIES)
    enzSys = EnzymeSystem()
    enzSys.get_rate(c, rates)
    print(rates)



