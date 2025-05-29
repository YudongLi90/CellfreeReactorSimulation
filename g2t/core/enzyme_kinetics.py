"""
Author: Yudong Li (Yudong.Li@nrel.gov)
"""
from g2t.core.kinetic_rate_eqs import *
from g2t.utils.returndecs import non_negative_return
import numpy as np

class Enzyme:

    
    def __init__(self, c=None, packedBed=True):
        self.c = c
        self.name = self.__class__.__name__
        self.conc_U_mg = None
        self.Kcat = None
        self.reversible = True
    # def set_concentration(self, c):
    #     self.U = self.U/(c/self.M)
    
    def set_enzyme_concentration_with_mg(self, mg, microLiter):
        # enzyme concentration with unit mol/L
        self.c = mg*1e-3/(self.M*1e3)/(microLiter*1e-6)
        # print(f"c: {self.c}")
    
    def set_enzyme_weightconcentration_with_mg(self, mg, microLiter):
        self.c = mg*1e-3/(microLiter*1e-6)

    def set_enzyme_conc_with_U(self, U, microLiter):
        mg = U/self.conc_U_mg
        self.set_enzyme_concentration_with_mg(mg, microLiter)

    def set_Kcat(self, Kcat):
        # literature doesn't give unit for Kcat
        self.Kcat = Kcat
    
    def get_U(self):
        if self.c:
            return self.Kcat*self.c
        else:
            return self.U 

class Hex(Enzyme):
    # Biologically irreversible
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 0.05 # mM/s
        self.Kcat = 49615.3784209611 #4e4 #4e6#2.6e6
        self.M = 100 # 1e3 g/mol
        self.ATP=0.1
        self.Glc=0.15
        self.ADP=0.1
        self.G6P=0.1
        self.Keq=1000
        self.conc_U_mg = 7.2 
        if c:
            self.K = self.U/(c/self.M)  
        else:
            self.K = None
        self.packedBed = packedBed
        self.reversible = False

    # @non_negative_return
    def get_rate(self, c_Glc, c_ATP, c_G6P, c_ADP):
        if not self.reversible:
            c_G6P = 0
            c_ADP = 0
        if self.packedBed:
            if self.Kcat:
                return rate_2_2(c_Glc, c_ATP, c_G6P, c_ADP, self.Kcat*self.c, self.Keq, self.Glc, self.ATP, self.G6P, self.ADP)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_2_2(c_Glc, c_ATP, c_G6P, c_ADP, U, self.Keq, self.Glc, self.ATP, self.G6P, self.ADP)

class Pgi(Enzyme):
     
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 1
        self.Kcat = 7.3e4 #1.68e6
        self.M = 63
        self.G6P = 0.3
        self.F6P = 0.15
        self.Keq = 0.5
        self.conc_U_mg = 19.2
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_G6P, c_F6P):
        if self.packedBed:
            if self.Kcat:
                return rate_1_1(c_G6P, c_F6P, self.Kcat*self.c, self.Keq, self.G6P, self.F6P)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_1_1(c_G6P, c_F6P, U, self.Keq, self.G6P, self.F6P)

class Pfk(Enzyme):
    # reaction not reversible https://en.wikipedia.org/wiki/Phosphofructokinase
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 0.5
        self.Kcat = 1.06e5 #1.66e6
        self.M = 83
        self.ATP = 0.1
        self.F6P = 0.03
        self.ADP = 0.1
        self.FBP = 0.1
        self.Keq = 300
        self.conc_U_mg = 9.6
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed
        self.reversible = False
        
    def get_rate(self, c_F6P, c_ATP, c_FBP, c_ADP):
        if not self.reversible:
            c_FBP = 0
            c_ADP = 0
        if self.packedBed:
            if self.Kcat:
                return rate_2_2(c_F6P, c_ATP, c_FBP, c_ADP, self.Kcat*self.c, self.Keq, self.F6P, self.ATP, self.FBP, self.ADP)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_2_2(c_F6P, c_ATP, c_FBP, c_ADP, U, self.Keq, self.F6P, self.ATP, self.FBP, self.ADP)

class Fba(Enzyme):
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 1
        self.Kcat = 1.39e5 #6.6e5
        self.M = 33
        self.FBP = 0.015
        self.DHAP = 0.1
        self.G3P = 0.1
        self.Keq = 0.2
        self.conc_U_mg = 21.1
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_FBP, c_DHAP, c_G3P):
        if self.packedBed:
            if self.Kcat:
                return rate_1_2(c_FBP, c_DHAP, c_G3P, self.Kcat*self.c, self.Keq, self.FBP, self.DHAP, self.G3P)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_1_2(c_FBP, c_DHAP, c_G3P, U, self.Keq, self.FBP, self.DHAP, self.G3P)
    
class Tpi(Enzyme):

    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 0.25
        self.Kcat = 1.72e5 #2.7e6
        self.M = 54
        self.DHAP = 1
        self.G3P = 1
        self.Keq = 0.04
        self.conc_U_mg = 96.6
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_DHAP, c_G3P):
        if self.packedBed:
            if self.Kcat:
                return rate_1_1(c_DHAP, c_G3P, self.Kcat*self.c, self.Keq, self.DHAP, self.G3P)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_1_1(c_DHAP, c_G3P, U, self.Keq, self.DHAP, self.G3P)

class Gap(Enzyme):
    
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 0.5
        self.Kcat = 2.05e5 #8e6
        self.M = 140
        self.G3P = 1
        self.NAD = 0.05
        self.PI = 0.5
        self.BPG13 = 2
        self.NADH = 0.1
        self.Keq = 0.07
        self.conc_U_mg = 9.3
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_G3P, c_NAD, c_PI, c_BPG13, c_NADH):
        if self.packedBed:
            if self.Kcat:
                return rate_3_2(c_G3P, c_NAD, c_PI, c_BPG13, c_NADH, self.Kcat*self.c, self.Keq, self.G3P, self.NAD, self.PI, self.BPG13, self.NADH)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_3_2(c_G3P, c_NAD, c_PI, c_BPG13, c_NADH, U, self.Keq, self.G3P, self.NAD, self.PI, self.BPG13, self.NADH)

class mGap(Enzyme):

    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 10
        self.Kcat = 2.38e5 #3.73e6
        self.M = 14 # this is assumed to be 1/10th of Gap
        self.Keq = 0.07
        
        self.G3P = 1
        self.NADP = 0.05
        self.PI = 0.5

        self.BPG13 = 2
        self.NADPH = 0.1

        self.conc_U_mg = 12.9
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_G3P, c_NADP, c_PI, c_BPG13, c_NADPH):
        if self.packedBed:
            if self.Kcat:
                return rate_3_2(c_G3P, c_NADP, c_PI, c_BPG13, c_NADPH, self.Kcat*self.c, self.Keq, self.G3P, self.NADP, self.PI, self.BPG13, self.NADPH)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_3_2(c_G3P, c_NADP, c_PI, c_BPG13, c_NADPH, U, self.Keq, self.G3P, self.NADP, self.PI, self.BPG13, self.NADPH)
    
class Pgk(Enzyme):

    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 2
        self.Kcat = 2.71e5 #7.2e6
        self.M = 45
        self.Keq = 3200
            
        self.BPG13 = 0.1
        self.ADP = 0.1

        self.PG3 = 2.2
        self.ATP = 3

        self.conc_U_mg = 102.3
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_BPG13, c_ADP, c_PG3, c_ATP):
        if self.packedBed:
            if self.Kcat:
                return rate_2_2(c_BPG13, c_ADP, c_PG3, c_ATP, self.Kcat*self.c, self.Keq, self.BPG13, self.ADP, self.PG3, self.ATP)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_2_2(c_BPG13, c_ADP, c_PG3, c_ATP, U, self.Keq, self.BPG13, self.ADP, self.PG3, self.ATP)

class Pgm(Enzyme):

    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 1
        self.Kcat = 3.04e5 #1.05e6
        self.M = 65.69
        self.Keq = 0.19

        self.PG3 = 0.7 #0.5

        self.PG2 = 0.04 #0.1

        self.conc_U_mg=97.2
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_PG3, c_PG2):
        if self.packedBed:
            if self.Kcat:
                return rate_1_1(c_PG3, c_PG2, self.Kcat*self.c, self.Keq, self.PG3, self.PG2)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_1_1(c_PG3, c_PG2, U, self.Keq, self.PG3, self.PG2)
    
class Eno(Enzyme):

    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 1
        self.Kcat = 3.37e5 #1.312e6
        self.M = 82
        self.Keq = 6.7

        self.PG2 = 0.1

        self.PEP = 1

        self.conc_U_mg = 82.5
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_PG2, c_PEP):
        if self.packedBed:
            if self.Kcat:
                return rate_1_1(c_PG2, c_PEP, self.Kcat*self.c, self.Keq, self.PG2, self.PEP)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_1_1(c_PG2, c_PEP, U, self.Keq, self.PG2, self.PEP)

class Pyk(Enzyme):
    # irreversible reaction
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 1
        self.Kcat = 3.7e5 #1.896e7
        self.M = 237
        self.Keq = 6500 #5000

        self.PEP = 0.3
        self.ADP = 0.3

        self.PYR = 1
        self.ATP = 1

        self.conc_U_mg = 365.4
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_PEP, c_ADP, c_PYR, c_ATP):
        if self.packedBed:
            if self.Kcat:
                return rate_2_2(c_PEP, c_ADP, c_PYR, c_ATP, self.Kcat*self.c, self.Keq, self.PEP, self.ADP, self.PYR, self.ATP)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_2_2(c_PEP, c_ADP, c_PYR, c_ATP, U, self.Keq, self.PEP, self.ADP, self.PYR, self.ATP)

class Pdh(Enzyme):
    # irreversible reaction
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 10
        self.Kcat = 4.03e5 #1.4e8
        self.M = 175
        self.Keq = 6500 #5000

        self.CoA = 0.1
        self.NAD = 0.1
        self.PYR = 0.5

        self.AcCoA = 0.1
        self.NADH = 0.1
        self.CO2 = 0.1

        self.conc_U_mg = 1.4
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed
        self.reversible = False

    # @non_negative_return
    def get_rate(self, c_CoA, c_NAD, c_PYR, c_AcCoA, c_NADH, c_CO2):
        if not self.reversible:
            c_AcCoA = 0
            c_NADH = 0
            c_CO2 = 0
        if self.packedBed:
            if self.Kcat:
                return rate_3_3(c_CoA, c_NAD, c_PYR, c_AcCoA, c_NADH, c_CO2, self.Kcat*self.c, self.Keq, self.CoA, self.NAD, self.PYR, self.AcCoA, self.NADH, self.CO2)

            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_3_3(c_CoA, c_NAD, c_PYR, c_AcCoA, c_NADH, c_CO2, U, self.Keq, self.CoA, self.NAD, self.PYR, self.AcCoA, self.NADH, self.CO2)

class PhaA(Enzyme):
    
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 1
        self.Kcat = 4.36e5 #6.08e6
        self.M = 152
        self.Keq = 0.05

        self.AcCoA = 0.4

        self.AcAcCoA = 0.1
        self.CoA = 0.1

        self.conc_U_mg = 81.1
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_AcCoA, c_AcAcCoA, c_CoA):
        if self.packedBed:
            if self.Kcat:
                return rate_2_2(c_AcCoA, c_AcCoA, c_AcAcCoA, c_CoA, self.Kcat*self.c, self.Keq,
                        self.AcCoA, self.AcCoA, self.AcAcCoA, self.CoA)

            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_2_2(c_AcCoA, c_AcCoA, c_AcAcCoA, c_CoA, U, self.Keq,
                        self.AcCoA, self.AcCoA, self.AcAcCoA, self.CoA)
    
class Hmgs(Enzyme):
    # riversible https://en.wikipedia.org/wiki/Hydroxymethylglutaryl-CoA_synthase
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 1
        self.Kcat = 4.69e5 #2.82e6
        self.M = 105.8
        self.Keq = 50
        self.AcCoA = 0.35 #0.4
        self.AcAcCoA = 0.01

        self.HMGCoA = 0.5
        self.CoA = 0.5

        self.conc_U_mg = 1.5
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_AcCoA, c_AcAcCoA, c_HMGCoA, c_CoA):
        if self.packedBed:
            if self.Kcat:
                return rate_2_2(c_AcCoA, c_AcAcCoA, c_HMGCoA, c_CoA, self.Kcat*self.c, self.Keq, self.AcCoA, self.AcAcCoA, self.HMGCoA, self.CoA)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_2_2(c_AcCoA, c_AcAcCoA, c_HMGCoA, c_CoA, U, self.Keq, self.AcCoA, self.AcAcCoA, self.HMGCoA, self.CoA)

class Hmgr(Enzyme):
    # reversible https://www.sciencedirect.com/science/article/pii/B9780080912837000357
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 1
        self.Kcat = 5.02e5 #1.014e6
        self.M = 101.4
        self.Keq = 50
        self.HMGCoA = 0.02
        self.NADPH = 0.1

        self.MEV = 1
        self.NADP = 0.3
        self.CoA = 0.3

        self.conc_U_mg = 4.2
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_HMGCoA, c_NADPH, c_MEV, c_NADP, c_CoA):
        if self.packedBed:
            if self.Kcat:
                return rate_3_4(c_HMGCoA, c_NADPH, c_NADPH, c_MEV, c_CoA, c_NADP, c_NADP, self.Kcat*self.c, self.Keq, 
                        self.HMGCoA, self.NADPH, self.NADPH, self.MEV, self.CoA, self.NADP, self.NADP)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_3_4(c_HMGCoA, c_NADPH, c_NADPH, c_MEV, c_CoA, c_NADP, c_NADP, U, self.Keq, 
                        self.HMGCoA, self.NADPH, self.NADPH, self.MEV, self.CoA, self.NADP, self.NADP)
    
class Mvk(Enzyme):
    # reversibility unclear, assume reversible
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 5
        self.Kcat = 5.35e5 #6.62e6
        self.M = 33.1
        self.Keq = 653 #10
        self.MEV = 0.07
        self.ATP = 0.5

        self.MEVP = 1.5
        self.ADP = 0.1

        self.conc_U_mg = 8.1
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed
        self.reversible = False

    # @non_negative_return
    def get_rate(self, c_MEV, c_ATP, c_MEVP, c_ADP):
        if not self.reversible:
            c_MEVP = 0
            c_ADP = 0
        if self.packedBed:
            if self.Kcat:
                return rate_mvk(c_MEV, c_ATP, c_MEVP, c_ADP, self.Kcat*self.c, self.Keq, self.MEV, self.ATP, self.MEVP, self.ADP)
                # return rate_2_2(c_MEV, c_ATP, c_MEVP, c_ADP, self.Kcat*self.c, self.Keq, self.MEV, self.ATP, self.MEVP, self.ADP)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_mvk(c_MEV, c_ATP, c_MEVP, c_ADP, U, self.Keq, self.MEV, self.ATP, self.MEVP, self.ADP)
            # return rate_2_2(c_MEV, c_ATP, c_MEVP, c_ADP, U, self.Keq, self.MEV, self.ATP, self.MEVP, self.ADP)


class Pmvk(Enzyme):
    # reversible
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 5
        self.Kcat = 5.68e5 #2.2e6
        self.M = 22
        self.Keq = 2
        self.MEVP = 0.0077 #0.008
        self.ATP = 0.137 #0.14

        self.MEVPP = 0.014 #0.01
        self.ADP = 0.41 #0.4

        self.conc_U_mg = 14.8
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed
        self.reversible = False

    # @non_negative_return
    def get_rate(self, c_MEVP, c_ATP, c_MEVPP, c_ADP):
        if not self.reversible:
            c_MEVPP = 0
            c_ADP = 0
        if self.packedBed:
            if self.Kcat:
                return rate_mvk(c_MEVP, c_ATP, c_MEVPP, c_ADP, self.Kcat*self.c, self.Keq, self.MEVP, self.ATP, self.MEVPP, self.ADP)
                # return rate_2_2(c_MEVP, c_ATP, c_MEVPP, c_ADP, self.Kcat*self.c, self.Keq, self.MEVP, self.ATP, self.MEVPP, self.ADP)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_mvk(c_MEVP, c_ATP, c_MEVPP, c_ADP, U, self.Keq, self.MEVP, self.ATP, self.MEVPP, self.ADP)
            # return rate_2_2(c_MEVP, c_ATP, c_MEVPP, c_ADP, self.Kcat*self.c, self.Keq, self.MEVP, self.ATP, self.MEVPP, self.ADP)

class Mdc(Enzyme):
    # assume not reversible though wiki shows with reversible sign https://en.wikipedia.org/wiki/Diphosphomevalonate_decarboxylase
    # because of the involvement of CO2
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 5
        self.Kcat = 6.01e5 #2.17e6
        self.M = 43.4
        self.Keq = 10
        self.MEVPP = 0.1
        self.ATP = 0.06

        self.IPP = 0.5
        self.ADP = 0.3
        self.CO2 = 5
        self.PI = 0.5

        self.conc_U_mg = 4.1
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed
        self.reversible = False
    
    # @non_negative_return
    def get_rate(self, c_MEVPP, c_ATP, c_IPP, c_ADP, c_CO2, c_PI):
        if not self.reversible:
            c_IPP = 0
            c_ADP = 0
            c_CO2 = 0
            c_PI = 0
        if self.packedBed:
            if self.Kcat:
                return rate_2_4(c_MEVPP, c_ATP, c_IPP, c_ADP, c_CO2, c_PI, self.Kcat*self.c, self.Keq,
                        self.MEVPP, self.ATP, self.IPP, self.ADP, self.CO2, self.PI)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_2_4(c_MEVPP, c_ATP, c_IPP, c_ADP, c_CO2, c_PI, U, self.Keq,
                        self.MEVPP, self.ATP, self.IPP, self.ADP, self.CO2, self.PI)

class Idi(Enzyme):
    # reversible https://www.sciencedirect.com/science/article/abs/pii/B9780128239759000055
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 1
        self.Kcat = 6.34e5 #3.71257e5
        self.M = 32.485 # https://hmdb.ca/proteins/HMDBP01541
        self.Keq = 0.77
        self.IPP = 0.0035 #0.004

        self.DMAPP = 0.05

        self.conc_U_mg = 4.3
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    def get_rate(self, c_IPP, c_DMAPP):
        if self.packedBed:
            if self.Kcat:
                return rate_1_1(c_IPP, c_DMAPP, self.Kcat*self.c, self.Keq, self.IPP, self.DMAPP)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_1_1(c_IPP, c_DMAPP, U, self.Keq, self.IPP, self.DMAPP)
    
class Fpps_S82F(Enzyme):

    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 1
        self.Kcat = 6.67e5 #1.6e6
        self.M = 40 #https://www.frontiersin.org/articles/10.3389/fpls.2022.946629
        self.Keq = 200
        self.IPP = 0.005
        self.DMAPP = 0.005

        self.GPP = 4
        self.PPi = 5

        self.conc_U_mg = 7.4
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed
        self.reversible = False

    # @non_negative_return
    def get_rate(self, c_IPP, c_DMAPP, c_GPP, c_PPi):
        if not self.reversible:
            c_GPP = 0
            c_PPi = 0
        if self.packedBed:
            if self.Kcat:
                return rate_2_2(c_IPP, c_DMAPP, c_GPP, c_PPi, self.Kcat*self.c, self.Keq, self.IPP, self.DMAPP, self.GPP, self.PPi)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_2_2(c_IPP, c_DMAPP, c_GPP, c_PPi, U, self.Keq, self.IPP, self.DMAPP, self.GPP, self.PPi)

class LimSyn(Enzyme):

    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 5
        self.Kcat = 7e5 #4.827e5
        self.M = 72.4 #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC59327/#:~:text=The%20limonene%20synthase%20cDNA%20encodes,et%20al.%2C%201993).
        self.Keq = 5000
        self.GPP = 0.05 #0.002

        self.LIM = 2
        self.PPi = 0.5

        self.conc_U_mg = 0.01
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed
        self.reversible = False

    # @non_negative_return
    def get_rate(self, c_GPP, c_LIM, c_PPi):
        if not self.reversible:
            c_LIM = 0
            c_PPi = 0
        if self.packedBed:
            if self.Kcat:
                return rate_1_2(c_GPP, c_LIM, c_PPi, self.Kcat*self.c, self.Keq, self.GPP, self.LIM, self.PPi)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_1_2(c_GPP, c_LIM, c_PPi, U, self.Keq, self.GPP, self.LIM, self.PPi)

class PPase(Enzyme):

    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 20
        self.Kcat = 2206.3524049873267 #3828.2 #2.664e6
        self.M = 33.3 #https://www.genscript.com/enzyme/E00071-Inorganic_Pyrophosphatase_yeast_.html
        self.Keq = 10
        self.PPi = 0.1

        self.PI = 0.1

        self.conc_U_mg = 4.4
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed
        self.reversible = False

    # @non_negative_return
    def get_rate(self, c_PPi, c_PI):
        if not self.reversible:
            c_PI = 0
        if self.packedBed:
            if self.Kcat:
                return rate_1_2(c_PPi, c_PI, c_PI, self.Kcat*self.c, self.Keq, self.PPi, self.PI, self.PI)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_1_2(c_PPi, c_PI, c_PI, U, self.Keq, self.PPi, self.PI, self.PI)

class NoxE(Enzyme):
    def __init__(self, c=None, packedBed=True):
        super().__init__(c, packedBed)
        self.U = 1
        self.Kcat = 1136.3479611885175 #3257.3 #4.32e5
        self.M = 27 #https://pubmed.ncbi.nlm.nih.gov/8867896/#:~:text=S.,35%20kDa%20per%20subunit%3B%20S.
        self.NADH = 0.1
        self.conc_U_mg = 43.4
        if c:
            self.K = self.U/(c/self.M)
        else:
            self.K = None
        self.packedBed = packedBed

    # @non_negative_return
    def get_rate(self, c_NADH):
        if self.packedBed:
            if self.Kcat:
                return rate_1_H(c_NADH, self.Kcat*self.c, self.NADH)
            else:
                raise ValueError("Set K first by passing in concentration of enzyme")
        else:
            if self.c:
                U = self.Kcat*self.c
            else:
                U = self.U
            return rate_1_H(c_NADH, U, self.NADH)


if __name__=="__main__":
    # hex = Hex(packedBed=False)
    # print(hex.conc_U_mg)
    # # hex.set_enzyme_concentration_for_batch(0.001, 800)
    # hex.set_enzyme_conc_with_U(0.1, 800)
    # print(hex.c)
    # print(hex.c*3.3e6)
    l_enzyme_cls = [Hex, Pgi, Pfk, Fba, Tpi, Gap, mGap, Pgk, Pgm, Eno, Pyk, Pdh, PhaA, Hmgs, Hmgr, Mvk, Pmvk, Mdc, Idi, Fpps_S82F, LimSyn, PPase, NoxE]
    print(f"enzymes: {len(l_enzyme_cls)}")
    l_mg = [0.001,0.03,0.02,0.04,0.004,0.007,0.03,0.01,0.05,0.05,0.01,0.01,0.02,0.03,0.08,0.02,0.04,0.08,0.07,0.02,0.6,0.2,0.05]
    l_U = [0.01,0.6,0.2,0.9,0.4,0.1,0.4,1.1,4.6,1.4,2.9,0.2,1.3,0.04,0.4,0.2,0.6,0.3,0.3,0.1,0.006,0.1,2.1]
    enzyme = []
    conc = []
    for e, mg, U in zip(l_enzyme_cls, l_mg, l_U):
        en = e(packedBed=False)
        en.set_enzyme_concentration_with_mg(mg, 800)
        # en.set_enzyme_conc_with_U(U, 800)
        enzyme.append(en)
        conc.append(en.c)
        print(f"Enzyme {en.name} has U computed from Kcat[E]={en.c*en.Kcat:0.2f} compared with table data: {en.U}; k_cat should set to: {en.U/en.c} ")
    
    print(f"conc: {conc}")
    c = np.array(conc)
    np.savetxt("concentration_enzyme.txt", c)
