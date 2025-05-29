"""
Author: Yudong Li (Yudong.Li@nrel.gov)
"""
from g2t.utils.returndecs import non_negative_return


#@non_negative_return
def rate_2_2(A, B, P, Q, Vmax, Keq, Kf_A, Kf_B, Kr_P, Kr_Q):
    denominator = ((1+A/Kf_A)*(1+B/Kf_B) + (1+P/Kr_P)*(1+Q/Kr_Q) -1)
    return Vmax * ((A*B - P*Q/Keq)/(Kf_A*Kf_B)/denominator)

#@non_negative_return
def rate_1_1(A, P, Vmax, Keq, Kf_A, Kr_P):
    denominator = (1+A/Kf_A + 1+P/Kr_P - 1)
    return Vmax*((A-P/Keq)/Kf_A)/(1+A/Kf_A + 1+P/Kr_P - 1)

#@non_negative_return
def rate_1_2(A, P, Q, Vmax, Keq, Kf_A, Kr_P, Kr_Q):
    denominator = ((1+A/Kf_A) + (1+P/Kr_P)*(1+Q/Kr_Q) -1)
    return Vmax * ((A - P*Q/Keq)/(Kf_A)/denominator)

#@non_negative_return
# def rate_3_2(A, B, C, P, Q, Vmax, Keq, Kf_A, Kf_B, Kf_C, Kr_P, Kr_Q):
#     return Vmax * ((A*B *C- P*Q/Keq)/(Kf_A*Kf_B*Kf_C)/((1+A/Kf_A)*(1+B/Kf_B)*(1+C/Kf_C) + (1+P/Kr_P)*(1+Q/Kr_Q) -1))

def rate_3_2(A, B, C, P, Q, Vmax, Keq, Kf_A, Kf_B, Kf_C, Kr_P, Kr_Q):
    return Vmax * ((A*B *C- P*Q/Keq)/(Kf_A*Kf_B*Kf_C)/((1+A/Kf_A)*(1+B/Kf_B)*(1+C/Kf_C) + (1+P/Kr_P)*(1+Q/Kr_Q) -1))

#@non_negative_return
def rate_3_3(A, B, C, P, Q, Z, Vmax, Keq, Kf_A, Kf_B, Kf_C, Kr_P, Kr_Q, Kr_Z):
    return Vmax * ((A*B *C- P*Q*Z/Keq)/(Kf_A*Kf_B*Kf_C)/((1+A/Kf_A)*(1+B/Kf_B)*(1+C/Kf_C) + (1+P/Kr_P)*(1+Q/Kr_Q)*(1+Z/Kr_Z) -1))

#@non_negative_return
def rate_3_4(A, B, C, P, Q, Z, X, Vmax, Keq, Kf_A, Kf_B, Kf_C, Kr_P, Kr_Q, Kr_Z, Kr_X):
    return Vmax * ((A*B *C- P*Q*Z*X/Keq)/(Kf_A*Kf_B*Kf_C)/((1+A/Kf_A)*(1+B/Kf_B)*(1+C/Kf_C) + (1+P/Kr_P)*(1+Q/Kr_Q)*(1+Z/Kr_Z)*(1+X/Kr_X)-1))

#@non_negative_return
def rate_2_4(A, B, P, Q, Z, X, Vmax, Keq, Kf_A, Kf_B, Kr_P, Kr_Q, Kr_Z, Kr_X):
    return Vmax * ((A*B-P*Q*Z*X/Keq)/(Kf_A*Kf_B)/((1+A/Kf_A)*(1+B/Kf_B)+(1+P/Kr_P)*(1+Q/Kr_Q)*(1+Z/Kr_Z)*(1+X/Kr_X) -1))

#@non_negative_return
def rate_1_H(A, Vmax, Kf_A):
    return Vmax*A/(Kf_A + A)

#@non_negative_return
def rate_mvk(A, B, P, Q, Vmax, Keq, Kf_A, Kf_B, Kr_P, Kr_Q):
    # return Vmax * ((A*B - P*Q/Keq)/(Kf_A*Kf_B)/((1+A/Kf_A+P/Kr_P)*(1+B/Kf_B+Q/Kr_Q)))
    return rate_2_2(A, B, P, Q, Vmax, Keq, Kf_A, Kf_B, Kr_P, Kr_Q)
