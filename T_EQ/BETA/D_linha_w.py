import numpy as np

def calculate_Dw_prime(T_gota_em_k, R,P, r_i, M_H2O , T0 , P0 , alpha_c = 0.036 , Delta_w = 8e-8):
    D_w = 2.11e-5 * (T_gota_em_k / T0)**1.94 * (P0 / P)
    denominator = r_i / (r_i + Delta_w) + (D_w / (r_i * alpha_c)) * ((2 * np.pi * M_H2O)/(R*T_gota_em_k))**0.5
    D_linha_water = D_w / denominator
    return D_linha_water

