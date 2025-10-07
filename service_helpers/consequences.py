import math


def S(p, vi):
    return p * vi * vi


def dynamic_pressure(p0h, vi):
    return p0h / 2 * vi * vi

def fragmentation_density(S, vi):
    return 2 * S / (vi * vi)

def fragmentation_height(p):
    return -7,5 * math.log(p/1,225)

def kinetic_energy(m, v):
    return m * v * v / 2

def amplitude_magnitude(e):
    return math.log10(e) - 5.87

def tnt_equivalent(e):
    return e / (4,184 * 10e15)

def transition_d(d, pi, pt, v, g, o):
    return 1.161 * ((pi / pt) **  1/3) * (d ** 0,78) * (v ** 0,44) * (g ** -0.22) * (math.sin(o) ** 1/3)

def main_d(d):
    return d * 1,25

def transition_depth(k, d): # k = 0.2 or 0.15
    return k * d

def main_depth(k, d): # k = 0.15 or 0.1
    return k * d

def blast_effects(d1, w):
    return d1 * w  (1/3)

def initial_wave_height(Dtr, h_sea):
    return min(0.14 * Dtr, h_sea)

def wave_amplitude(A0, Dtr, l):
    return A0 * ((Dtr / l) ** 1.5)

def wave_length(Dtr, h_sea):
    return 2 * min(0.14 * Dtr, h_sea)

def wave_coefficient(s, omega, J):
    return s * (omega / (2 * J))

def wave_height_on_shore(k, J, xi):
    return k * xi * J
