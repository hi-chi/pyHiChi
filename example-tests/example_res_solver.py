from math import *
import numpy as np
import matplotlib.pyplot as plt 
import sys
sys.path.append("C:/Codes/hi-chi-new-module/pyHiChi/build/visual_studio/src/pyHiChi/Release") #change the path accordingly
from pyHiChi import *
import pyHiChi as hichi

def example_gonoskov_pre_2011_fig3():
    gridSize = 1024
    wavelength = 1e-4
    L = 3 * wavelength
    L_plasma = wavelength
    Solver = hichi.res_solver(gridSize, gridSize, gridSize, L, L, L_plasma)
    Ey_in = [0]*gridSize
    Ez_in = [0]*gridSize
    Ey_out = [0]*gridSize
    Ez_out = [0]*gridSize
    N_plasma = [0]*gridSize
    X = [0]*gridSize
    X_plasma = [0]*gridSize

    E_rel = 2 * pi * hichi.ELECTRON_MASS * (hichi.LIGHT_VELOCITY**2) / (-hichi.ELECTRON_CHARGE * wavelength)
    N_cr = pi * hichi.ELECTRON_MASS * (hichi.LIGHT_VELOCITY**2) / ((hichi.ELECTRON_CHARGE * wavelength) ** 2)
    print("E_rel = ", E_rel, ", N_cr = ", N_cr)
    for ind in range(gridSize):
        X[ind] = (-L + L*(ind/gridSize))/wavelength
        Ey_in[ind] = 191 * E_rel * sin(2 * pi * X[ind])
        Ez_in[ind] = 0
        Solver.E_in_set(ind, Ey_in[ind], Ez_in[ind])
        X_plasma[ind] = (L_plasma * (ind + 0.5)/gridSize)/wavelength
        N_plasma[ind] = 360 * N_cr
        Solver.N_plasma_set(ind, N_plasma[ind])
    Solver.incidence_angle_set(pi * 11.25 / 180.0)
    Solver.compute()
    for ind in range(gridSize):
        Ey_out[ind] = Solver.Ey_out_get(ind)
        Ez_out[ind] = Solver.Ez_out_get(ind)
    fig, ax = plt.subplots(2,1)
    im0 = ax[0].plot(X, Ey_in)
    im1 = ax[0].plot(X, Ez_in)
    im2 = ax[1].plot(X, Ey_out)
    im3 = ax[1].plot(X, Ez_out)
    ax[0].set_ylim([-3 * 191 * E_rel, 3 * 191 * E_rel])
    ax[1].set_ylim([-3 * 191 * E_rel, 3 * 191 * E_rel])
    plt.show()

def example_gonoskov_pop_2018_fig2():
    gridSize = 1024
    wavelength = 1e-4
    L_in = 3 * wavelength
    L_out = 3 * wavelength
    L_plasma = wavelength
    Solver = hichi.res_solver(gridSize, gridSize, gridSize, L_in, L_out, L_plasma)
    Ey_in = [0]*gridSize
    Ez_in = [0]*gridSize
    Ey_out = [0]*gridSize
    Ez_out = [0]*gridSize
    N_plasma = [0]*gridSize
    X_in = [0]*gridSize
    X_out = [0]*gridSize
    X_plasma = [0]*gridSize
    E_rel = 2 * pi * hichi.ELECTRON_MASS * (hichi.LIGHT_VELOCITY**2) / (-hichi.ELECTRON_CHARGE * wavelength)
    N_cr = pi * hichi.ELECTRON_MASS * (hichi.LIGHT_VELOCITY**2) / ((hichi.ELECTRON_CHARGE * wavelength) ** 2)
    print("E_rel = ", E_rel, ", N_cr = ", N_cr)
    for ind in range(gridSize):
        X_in[ind] = (-L_in + L_in*(ind/gridSize))/wavelength
        Ey_in[ind] = 300 * E_rel * sin(2 * pi * X_in[ind])
        Ez_in[ind] = 150 * E_rel * sin((7/4) * 2 * pi * X_in[ind])
        if(X_in[ind] < -2):
            Ey_in[ind] = 0
            Ez_in[ind] = 0
        Solver.E_in_set(ind, Ey_in[ind], Ez_in[ind])
        X_plasma[ind] = (L_plasma * (ind + 0.5)/gridSize)/wavelength # 0.5 is added to make > 0 everywhere
        N_plasma[ind] = 500 * N_cr
        if(X_plasma[ind] < 1/3.0):
            N_plasma[ind] = 500 * N_cr * X_plasma[ind]/ (1/3.0)
        if(X_plasma[ind] > 2/3.0):
            N_plasma[ind] = 500 * N_cr * (1 - (X_plasma[ind] - 2/3.0) / (1/3.0))
        Solver.N_plasma_set(ind, N_plasma[ind])
    Solver.incidence_angle_set(pi/7.0)
    Solver.compute()
    for ind in range(gridSize):
        X_out[ind] = (-L_out + L_out*(ind/gridSize))/wavelength
        Ey_out[ind] = Solver.Ey_out_get(ind)
        Ez_out[ind] = Solver.Ez_out_get(ind)
    fig, ax = plt.subplots(3,1)
    im0 = ax[0].plot(X_in, Ey_in)
    im1 = ax[0].plot(X_in, Ez_in)
    imn = ax[2].plot(X_plasma, N_plasma)
    im2 = ax[1].plot(X_out, Ey_out)
    im3 = ax[1].plot(X_out, Ez_out)
    plt.show()

example_gonoskov_pre_2011_fig3()
#example_gonoskov_pop_2018_fig2()