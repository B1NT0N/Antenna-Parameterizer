import numpy as np
import scipy
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math
import streamlit as st
from millify import millify

def ID(x):
    theta = np.linspace(0, cmath.pi, 100)
    phi = np.linspace(0, 2*cmath.pi, 100)
    I_0 = 1
    eta = 120*cmath.pi
    C = 0.5772
    v = 3 * 10**8
    f = 10 * 10**6
    wl = v/f
    k = (2*cmath.pi)/wl
    l = x*wl
    k=(2*cmath.pi)/wl

    U = ((eta/2)*((k*I_0*l)/(4*cmath.pi))**2)*np.sin(theta)

    D_0 = 3/2

    A_em = ((3*wl**2)/(8*cmath.pi))

    R_r = (80*cmath.pi**2)*(l/wl)**2

    R_in = R_r

    P_rad = 0.5*(abs(I_0)**2)*R_r

    new_theta, new_phi = np.meshgrid(theta, phi)
    new_U_n = (eta/2)*((k*I_0*l)/(4*cmath.pi))**2*np.sin(new_theta)

    x = new_U_n*np.sin(new_theta)*np.cos(new_phi)
    y = new_U_n*np.sin(new_theta)*np.sin(new_phi)
    z = new_U_n*np.cos(new_theta)

    # Plot of the 3D Mesh
    fig = plt.figure()
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    ax.plot_surface(x, y, z, cmap="magma")

    #Diagrama de Radiacao 2D

    #Plano XY
    new2_theta = math.pi/2
    new2_U_n = ((eta/2)*((k*I_0*l)/(4*cmath.pi))**2)*np.sin(new2_theta)
    ax = fig.add_subplot(2, 2, 2, projection='polar')
    ax.plot(phi, np.ones(phi.size)*new2_U_n)


    new2_phi = math.pi;
    new3_U_n = ((eta/2)*((k*I_0*l)/(4*cmath.pi))**2)*np.sin(theta)
    ax = fig.add_subplot(2, 2, 3, projection='polar')
    ax.plot(theta, new3_U_n)

    new2_phi = 0;
    new3_U_n = ((eta/2)*((k*I_0*l)/(4*cmath.pi))**2)*np.sin(theta)
    ax = fig.add_subplot(2, 2, 4, projection='polar')
    ax.plot(theta, new3_U_n)

    return D_0, P_rad, A_em, R_r, R_in, fig



print(ID(0.02))