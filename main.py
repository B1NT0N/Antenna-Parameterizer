import numpy as np
import scipy
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math
import streamlit as st
from millify import millify

def FLD(x):
    theta = np.linspace(0, cmath.pi, 100)
    phi = np.linspace(0, 2*cmath.pi, 100)
    I_0 = 1
    eta = 120*cmath.pi
    C = 0.5772
    v = 3 * 10**8
    f = 10 * 10**6
    wl = v/f
    k = (2*cmath.pi)/wl
    l = x
    k=(2*cmath.pi)/wl


    F_theta = ((np.cos(((k*l)/2)*np.cos(theta))-np.cos(((k*l)/2)))/np.sin(theta))**2
    Q = (C+cmath.log(k*l)-scipy.special.sici(k*l)[1]+0.5*np.sin(k*l)*(scipy.special.sici(2*k*l)[0]-2*scipy.special.sici(k*l)[0])+0.5*np.cos(k*l)*(C+cmath.log(k*l/2)+scipy.special.sici(2*k*l)[1]-2*scipy.special.sici(k*l)[1]))

    #Radiation Intensity
    U = eta * ((abs(I_0)**2)/(8*cmath.pi**2))*F_theta

    #Radiation Power
    P_rad = eta * ((abs(I_0)**2)/(4*cmath.pi))*Q

    #Radiation Resistance
    R_r = (2*P_rad)/(abs(I_0)**2)

    #Maximum Directivity
    D_0 = 2*np.nanmax(F_theta)/Q

    #Maximum Effective Area
    A_em = ((wl**2)/(4*cmath.pi))*D_0

    #Input Radiation Resistance
    R_in = R_r/(np.sin(k*l/2)**2)

    new_theta, new_phi = np.meshgrid(theta, phi)
    new_U_n = eta * ((abs(I_0)**2)/(8*cmath.pi**2))*((np.cos(((k*l)/2)*np.cos(new_theta))-np.cos(((k*l)/2)))/np.sin(new_theta))**2

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
    new2_U_n = eta * ((abs(I_0)**2)/(8*cmath.pi**2))*((np.cos(((k*l)/2)*np.cos(new2_theta))-np.cos(((k*l)/2)))/np.sin(new2_theta))**2
    ax = fig.add_subplot(2, 2, 2, projection='polar')
    ax.plot(phi, np.ones(phi.size)*new2_U_n)

    new2_phi = math.pi;
    new3_U_n = eta * ((abs(I_0)**2)/(8*cmath.pi**2))*((np.cos(((k*l)/2)*np.cos(new_theta))-np.cos(((k*l)/2)))/np.sin(new_theta))**2
    ax = fig.add_subplot(2, 2, 3, projection='polar')
    ax.plot(theta, new3_U_n)

    new2_phi = 0;
    new3_U_n = eta * ((abs(I_0)**2)/(8*cmath.pi**2))*((np.cos(((k*l)/2)*np.cos(new_theta))-np.cos(((k*l)/2)))/np.sin(new_theta))**2
    ax = fig.add_subplot(2, 2, 4, projection='polar')
    ax.plot(theta, new3_U_n)

    return D_0, P_rad, A_em, R_r, R_in, fig

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
    l = x
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

def SD(x):

    theta = np.linspace(0, cmath.pi, 100)
    phi = np.linspace(0, 2*cmath.pi, 100)
    I_0 = 1
    eta = 120*cmath.pi
    C = 0.5772
    v = 3 * 10**8
    f = 10 * 10**6
    wl = v/f
    k = (2*cmath.pi)/wl
    l = x
    k=(2*cmath.pi)/wl

    r = (k*l)/2

    C1 = 1j* eta*((k*I_0*l*cmath.exp(-1j*k*r))/(8*cmath.pi*r))

    R_r = (20*cmath.pi**2)*(l/wl)**2

    R_in = R_r

    D_0 = 3/2

    A_em = ((3*wl**2)/(8*cmath.pi))

    new_theta, new_phi = np.meshgrid(theta, phi)
    new_U_n = np.sin(new_theta)

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
    new2_U_n = C1*np.sin(new2_theta)
    ax = fig.add_subplot(2, 2, 2, projection='polar')
    ax.plot(phi, np.ones(phi.size)*new2_U_n)


    new2_phi = math.pi;
    new3_U_n = C1*np.sin(new_theta)
    ax = fig.add_subplot(2, 2, 3, projection='polar')
    ax.plot(theta, new3_U_n)

    new2_phi = 0;
    new3_U_n = C1*np.sin(new_theta)
    ax = fig.add_subplot(2, 2, 4, projection='polar')
    ax.plot(theta, new3_U_n)

    return D_0, 0, A_em, R_r, R_in, fig

st.set_page_config(
    page_title='Antenna Parameterizer', 
    page_icon="https://cdn-icons-png.flaticon.com/512/1025/1025167.png",
    layout="wide"
    )
    
st.title('Antenna Parameterizer')
st.write('')

col1,col2 = st.columns([0.3,1])

with col1:
    x = st.slider(label='Choose a value for $l$', min_value=0.0, max_value=4.0, value=0.5, step=0.001,)

    v = 3 * 10**8
    f = 10 * 10**6
    wl = v/f
    k = (2*cmath.pi)/wl
    l = x


    if l>wl/50 and l<wl/10:
        st.info("SMALL DIPOLE",)
    elif l < wl/50:
        st.info("INFINITESIMAL DIPOLE",)
    elif l > wl/10:
        st.info("FINITE LENGTH DIPOLE",)
        
   
with col2:
    v = 3 * 10**8
    f = 10 * 10**6
    wl = v/f
    k = (2*cmath.pi)/wl
    l = x

    if l>wl/50 and l<wl/10:
        D_0, P_rad, A_em, R_r, R_in, fig = SD(x)
    elif l < wl/50:
        D_0, P_rad, A_em, R_r, R_in, fig = ID(x)
    elif l>wl/10:
        D_0, P_rad, A_em, R_r, R_in, fig = FLD(x)
    
    print(wl)
    col1, col2, col3, col4, col5 = st.columns(5)
    col1.metric("$D_0$", f"{millify(D_0, precision=2)}")
    col2.metric("$A_{em}$", f"{millify(A_em, precision=2)} m2")
    col3.metric("$R_r$", f"{millify(R_r, precision=2)} Ω")
    col4.metric("$R_{in}$", f"{millify(R_in, precision=2)} Ω")
    col5.metric("$P_{rad}$", f"{millify(P_rad, precision=2)} W")
    st.pyplot(fig=fig)

