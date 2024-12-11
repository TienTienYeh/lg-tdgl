# import logging
# import os
# import shutil
# from contextlib import nullcontext
# from logging import Logger
from typing import Any, Dict, Literal, Optional, Sequence, Tuple, Union

# from matplotlib import animation
from matplotlib import pyplot as plt
# from tqdm import tqdm

# from tdgl.device.device import Device
# from tdgl.solution.data import get_data_range
# from tdgl.visualization.common import DEFAULT_QUANTITIES, PLOT_DEFAULTS, Quantity, auto_grid
# from tdgl.visualization.io import get_plot_data, get_state_string

### LG beam

import numpy as np
from scipy.special import genlaguerre
import scipy.constants as cons
from tdgl.parameter import Parameter
import pint # https://pint.readthedocs.io/en/0.10.1/tutorial.html
from pint import UnitRegistry
ureg = UnitRegistry()
# %pip install pint==0.23

def uniform_Bz_vector_potential(
    positions: np.ndarray,
    Bz: Union[float, str, pint.Quantity],
) -> np.ndarray:
    """Calculates the magnetic vector potential [Ax, Ay, Az] at ``positions``
    due uniform magnetic field along the z-axis with strength ``Bz``.

    Args:
        positions: Shape (n, 3) array of (x, y, z) positions in meters at which to
            evaluate the vector potential.
        Bz: The strength of the uniform field, as a pint-parseable string,
            a pint.Quantity, or a float with units of Tesla.

    Returns:
        Shape (n, 3) array of the vector potential [Ax, Ay, Az] at ``positions``
        in units of Tesla * meter.
    """
    positions = np.atleast_2d(positions)
    xs = positions[:, 0]
    ys = positions[:, 1]
    dx = np.ptp(xs)
    dy = np.ptp(ys)
    xs = xs - (xs.min() + dx / 2)
    ys = ys - (ys.min() + dy / 2)
    Ax = -Bz * ys / 2
    Ay = Bz * xs / 2
    A = np.stack([Ax, Ay, np.zeros_like(Ax)], axis=1)
    return A

def constant_field_vector_potential(
    x,
    y,
    z,
    *,
    Bz: float,
    field_units: str = "mT",
    length_units: str = "um",
):
    if z.ndim == 0:
        z = z * np.ones_like(x)
    positions = np.array([x.squeeze(), y.squeeze(), z.squeeze()]).T
    # Bz = Bz * ureg(field_units)
    A = uniform_Bz_vector_potential(positions, Bz)
    return A

def findval(X,x_value):
    return np.argmin(abs(X-x_value))

### ------------------------------------------------------------------------------------------ ###

### A of LG beam (equip the functions of "Step structured Bz" and "constant uniform Bz")

def Cpl(p,l):
  if abs(p)==0: p=1
  factorial_p = 1
  factorial_pl = 1
  for i in range(abs(int(p))): factorial_p = factorial_p*(i+1)
  for i in range(abs(int(p))+abs(int(l))): factorial_pl = factorial_pl*(i+1)
  return np.sqrt(2*factorial_p/np.pi/factorial_pl)

def A_LG_t(x, y, z, *, t,
                   w: float = 1.0,
                   E0: float = 1.0,
                   w0: float = 1.0,
                   xc_Gauss: float = 0.0,
                   yc_Gauss: float = 0.0,
                   z0: float = 0.0,
                   n: float = 1.0,
                   phi0_t: float = 0,
                   phi0_xy: float = 0,
                   tau: float = 1.0,
                   p: int = 0.0,
                   l: int = 0.0,
                   s: float = 0.0,
                   c: float = 1.0,
                   t_on: float = 0.0,
                   t_off: float = 1.0,
                   Bz: float = 0,
                   polarization_modulation: bool = False,
                   polarization: str = 'none',
                   angular_freq_units: str = "THz",
                   length_units: str = "um",
                   field_units: str = "mT",
                   E_field_units: str = "newton / coulomb",
                   time_units: str = "ps",
                   take_complex=False,
                   time_evolute: bool = True,

):
    """ Vector potential of Laguerre-Gaussian beam of p-th degree of mode and l-th order of mode
    From E to A, A = -iE/w
    """
    X = (x-xc_Gauss) * 1 # ureg(length_units).to("m").magnitude   # mesh of x-axis [um]
    Y = (y-yc_Gauss) * 1 # ureg(length_units).to("m").magnitude   # mesh of y-axis [um]
    Z = (z-z0) * 1 # ureg(length_units)).to("m").magnitude   # mesh of y-axis [um]
    w0 = w0 * 1 # ureg(length_units) .to("m").magnitude   # Beam waist (suppose w0_x =w0_y) [um]
    k = w/c   # Wavenumber
    ti = np.copy(t)#(t * ureg(time_units)).to("s").magnitude   # Input time
    # E0 = (E0 * ureg(E_field_units)).to("newton / coulomb").magnitude   # Electrical field |E0|
    zR = n*k*w0**2/2   # Rayleigh length, n: refractive index
    wz = w0 * np.sqrt(1+(Z/zR)**2)   # Spot radius
    zeta = Z/zR

    r = np.sqrt(X**2+Y**2)   # Radius
    phi = np.angle(X+1j*Y)   # Azimuthal angle
    if polarization.lower()=='x' or polarization.lower()=='linear x': s, phi0_xy = [0.0,0.0]
    if polarization.lower()=='y' or polarization.lower()=='linear y': s, phi0_xy = [0.0,np.pi/2]
    if polarization.lower()=='lc' or polarization.lower()=='left circular': s, phi0_xy = [1,np.pi/4]
    if polarization.lower()=='rc' or polarization.lower()=='right circular': s, phi0_xy = [-1,np.pi/4]
    phi_t = np.copy(phi0_t) # ONLY CONSTANT PHASE REMAINS if time_evolute==False
    if time_evolute: phi_t = phi_t + w*ti
    phiGouy = (2*p+np.abs(l)+1)*np.arctan(Z/zR) # Gouy phase
    u = E0 * Cpl(p,l)/w0 * (np.sqrt(2)*r/wz)**l * genlaguerre(p,l)(2*r**2/wz**2) * w0/wz * np.exp(-1j*phiGouy +1j*(l*phi +k*z +k*r**2/2/(z-1j*zR)))


    if t>t_off or t<t_on: t_step = 0
    else: t_step = 1

    pol_m_x, pol_m_y = [np.cos(phi0_xy),np.sin(phi0_xy)]
    if polarization_modulation: pol_m_x, pol_m_y = (np.abs([np.cos(l*phi+phi0_xy), np.sin(l*phi+phi0_xy)]))
    Ex = u * np.exp(-1j*(phi_t)) * pol_m_x * t_step
    Ey = u * np.exp(-1j*(phi_t+s*np.pi/2)) * pol_m_y * t_step

    Ax = -1j/w*Ex
    Ay = -1j/w*Ey
    Az = np.zeros_like(Ax)
    Ax[np.isnan(Ax)] = 0
    Ax[np.isinf(Ax)] = 0
    Ay[np.isnan(Ay)] = 0
    Ay[np.isinf(Ay)] = 0
    A_constBz = constant_field_vector_potential(x, y, z, Bz=Bz, field_units=field_units, length_units=length_units)
    A = np.stack([np.real(Ax), np.real(Ay), np.real(Az)], axis=1) + A_constBz
    if take_complex: A = np.stack([(Ax), (Ay), (Az)], axis=1) + A_constBz
    # A = np.array([np.real(Ax), np.real(Ay), np.real(Az)]).T
    return A#.to(f"{field_units} * {length_units}").magnitude

# callable(A_LG_t_xy)

def A_LG(*,
        w: float = 1.0,
        E0: float = 1.0,
        w0: float = 1.0,
        xc_Gauss: float = 0.0,
        yc_Gauss: float = 0.0,
        z0: float = 0.0,
        n: float = 1.0,
        phi0_t: float = 0,
        phi0_xy: float = 0,
        tau: float = 1.0,
        c: float = 1.0,
        p: int = 0.0,
        l: int = 0.0,
        s: float = 0.0,
        t_on: float = 0.0,
        t_off: float = 1.0,
        Bz: float = 0,
        polarization_modulation: bool = False,
        polarization: str = 'none',
        angular_freq_units: str = "THz",
        length_units: str = "um",
        E_field_units: str = "newton / coulomb",
        field_units: str = "mT",
        time_units: str = "ps",
        time_evolute: bool = True,
        time_dependent=True,
)-> Parameter:
    """Vector potential of Laguerre-Gaussian beam  LG(p)(l)
    # for linear polarization, LG00, phi0_xy could be any number and s=0
    # for Circular polarization, LG00, phi0_xy=pi/4 and s=+-1
    # for linear polarization, LG01, phi0_xy could be any number and s=0
    # for Circular polarization, LG01, phi0_xy=pi/4 and s=+-1
    # for Radial polarization, LG01, s=+1, phi0_xy=0, polarization_modulation = True
    # for Azimuthal polarization, LG01, s=+1, phi0_xy=np.pi/2, polarization_modulation = True

    Equip the function "Step structured Bz" and "constant uniform Bz"
    # Step structured Bz:
        # Step time: t_on, t_off, time_evolute = [t_on, t_off, False]
        # Continuous case: t_on, t_off, time_evolute = [0, solve_time, True]
    # constant Bz:
        # Bz = Bz

    Note of useful relation: f=w/2p, c=fL, k=2p/L, k=w/c, 1/w=1/kc=L/2pc

    Args:
        w: angular frequency ( w = 2 pi f ) ,
        k_dir: prapagation direction
        E0: amplitude of electrical field
        phi0_t: initial phase of time
        phi0_xy: initial angle of xy plane azimuthal angle
        tau: Unit time (SC dissipation time)
        p: Degree of mode
        l: Order of mode, or orbital angular momentum of LG beam
        s: spin angular momentum of LG beam
    Returns:
        A :class:`tdgl.Parameter` that produces a linear ramp.
    """
    return Parameter(
        A_LG_t,
        w=w,
        w0=w0,
        E0=E0,
        phi0_t=phi0_t,
        phi0_xy=phi0_xy,
        xc_Gauss=xc_Gauss, yc_Gauss=yc_Gauss,
        p=p, l=l, s=s, c=c, Bz=Bz, z0=z0, n=n,
        tau=tau, t_on=t_on, t_off=t_off,
        polarization=polarization,
        polarization_modulation=polarization_modulation,
        angular_freq_units=angular_freq_units,
        length_units=length_units,
        E_field_units=E_field_units,
        field_units=field_units,
        time_units=time_units,
        time_evolute=time_evolute,
        time_dependent=True,
    )

### ------------------------------------------------------------------------------------------ ###

### B and E component EM wave

# B = curl(A) = (dyAz - dzAy) ex + (dzAx - dxAz) ey + (dxAy - dyAx) ez
def A2B(x, y, z, A):

    ''' Calculate magnetic field B from vector potential A
    return B
    '''
    B.x = np.diff(A[:,2])/np.diff(y) - np.diff(A[:,1])/np.diff(z)
    B.y = np.diff(A[:,0])/np.diff(z) - np.diff(A[:,2])/np.diff(x)
    B.z = np.diff(A[:,1])/np.diff(x) - np.diff(A[:,0])/np.diff(y)
    return B

def E_LG_t(x, y, z, *, t,
                   w: float = 1.0,
                   E0: float = 1.0,
                   w0: float = 1.0,
                   xc_Gauss: float = 0.0,
                   yc_Gauss: float = 0.0,
                   z0: float = 0.0,
                   n: float = 1.0,
                   phi0_t: float = 0,
                   phi0_xy: float = 0,
                   tau: float = 1.0,
                   p: float = 0.0,
                   l: float = 0.0,
                   s: float = 0.0,
                   c: float = 1.0,
                   t_on: float = 0.0,
                   t_off: float = 1.0,
                   Bz: float = 0,
                   polarization_modulation: bool = False,
                   polarization: str = 'none',
                   angular_freq_units: str = "THz",
                   length_units: str = "um",
                   field_units: str = "mT",
                   E_field_units: str = "newton / coulomb",
                   time_units: str = "ps",
                   time_evolute: bool = True,

):
    """ Electric field x, y of Laguerre-Gaussian beam of p-th degree of mode and l-th order of mode
    From E to A, A = -iE/w
    """
    X = (x-xc_Gauss) * 1 # ureg(length_units).to("m").magnitude   # mesh of x-axis [um]
    Y = (y-yc_Gauss) * 1 # ureg(length_units).to("m").magnitude   # mesh of y-axis [um]
    Z = (z-z0) * 1 # ureg(length_units)).to("m").magnitude   # mesh of y-axis [um]
    w0 = w0 * 1 # ureg(length_units) .to("m").magnitude   # Beam waist (suppose w0_x =w0_y) [um]
    k = w/c   # Wavenumber
    ti = np.copy(t)#(t * ureg(time_units)).to("s").magnitude   # Input time
    # E0 = (E0 * ureg(E_field_units)).to("newton / coulomb").magnitude   # Electrical field |E0|
    zR = (n*k*w0**2/2)   # Rayleigh length, n: refractive index
    wz = w0 * np.sqrt(1+(Z/zR)**2)   # Spot radius
    zeta = Z/zR

    r = np.sqrt(X**2+Y**2)   # Radius
    phi = np.angle(X+1j*Y)   # Azimuthal angle
    if polarization.lower()=='x' or polarization.lower()=='linear x': s, phi0_xy = [0.0,0.0]
    if polarization.lower()=='y' or polarization.lower()=='linear y': s, phi0_xy = [0.0,np.pi/2]
    if polarization.lower()=='lc' or polarization.lower()=='left circular': s, phi0_xy = [1,np.pi/4]
    if polarization.lower()=='rc' or polarization.lower()=='right circular': s, phi0_xy = [-1,np.pi/4]
    phi_t = np.copy(phi0_t) # ONLY CONSTANT PHASE REMAINS if time_evolute==False
    if time_evolute: phi_t = phi_t + w*ti
    phiGouy = (2*p+np.abs(l)+1)*np.arctan(Z/zR) # Gouy phase
    u = E0 * w0/wz * np.exp(-r**2/wz**2 -1j*phiGouy +1j*(l*phi +k*z +k*r**2/2/(z-1j*zR)))

    if t>t_off or t<t_on: t_step = 0
    else: t_step = 1

    pol_m_x, pol_m_y = [np.cos(phi0_xy),np.sin(phi0_xy)]
    if polarization_modulation: pol_m_x, pol_m_y = (np.abs([np.cos(l*phi+phi0_xy), np.sin(l*phi+phi0_xy)]))
    Ex = u * np.exp(-1j*(phi_t)) * pol_m_x * t_step
    Ey = u * np.exp(-1j*(phi_t+s*np.pi/2)) * pol_m_y * t_step
    Ez = np.zeros_like(Ex)
    return Ex, Ey, Ez

def E2B(x,y,Ex,Ey,Bz_constant,c,w):
    By = Ex/c
    Bx = -Ey/c
    Ax = -1j/w*Ex
    Ay = -1j/w*Ey
    Bz_A = np.zeros_like(By)
    Bz_A[1:] = np.diff(Ay)/np.diff(x) - np.diff(Ax)/np.diff(y)
    Bz = Bz_A + Bz_constant
    # B = np.stack([np.real(Bx), np.real(By), np.real(Bz)], axis=1)
    return Bx, By, Bz #.to(f"{field_units}").magnitude

def E2Bv(xv,yv,Ex,Ey,Bz_constant,c,w):
    By = Ex/c
    Bx = -Ey/c
    Ax = -1j/w*Ex
    Ay = -1j/w*Ey
    Bz_A = np.zeros_like(By)
    dAydx = np.diff(np.real(Ay),axis=1)/np.diff(xv,axis=1)
    dAxdy = np.diff(np.real(Ax),axis=0)/np.diff(yv,axis=0)
    Bz_A[1:,1:] =  dAydx[1:,:] - dAxdy[:,1:]
    Bz = Bz_A + np.ones_like(Bz_A)*Bz_constant
    # B = np.stack([np.real(Bx), np.real(By), np.real(Bz)], axis=1)
    return Bx, By, Bz #.to(f"{field_units}").magnitude

def find_max_Bz(P):
    Ex, Ey = P.E_input_frame(0,take_real=False)
    Bx, By, Bz1 = E2Bv(P.Xv,P.Yv,Ex,P.E0i*Ey,0,P.c,P.w_GL)
    Ex, Ey = P.E_input_frame(2*np.pi/4/P.w_GL,take_real=False)
    Bx, By, Bz2 = E2Bv(P.Xv,P.Yv,P.E0i*Ex,P.E0i*Ey,0,P.c,P.w_GL)
    Ex, Ey = P.E_input_frame(2*np.pi/2/P.w_GL,take_real=False)
    Bx, By, Bz3 = E2Bv(P.Xv,P.Yv,P.E0i*Ex,P.E0i*Ey,0,P.c,P.w_GL)
    return max([abs(np.real(Bz1)).max(), abs(np.real(Bz2)).max(), abs(np.real(Bz3)).max()]) + abs(P.constant_Bz)


### ------------------------------------------------------------------------------------------ ###

### Other useful functions

def light_state_control_LaguerreGauss(*,keyword_of_state="None"):

    ''' Select the parameters for optical states
    options: 'lg00_l_x','lg00_l_y','lg00_c_l','lg00_c_r','lg01_l_x','lg24_l_y',etc.
    return p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head
    '''
    if keyword_of_state[:2].lower()!='lg': keyword_of_state = 'None'
    if keyword_of_state=='None': 
        output_file_head = 'LGbeam'
        return 0, 0, 0, 0, 0, False, 'Gaussian_beam_linear_x'
    else:
        p = int(keyword_of_state[2])
        l = int(keyword_of_state[3])
        if keyword_of_state[-1].lower()=='x': s, phi0_t, phi0_xy, polarization_modulation, output_file_head_suffix = [0,0,0, False,'_linear_x']
        if keyword_of_state[-1].lower()=='y': s, phi0_t, phi0_xy, polarization_modulation, output_file_head_suffix = [0,0,np.pi/2, False,'_linear_y']
        if keyword_of_state[-1].lower()=='l': s, phi0_t, phi0_xy, polarization_modulation, output_file_head_suffix = [1,0,np.pi/4, False,'_circular_l']
        if keyword_of_state[-1].lower()=='r': s, phi0_t, phi0_xy, polarization_modulation, output_file_head_suffix = [-1,0,np.pi/4, False,'_circular_r']
        output_file_head = 'LG'+str(p)+str(l)+output_file_head_suffix
        return p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head

class input_value:
    def __init__(self,*,
            # Unit:
                 length_units: str =  "um", # SI: m
                 time_units: str  = 'ps', # SI: s
                 current_units: str = "mA", # SI: A
                 mass_units: str = 'kg', # SI: kg
                 angular_freq_units: str = 'THz', # SI: 1/s # for demonstration
                 field_units: str = "mT", # SI: kg/s^2/A
                 E_field_units: str = "kvolt/meter", # newton per coulomb (N/C), or volt per meter (V/m), SI: kg*m^2/s^3/A
                 resistivity_units: str = 'ohm * cm',
            # [SC] Properties SC:
                xi: float = 100/1000, # Coherent length
                london_lambda: float = 100/1000, # London penetration depth
                gamma: float = 10, # Strength of inelastic scattering
                u: float = 5.79, # Time ratio of order parameter relaxation
            # [SC] Size of sample
                height: float = 2, # Height of SC sheet
                width: float = 2, # Width of SC sheet
                thickness: float = 2/1000, # Thickness of SC sheet
            # [SC] Conductivity of sample and others
                temperature: float = 0.5, # Temperature of SC (Unit of Tc)
                resistivity: float = 150e-6, # Resistivity of normal state
            # [EM] Properties of EM wave
                E_amp: float = 200, # Input amplitude of electric field (unit: A0*w_EM)
                w_EM: float = -1, # Angular frequency of light (choose one "w_EM" "f_EM" which >=0)
                f_EM: float = -1, # Frequency of light (choose one "w_EM" "f_EM" which >=0)
                w_0: float = 0.4,# Radius of spot of Gaussian beam
                xc_Gauss: float = 0.0, # Center position of Gaussian beam
                yc_Gauss: float = 0.0, # Center position of Gaussian beam
                light_source_type: str = 'None', # Setting of light source (autometically setting with keywords "linear_x", "linear_y", "x" , or "y")
                p: int = 0, # Quantum number of radial order: p
                l: int = 0, # Quantum number of orbital angular momentum: l
                s = 0, # Quantum number of spin angular momentum: s
            # [EM] Details of phase
                phi0_t: float = 0.0, # Phase shift of time, i.e. w*t +　phi0_t
                phi0_xy: float = 0.0, # Ａngle shift of polarization, e.g. 0:'x-pol', and pi/2:'y-pol'
                polarization_modulation: bool = False, # Polarization modulation with azimuthal angle around center
            # Setting of system
                solve_time: float = 5.0, # Total solving time (Unit: tau_GL)
                screenSet: bool = False, # Setting of screening effect
                constant_Bz: float = 0.0,
            # Settings of figures demonstration
                quiver_mesh_n: float = 20, # Plot of E (quiver): mesh number
                quiver_scale: float = 5, # Plot of E (quiver): quier scale
                width_quiver: float = 0.1, # Plot of E (quiver): quiver width
                dpi: int = 100, # DPI of figure
):

        ''' Function for manage the input parameters for LG_TDGL '''
        ureg = UnitRegistry()
    # Unit:
        self.length_units = length_units
        self.time_units = time_units
        self.current_units = current_units
        self.mass_units = mass_units
        self.angular_freq_units = angular_freq_units
        self.field_units = field_units
        self.E_field_units = E_field_units
        self.resistivity_units = resistivity_units
    ### Superconductor (SC) thin film ###
    # [SC] SC properties:
        self.xi = xi # Coherent length
        self.london_lambda = london_lambda # London penetration depth
        self.gamma = gamma
        self.u = u
    # [SC] Size of sample
        self.height = height
        self.width = width
        self.thickness = thickness
        self.pearl_lambda = london_lambda**2/thickness # Updated value
        xi_coherent = xi * ureg(length_units)
        lambdaL = london_lambda * ureg(length_units)
        d_thickness = thickness * ureg(length_units)
        lambdaL_eff = lambdaL**2/d_thickness
    # [SC] Size of sample
        self.temperature = temperature # Unit of Tc
        self.disorder_epsilon = 1 # Originally, the value in pytdgl is 1/temperature-1. Here we fixed it as 1 and change the unit as temperature dependent dimentions.
    # [SC] onductivity of sample
        self.resistivity = resistivity
        rn_resistivity = resistivity* ureg(resistivity_units)
        rho_conductivity = (1/rn_resistivity).to('1 / ohm / '+length_units)
        self.conductivity = rho_conductivity.to('1 / ('+resistivity_units+')').magnitude # Updated value
    # [SC] Unit of Time and c
        mu_0 = cons.mu_0 * ureg('newton/A**2')
        tau_0 = mu_0 * rho_conductivity * lambdaL**2
        tau = tau_0.to(time_units).magnitude # Unit time, tau_GL in the real physical unit
        speed_of_light = cons.speed_of_light * ureg('m/s')
        revised_speed_of_light = (speed_of_light.to(f"{length_units}/{time_units}") * tau).magnitude # Speed of light after unit transformation
        self.c = revised_speed_of_light # Updated value
        self.tau = tau # Updated value, tau_GL
        self.tau_GL = tau # Updated value, tau_GL
        self.w_GL = 2*np.pi/tau # Updated value, the real physical value of angular freq of tau_GL
        self.f_GL = 1/tau # Updated value, the real physical value of freq of tau_GL
    # [SC] Unit of critical fields and current (only for examing)
        self.Phi0 = cons.h/2/cons.e * ureg('J/A').to(field_units+'*'+length_units+'^2')
        self.Hc1 = self.Phi0/4/np.pi/mu_0/lambdaL**2*np.log(lambdaL/xi_coherent) # [ref] Gennes.P.D., pp.66, Eq (3-56), but no mu_0
        self.Hc2 = self.Phi0/2/np.pi/mu_0/xi_coherent**2 # [ref] Logan
        self.Hc  = self.Hc1/(np.pi/np.sqrt(24)*xi_coherent/lambdaL*np.log(lambdaL/xi_coherent)) # [ref] Gennes.P.D., pp.66, Eq (3-56), but no mu_0
        self.Bc1 = (mu_0*self.Hc1).to(field_units)
        self.Bc2 = (mu_0*self.Hc2).to(field_units)
        self.Bc  = (mu_0*self.Hc).to(field_units)
        self.A0 = (xi_coherent*self.Bc2)
        self.J0 = (4*xi_coherent*self.Bc2/mu_0/lambdaL**2).to(f"{current_units}/{length_units}^2")
    ### Applied electromagnetic (EM) wave ###
    # [EM] Frequency of light
        # Angular frequency of light & frequency of light
        w_EM_setting, f_EM_setting = [np.copy(w_EM), np.copy(f_EM)]
        if w_EM_setting<0:
            self.w_EM = f_EM*(2*np.pi)
            self.f_EM = f_EM
        elif f_EM_setting<0:
            self.w_EM = w_EM
            self.f_EM = w_EM/(2*np.pi)
        else:
            self.w_EM = w_EM
            self.f_EM = f_EM
        self.tau_EM = 2*np.pi/self.w_EM # Unit of 1/tau_GL
        self.f_EM_unit = self.f_EM * self.f_GL
        self.w_EM_unit = self.w_EM * self.f_GL
        self.tau_EM_unit = self.tau_EM * self.tau_GL
    # [EM] Fundamental properties
        B0_EM = (1/self.c) * ureg(field_units) # Output value to TDGL detedmined by unit
        A0_EM = (1/self.w_GL) * ureg(f"{field_units} * {length_units}") # Output value to TDGL detedmined by unit
        E0_EM = B0_EM * speed_of_light # Output value to TDGL detedmined by unit
        self.E0_GL = self.A0 * self.f_GL * ureg('1/'+time_units) # Output value to TDGL detedmined by unit
        self.E0 = (E_amp * self.E0_GL).to(E_field_units) # Input amplitude of electric field
        self.E_amp = E_amp
        # E0i = (self.E_amp/E0_EM.to(E_field_units)).magnitude # From E_amp with "E_field_units" to unitless TDGL
        self.E0i = E_amp # Updated value
        # * ureg("1/angular_freq_units").to("time_units").magnitude
        self.lambda_EM = (speed_of_light * self.tau_EM_unit * ureg(time_units).to("s")).to(length_units).magnitude
        self.w_0 = w_0 # Radius of spot of Gaussian beam
        self.xc_Gauss = xc_Gauss # Center position of Gaussian beam
        self.yc_Gauss = yc_Gauss # Center position of Gaussian beam
        self.light_source_type = light_source_type # Setting of light source (autometically setting with keywords "linear_x", "linear_y", "x" , or "y")
        self.B_EM_inplane = self.E0i*B0_EM
        self.A_EM_inplane = self.E0i*A0_EM
        self.phi0_t = phi0_t # Phase shift of time, i.e. w*t +　phi0_t
        self.phi0_xy = phi0_xy # Ａngle shift of polarization, e.g. 0:'x-pol', and pi/2:'y-pol'
        self.polarization_modulation = polarization_modulation
        self.t_on, self.t_off, self.time_evolute = [0, solve_time, 'True']
        self.p = p
        self.l = l
        self.s = s
        self.light_source_type = light_source_type
        keyword_list = ['_l_x','_l_y','_c_r','_c_l']
        if light_source_type[:2].lower()=='lg' and light_source_type[-4:].lower() in keyword_list:
            p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = light_state_control_LaguerreGauss(keyword_of_state=light_source_type)
            self.p = p
            self.l = l
            self.s = s
            self.phi0_t = phi0_t
            self.phi0_xy = phi0_xy
            self.polarization_modulation = polarization_modulation
            self.output_file_head = output_file_head
    # Settings of system
        self.solve_time = solve_time # Total solving time (Unit: tau_GL)
        self.screenSet = screenSet
        self.constant_Bz = constant_Bz
    # Settings of figures demonstration
        self.quiver_mesh_n = quiver_mesh_n # Plot of E (quiver): mesh number
        self.quiver_scale = quiver_scale # Plot of E (quiver): quier scale
        self.width_quiver = width_quiver # Plot of E (quiver): quiver width
        self.dpi = dpi # DPI of figure
    # Settings of mesh for figures
        X, Y, Xv, Yv, Zv = v_grid_generation(-width/2,width/2,-height/2,height/2,quiver_mesh_n)
        self.X = X
        self.Y = Y
        self.Xv = Xv
        self.Yv = Yv
        self.Zv = Zv

        def E_input_frame(self, ti,*,take_real: bool=True):
            Zv = np.zeros_like(self.Xv)
            Ex, Ey, Ez =  (E_LG_t(self.Xv, self.Yv, self.Zv, t=ti, w=self.w_EM, w0=self.w_0, E0=self.E0i, xc_Gauss=self.xc_Gauss, yc_Gauss=self.yc_Gauss,
                                     phi0_t=self.phi0_t, phi0_xy=self.phi0_xy, p=self.p, l=self.l, s=self.s,
                                     tau=self.tau, polarization_modulation=self.polarization_modulation,
                                     t_on=self.t_on, t_off=self.t_off, Bz=self.constant_Bz, time_evolute=self.time_evolute,
                                     angular_freq_units=self.angular_freq_units, length_units=self.length_units, E_field_units=self.E_field_units, time_units=self.time_units,))
            if take_real: return np.real(Ex)/self.E0i, np.real(Ey)/self.E0i
            else:         return Ex/self.E0i, Ey/self.E0i
        def find_max_Bz(self):
            Ex, Ey = E_input_frame(self,0,take_real=False)
            Bx, By, Bz1 = E2Bv(self.Xv,self.Yv,self.E0i*Ex,self.E0i*Ey,0,self.c,self.w_EM)
            Ex, Ey = E_input_frame(self,2*np.pi/4/self.w_EM,take_real=False)
            Bx, By, Bz2 = E2Bv(self.Xv,self.Yv,self.E0i*Ex,self.E0i*Ey,0,self.c,self.w_EM)
            Ex, Ey = E_input_frame(self,2*np.pi/2/self.w_EM,take_real=False)
            Bx, By, Bz3 = E2Bv(self.Xv,self.Yv,self.E0i*Ex,self.E0i*Ey,0,self.c,self.w_EM)
            return max([abs(np.real(Bz1)).max(), abs(np.real(Bz2)).max(), abs(np.real(Bz3)).max()]) + abs(self.constant_Bz)
        self.Bz_max = find_max_Bz(self)

    def print_properties(self,Kwargs_list):
        kwargs_list = Kwargs_list
        for i in range(len(Kwargs_list)): kwargs_list[i] = kwargs_list[i].lower()
        PrintAll = 'all' in kwargs_list
        if PrintAll: print("[1] Length scale of sample")
        if 'xi' in kwargs_list or PrintAll: print("Coherent length (xi): "+str(self.xi * ureg(self.length_units)))
        if 'london_lambda' in kwargs_list or PrintAll: print("London penetration depth (london_lambda): "+str(self.london_lambda * ureg(self.length_units)))
        if 'thickness' in kwargs_list or PrintAll: print("Thickness (thickness): "+str(self.thickness * ureg(self.length_units)))
        if 'pearl_lambda' in kwargs_list or PrintAll: print("Pearl penetration depth (lambdaL**2/thickness): "+str(self.pearl_lambda * ureg(self.length_units)))
        if 'kapa' in kwargs_list or PrintAll: print("Ratio of length kapa (lambdaL/xi): "+str(self.london_lambda/self.xi))
        if PrintAll: print(" ")
        if PrintAll: print("[2] Condictivity and of sample")
        if 'resistivity' in kwargs_list or PrintAll: print("Resistivity (resistivity): "+str(self.resistivity * ureg(self.resistivity_units)))
        if 'conductivity' in kwargs_list or PrintAll: print("Conductivity (1/resistivity): "+str(self.conductivity * ureg(f"1/({self.resistivity_units})")))
        if PrintAll: print(" ")
        if PrintAll: print("[3] Time scale of sample")
        if 'tau_gl' in kwargs_list or PrintAll: print("Characteristic timescale (tau_GL): "+str(self.tau * ureg(self.time_units))) # characteristic timescale for this TDGL model
        if 'f_gl' in kwargs_list or PrintAll: print("Characteristic rate (1/tau_GL): "+str((1/self.tau) * ureg('1/'+self.time_units).to(self.angular_freq_units)))
        # if 'w_GL' in kwargs_list or PrintAll: print('Characteristic rate: {!s}'.format(self.w_GL*ureg(angular_freq_units)))
        if 'c' in kwargs_list or PrintAll: print("c input into unitless system (c): "+str(self.c))
        if PrintAll: print(" ")
        if PrintAll: print("[4] Strength of inelastic scattering")
        if 'gamma' in kwargs_list or PrintAll: print("Strength of inelastic scattering (gamma=2*tau_eph*gap_0): "+str(self.gamma))
        if PrintAll: print(" ")
        if PrintAll: print("[5] Critical E and B of SC")
        if 'phi0' in kwargs_list or PrintAll: print("Quantum magnetic flux (Phi0 = h/2e): "+str(self.Phi0))
        if 'hc1' in kwargs_list or PrintAll: print("Hc1 of SC (Phi0/4pi/lambda^2*ln[lambda/xi]): "+str(self.Hc1)) # [ref] Gennes.P.D., pp.66, Eq (3-56)
        if 'hc2' in kwargs_list or PrintAll: print("Hc2 of SC (Phi0/2pi xi^2): "+str(self.Hc2))
        if 'hc' in kwargs_list or PrintAll: print("Hc of SC: "+str(self.Hc))
        if 'bc1' in kwargs_list or PrintAll: print("Bc1 of SC (mu0Hc1): "+str((self.Bc1))) # lower critical field
        if 'bc2' in kwargs_list or PrintAll: print("Bc2 of SC (mu0Hc2): "+str((self.Bc2))) # upper critical field
        if 'bc' in kwargs_list or PrintAll: print("Bc of SC: "+str(self.Bc))
        # if 'bc/phi0' in kwargs_list or PrintAll: print("Bc*pi*w0^2/Phi0: "+str(self.Bc*np.pi*(self.w_0*ureg(self.length_units))**2/self.Phi0))
        # if 'bc2/phi0' in kwargs_list or PrintAll: print("Bc2*pi*w0^2/Phi0: "+str(self.Bc2*np.pi*(self.w_0*ureg(self.length_units))**2/self.Phi0))
        if 'a0' in kwargs_list or PrintAll: print("Unit vector potential A0 (xi*Bc2): "+str(self.A0))
        if 'j0' in kwargs_list or PrintAll: print("Unit current density (J0 = 4*xi*Bc2/mu_0/lambdaL**2): "+str(self.J0))
        if PrintAll: print(" ")
        if PrintAll: print("[6] Parameters of light source")
        if 'w_0' in kwargs_list or PrintAll: print("Beam size (2w0): "+str(2*self.w_0 * ureg(self.length_units)))
        if 'lambda_em' in kwargs_list or PrintAll: print("Wavelength: "+str(self.lambda_EM * ureg(self.length_units)))
        if 'w_em' in kwargs_list or PrintAll: print("Angular frequency of light (w, unit of 1/tau_GL): "+str(self.w_EM))
        if 'w_em_unit' in kwargs_list or PrintAll: print("Angular frequency of light (w, value with unit): "+str((self.w_EM_unit) * ureg('1/'+self.time_units).to(self.angular_freq_units)))
        if 'f_em' in kwargs_list or PrintAll: print("Frequency of light (w/2pi, unit of tau_GL): "+str(self.f_EM))
        if 'f_em_unit' in kwargs_list or PrintAll: print("Frequency of light (w/2pi, real value with unit): "+str((self.f_EM_unit) * ureg('1/'+self.time_units).to(self.angular_freq_units)))
        if 'e0_gl' in kwargs_list or PrintAll: print("Unit of electric field E0_GL (A0*f_GL): "+str(self.E0_GL)) 
        if 'e_amp' in kwargs_list or PrintAll: print("|E0| of light: "+str(self.E0)) #  * ureg(self.E_field_units)
        if 'e0i' in kwargs_list or PrintAll: print("Value of dimentionless |E0|: "+str(self.E0i))
        if 'b_em_inplane' in kwargs_list or PrintAll: print("In-plane |B0| of light (|E0|/c): "+str(self.B_EM_inplane))
        if 'a_em_inplane' in kwargs_list or PrintAll: print("In-plane |A0| of light (|E0|/w): "+str(self.A_EM_inplane))
#         if 'tau_GL' or 'all' in kwargs_list: print("Check ratio of 2pi|A0|/|B0|: "+str(2*np.pi*A0/B0)+"\n")
        if 'bz_max' in kwargs_list or PrintAll: print("|Bz| of light: "+str(self.Bz_max * ureg(self.field_units)))
        if 'phi0_t' in kwargs_list or PrintAll: print("Initial phase of time (phi0_t): "+str(self.phi0_t))
        if 'phi0_xy' in kwargs_list or PrintAll: print("Initial azimuthal angle (phi0_xy): "+str(self.phi0_xy))
        if 'p' in kwargs_list or PrintAll: print("Quantum number of radial order (p): "+str(self.p))
        if 'l' in kwargs_list or PrintAll: print("Quantum number of orbital angular momentum (l): "+str(self.l))
        if 's' in kwargs_list or PrintAll: print("Quantum number of spin angular momentum (s): "+str(self.s))
        # if 'polarization_modulation' in kwargs_list or PrintAll: print("Polarization modulation (T/F): "+str(self.polarization_modulation))
        if PrintAll: print(" ")
        if PrintAll: print("[7] Others")
        # if 'constant_bz' in kwargs_list or PrintAll: print("Applied constant Bz: "+str(self.constant_Bz * ureg(field_units)))
        if 'solve_time' in kwargs_list or PrintAll: print("Solve time (unit of tau_0): "+str(self.solve_time))
        if 'screenset' in kwargs_list or PrintAll: print("Screen Set (T/F): "+str(self.screenSet))

    def E_input_frame(self, ti,*,take_real: bool=True):
        Zv = np.zeros_like(self.Xv)
        Ex, Ey, Ez =  (E_LG_t(self.Xv, self.Yv, self.Zv, t=ti, w=self.w_EM, w0=self.w_0, E0=self.E0i, xc_Gauss=self.xc_Gauss, yc_Gauss=self.yc_Gauss,
                                 phi0_t=self.phi0_t, phi0_xy=self.phi0_xy, p=self.p, l=self.l, s=self.s,
                                 tau=self.tau, polarization_modulation=self.polarization_modulation,
                                 t_on=self.t_on, t_off=self.t_off, Bz=self.constant_Bz, time_evolute=self.time_evolute,
                                 angular_freq_units=self.angular_freq_units, length_units=self.length_units, E_field_units=self.E_field_units, time_units=self.time_units,))
        if take_real: return np.real(Ex)/self.E0i, np.real(Ey)/self.E0i
        else:         return Ex/self.E0i, Ey/self.E0i

    def set_state(self, **kwargs):
        markFreq = []
        for k, v in kwargs.items():
            if k in dir(self):
                setattr(self, k, v)
                if str(k)=='w_EM': markFreq = 'w_EM'
                if str(k)=='f_EM': markFreq = 'f_EM'

        # Update the calculated parameters:

        self.pearl_lambda = self.london_lambda**2/self.thickness # Updated value
        rn_resistivity = self.resistivity* ureg(self.resistivity_units)
        rho_conductivity = (1/rn_resistivity).to('1 / ohm / '+self.length_units)
        self.conductivity = rho_conductivity.to('1 / ('+self.resistivity_units+')').magnitude # Updated value
        self.disorder_epsilon = 1 # 1/self.temperature-1 # Updated value. Not allowed to be updated,
        xi_coherent = self.xi * ureg(self.length_units)
        lambdaL = self.london_lambda * ureg(self.length_units)
        d_thickness = self.thickness * ureg(self.length_units)
        lambdaL_eff = lambdaL**2/d_thickness
        mu_0 = cons.mu_0 * ureg('newton/A**2')
        tau_0 = mu_0 * rho_conductivity * lambdaL**2
        tau = tau_0.to(self.time_units).magnitude # Unit time, tau_GL in the real physical unit
        speed_of_light = cons.speed_of_light * ureg('m/s')
        revised_speed_of_light = (speed_of_light.to(f"{self.length_units}/{self.time_units}") * tau).magnitude # Speed of light after unit transformation
        self.c = revised_speed_of_light # Updated value
        self.tau = tau # Updated value, tau_GL
        self.tau_GL = tau # Updated value, tau_GL
        self.w_GL = 2*np.pi/tau # Updated value, the real physical value of angular freq of tau_GL
        self.f_GL = 1/tau # Updated value, the real physical value of freq of tau_GL
        self.Phi0 = cons.h/2/cons.e * ureg('J/A').to(self.field_units+'*'+self.length_units+'^2')
        self.Hc1 = self.Phi0/4/np.pi/mu_0/lambdaL**2*np.log(lambdaL/xi_coherent) # [ref] Gennes.P.D., pp.66, Eq (3-56), but no mu_0
        self.Hc2 = self.Phi0/2/np.pi/mu_0/xi_coherent**2 # [ref] Logan
        self.Hc  = self.Hc1/(np.pi/np.sqrt(24)*xi_coherent/lambdaL*np.log(lambdaL/xi_coherent)) # [ref] Gennes.P.D., pp.66, Eq (3-56), but no mu_0
        self.Bc1 = (mu_0*self.Hc1).to(self.field_units)
        self.Bc2 = (mu_0*self.Hc2).to(self.field_units)
        self.Bc  = (mu_0*self.Hc).to(self.field_units)
        self.A0 = (xi_coherent*self.Bc2)
        self.J0 = (4*xi_coherent*self.Bc2/mu_0/lambdaL**2).to(f"{self.current_units}/{self.length_units}^2")
        keyword_list = ['_l_x','_l_y','_c_r','_c_l']
        if self.light_source_type[:2].lower()=='lg' and self.light_source_type[-4:].lower() in keyword_list:
            p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = light_state_control_LaguerreGauss(keyword_of_state=self.light_source_type)
            self.p = p
            self.l = l
            self.s = s
            self.phi0_t = phi0_t
            self.phi0_xy = phi0_xy
            self.polarization_modulation = polarization_modulation
            self.output_file_head = output_file_head
        if markFreq == 'f_EM':
            self.w_EM = f_EM*(2*np.pi)
            self.f_EM = f_EM
        if markFreq == 'w_EM':
            self.w_EM = w_EM
            self.f_EM = w_EM/(2*np.pi)
        self.tau_EM = 1/self.f_EM # Unit of 1/tau_GL
        self.f_EM_unit = self.f_EM * self.f_GL
        self.w_EM_unit = self.w_EM * self.f_GL
        self.tau_EM_unit = self.tau_EM * self.tau_GL
    # [EM] Fundamental properties
        B0_EM = (1/self.c) * ureg(self.field_units) # Output value to TDGL detedmined by unit
        A0_EM = (1/self.w_EM) * ureg(f"{self.field_units} * {self.length_units}") # Output value to TDGL detedmined by unit
        E0_EM = B0_EM * speed_of_light # Output value to TDGL detedmined by unit
        self.E0_GL = self.A0 * self.f_GL * ureg('1/'+self.time_units) # Output value to TDGL detedmined by unit
        self.E0 = (self.E_amp * self.E0_GL).to(self.E_field_units) # Input amplitude of electric field
        self.E0i = self.E_amp # Updated value
        # * ureg("1/angular_freq_units").to("time_units").magnitude
        self.lambda_EM = (speed_of_light * self.tau_EM_unit * ureg(self.time_units).to("s")).to(self.length_units).magnitude
        self.B_EM_inplane = self.E0i*B0_EM
        self.A_EM_inplane = self.E0i*A0_EM

def printTAB(self,*,figsize=(5,2.5),dpi=100,fontsize=10,tab_scale=(1.5,2),round_num=2) -> plt.Figure:
        eps = 1-self.temperature
        fig, axs = plt.subplots(2,1,figsize=figsize,dpi=dpi)
    
        collabel=("$\\xi_0$\n ("+self.length_units+")", 
                  "$\\lambda_{L,0}$ or $\\Gamma$\n ("+self.length_units+")", 
                  "$r_n=1/\\rho$\n ("+self.resistivity_units+")",
                  "$\\tau_{GL,0}$\n ("+self.time_units+")",
                  "$\\omega_{GL,0}/2\\pi$\n ("+self.angular_freq_units+")",
                  "$\\xi_0/\\tau_{GL,0}$\n ("+self.length_units+'/'+self.time_units+")",)
        clust_data = ([round(self.xi*np.sqrt(eps),round_num),
                      round(max([self.london_lambda,self.pearl_lambda])*np.sqrt(eps),round_num),
                      round(self.resistivity,round_num),
                      round(self.tau*eps,round_num),
                      round(self.f_GL/eps,round_num),
                      round((self.xi*np.sqrt(eps))/(self.tau*eps),round_num)],)
        axs[0].axis('off')
        axs[0].set_title('Temperature: '+str(0)+' $T_c$',fontsize=fontsize,loc='left')
        # axs[0].set_fontsize(fontsize)
        TAB = axs[0].table(cellText=clust_data,colLabels=collabel,loc='center',fontsize=fontsize)
        # TAB.axis('off')
        TAB.scale(tab_scale[0], tab_scale[1])

        collabel=("$\\xi(T)$\n ("+self.length_units+")", 
                  "$\\lambda_{L}(T)$ or $\\Gamma(T)$\n ("+self.length_units+")", 
                  "$r_n=1/\\rho$\n ("+self.resistivity_units+")",
                  "$\\tau_{GL}(T)$\n ("+self.time_units+")",
                  "$\\omega_{GL}(T)/2\\pi$\n ("+self.angular_freq_units+")",
                  "$\\xi(T)/\\tau_{GL}(T)$\n ("+self.length_units+'/'+self.time_units+")",)
        clust_data = ([round(self.xi,round_num),
                      round(max([self.london_lambda,self.pearl_lambda]),round_num),
                      round(self.resistivity,round_num),
                      round(self.tau,round_num),
                      round(self.f_GL,round_num),
                      round(self.xi/self.tau,round_num)],)
        axs[1].axis('off')
        TAB = axs[1].set_title('Temperature: '+str(self.temperature)+' $T_c$, $\\epsilon=1-T/T_c$='+str(round(1-self.temperature,round_num)),fontsize=fontsize,loc='left')
        TAB = axs[1].table(cellText=clust_data,colLabels=collabel,loc='center',fontsize=fontsize)
        TAB.scale(tab_scale[0], tab_scale[1])

        plt.show()
        return fig, axs

  

        
def dump(obj):
  for attr in dir(obj):
    print("obj.%s = %r" % (attr, getattr(obj, attr)))


def plot_polarization(X,Y,E_x,E_y,*,E0i:float=1.0,title:str='',figsize:(3, 3),scale:float=12,dpi:float=100,field_units:str='mT'):
    fig = plt.figure(figsize=figsize,constrained_layout=True,dpi=dpi)
    plt.title(title)
    plt.quiver(X,Y,E_x/E0i,E_y/E0i,scale=12, scale_units='x',width=0.1*abs(X[2]-X[1]))
    plt.xlabel('x ($\mu$m)')
    plt.ylabel('y ($\mu$m)')

def plot_EM(X,Y,E_x,E_y,B_z,*,E0i:float=1.0,title:str='',figsize:(6, 3),scale:float=12,dpi:float=100,width_quiver:float=0.01):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize,dpi=dpi) #constrained_layout=True,
    fig.suptitle(title)
    ax1.quiver(X,Y,E_x/E0i,E_y/E0i,scale=scale, scale_units='x',width=width_quiver*abs(X[2]-X[1]))
    ax1.set_xlabel('x ($\mu$m)')
    ax1.set_ylabel('y ($\mu$m)')
    # ax1.text(min(X)*.95, max(Y)*.85, '$|E_{0}|$: '+str(E0i), horizontalalignment='left', fontsize='large')
    ax1.set_aspect("equal")
    Xv, Yv = np.meshgrid(X, Y)
    Bzmax, Bzmin = [B_z.max(), B_z.min()]
    contour_Bz = ax2.contourf(Xv, Yv, B_z, levels=50, linewidths=0.0, cmap="PRGn",vmin=Bzmin,vmax=Bzmax)
    cbar = plt.colorbar(contour_Bz)
    cbar.set_label('Normalized $\Phi_{B}$')
    ax2.set_xlabel('x ($\mu$m)')
    ax2.set_ylabel('y ($\mu$m)')
    # ax2.text(min(X)*.95, max(Y)*.85, '$|B_{z,0}|$: '+str(E0i), horizontalalignment='left', fontsize='large')
    ax2.set_aspect("equal")

def v_grid_generation(p1,p2,p3,p4,quiver_mesh_n): # (-width/2,width/2,-height/2,height/2,quiver_mesh_n)
    X = np.linspace(p1,p2,quiver_mesh_n)
    Y = np.linspace(p3,p4,quiver_mesh_n)
    Xv, Yv = np.meshgrid(X, Y)
    Zv = np.zeros_like(Xv)
    return X, Y, Xv, Yv, Zv
