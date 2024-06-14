
### LG beam

from scipy.special import genlaguerre
import scipy.constants as cons
from tdgl.parameter import Parameter
import pint # https://pint.readthedocs.io/en/0.10.1/tutorial.html

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
    # assert isinstance(Bz, (float, str, pint.Quantity)), type(Bz)
    positions = np.atleast_2d(positions)
    # assert positions.shape[1] == 3, positions.shape
    # if not isinstance(positions, pint.Quantity):
    #     positions = positions * ureg("meter")
    # if isinstance(Bz, str):
    #     Bz = ureg(Bz)
    # if isinstance(Bz, float):
    #     Bz = Bz * ureg("tesla")
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
    phi_t = np.copy(phi0_t) # !!!!!! ONLY CONSTANT PHASE REMAINS !!!!!!
    if time_evolute: phi_t = phi_t + w*ti
    phiGouy = (2*p+np.abs(l)+1)*np.arctan(Z/zR) # Gouy phase
    u = E0 * (np.sqrt(2)*r/wz)**l * genlaguerre(p,l)(2*r**2/wz**2) * w0/wz * np.exp(-r**2/wz**2 -1j*phiGouy +1j*(l*phi +k*z +k*r**2/2/(z-1j*zR)))

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
        p: float = 0.0,
        l: float = 0.0,
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
        E0=E0i,
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

### A of LG beam on GPU (cupy)

# import cupy as cp

# def uniform_Bz_vector_potential_cupy(
#     positions: np.ndarray,
#     Bz: Union[float, str, pint.Quantity],
# ) -> cp.ndarray:
#     """Calculates the magnetic vector potential [Ax, Ay, Az] at ``positions``
#     due uniform magnetic field along the z-axis with strength ``Bz``.

#     Args:
#         positions: Shape (n, 3) array of (x, y, z) positions in meters at which to
#             evaluate the vector potential.
#         Bz: The strength of the uniform field, as a pint-parseable string,
#             a pint.Quantity, or a float with units of Tesla.

#     Returns:
#         Shape (n, 3) array of the vector potential [Ax, Ay, Az] at ``positions``
#         in units of Tesla * meter.
#     """
#     # assert isinstance(Bz, (float, str, pint.Quantity)), type(Bz)
#     positions = cp.atleast_2d(positions)
#     # assert positions.shape[1] == 3, positions.shape
#     # if not isinstance(positions, pint.Quantity):
#     #     positions = positions * ureg("meter")
#     # if isinstance(Bz, str):
#     #     Bz = ureg(Bz)
#     # if isinstance(Bz, float):
#     #     Bz = Bz * ureg("tesla")
#     xs = cp.array(positions[:, 0])
#     ys = cp.array(positions[:, 1])
#     dx = cp.ptp(xs)
#     dy = cp.ptp(ys)
#     xs = xs - (xs.min() + dx / 2)
#     ys = ys - (ys.min() + dy / 2)
#     Ax = -Bz * ys / 2
#     Ay = Bz * xs / 2
#     A = cp.stack([Ax, Ay, np.zeros_like(Ax)], axis=1)
#     return A

# def constant_field_vector_potential_cupy(
#     x,
#     y,
#     z,
#     *,
#     Bz: float,
#     field_units: str = "mT",
#     length_units: str = "um",
# ):
#     if z.ndim == 0:
#         z = z * cp.ones_like(x)
#     positions = cp.array([x.squeeze(), y.squeeze(), z.squeeze()]).T
#     # Bz = Bz * ureg(field_units)
#     A = uniform_Bz_vector_potential_cupy(positions, Bz)
#     return A

def A_LG_t_cupy(x, y, z, *, t,
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
                   take_complex=False,
                   time_evolute: bool = True,

):
    """ Vector potential of Laguerre-Gaussian beam of p-th degree of mode and l-th order of mode
    From E to A, A = -iE/w
    """
    X = cp.array(x-xc_Gauss) * 1 # ureg(length_units).to("m").magnitude   # mesh of x-axis [um]
    Y = cp.array(y-yc_Gauss) * 1 # ureg(length_units).to("m").magnitude   # mesh of y-axis [um]
    Z = cp.array(z-z0) * 1 # ureg(length_units)).to("m").magnitude   # mesh of y-axis [um]
    w0 = w0 * 1 # ureg(length_units) .to("m").magnitude   # Beam waist (suppose w0_x =w0_y) [um]
    k = w/c   # Wavenumber
    ti = (t)#(t * ureg(time_units)).to("s").magnitude   # Input time
    # E0 = (E0 * ureg(E_field_units)).to("newton / coulomb").magnitude   # Electrical field |E0|
    # print(ti)
    zR = cp.array(n*k*w0**2/2)   # Rayleigh length, n: refractive index
    # print(zR)
    # print(w0)
    wz = cp.array(w0 * cp.sqrt(1+(Z/zR)**2))   # Spot radius
    zeta = cp.array(Z/zR)

    r = cp.sqrt(X**2+Y**2)   # Radius
    phi = cp.angle(X+1j*Y)   # Azimuthal angle
    if polarization.lower()=='x' or polarization.lower()=='linear x': s, phi0_xy = [0.0,0.0]
    if polarization.lower()=='y' or polarization.lower()=='linear y': s, phi0_xy = [0.0,np.pi/2]
    if polarization.lower()=='lc' or polarization.lower()=='left circular': s, phi0_xy = [1,np.pi/4]
    if polarization.lower()=='rc' or polarization.lower()=='right circular': s, phi0_xy = [-1,np.pi/4]
    phi_t = (phi0_t) # !!!!!! ONLY CONSTANT PHASE REMAINS !!!!!!
    if time_evolute: phi_t = phi_t + w*ti
    phiGouy = cp.array((2*p+cp.abs(l)+1)*cp.arctan(Z/zR)) # Gouy phase
    # genlaguerre: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.genlaguerre.html
    # u = E0 * (cp.sqrt(2)*r/(wz))**l * cp.array(genlaguerre(p,l)(2*np.array(r**2)/np.array(wz**2))) * w0/wz * cp.exp(-r**2/wz**2 -1j*phiGouy +1j*(l*phi +k*Z +k*r**2/2/(Z-1j*zR)))
    u = E0 * (cp.sqrt(2)*r/(wz))**l * 1 * w0/wz * cp.exp(-r**2/wz**2 -1j*phiGouy +1j*(l*phi +k*Z +k*r**2/2/(Z-1j*zR)))
    if t>t_off or t<t_on: t_step = 0
    else: t_step = 1

    pol_m_x, pol_m_y = [cp.cos(phi0_xy),cp.sin(phi0_xy)]
    if polarization_modulation: pol_m_x, pol_m_y = (cp.abs([cp.cos(l*phi+phi0_xy), cp.sin(l*phi+phi0_xy)]))
    Ex = u * cp.exp(-1j*(phi_t)) * pol_m_x * t_step
    Ey = u * cp.exp(-1j*(phi_t+s*cp.pi/2)) * pol_m_y * t_step

    Ax = -1j/w*Ex
    Ay = -1j/w*Ey
    Az = cp.zeros_like(Ax)
    Ax[np.isnan(Ax)] = 0
    Ax[np.isinf(Ax)] = 0
    Ay[np.isnan(Ay)] = 0
    Ay[np.isinf(Ay)] = 0
    A_constBz = 0#constant_field_vector_potential_cupy(x, y, z, Bz=Bz, field_units=field_units, length_units=length_units)
    A = cp.stack([cp.real(Ax), cp.real(Ay), cp.real(Az)], axis=1) + A_constBz
    if take_complex: A = cp.stack([(Ax), (Ay), (Az)], axis=1) + A_constBz
    # A = np.array([np.real(Ax), np.real(Ay), np.real(Az)]).T
    return cp.asnumpy(A)#.to(f"{field_units} * {length_units}").magnitude

def A_LG_cupy(*,
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
        p: float = 0.0,
        l: float = 0.0,
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
        A_LG_t_cupy,
        w=w,
        w0=w0,
        E0=E0i,
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

### Bz component

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
    phi_t = np.copy(phi0_t) # !!!!!! ONLY CONSTANT PHASE REMAINS !!!!!!
    if time_evolute: phi_t = phi_t + w*ti
    phiGouy = (2*p+np.abs(l)+1)*np.arctan(Z/zR) # Gouy phase
    u = E0 * (np.sqrt(2)*r/wz)**l * genlaguerre(p,l)(2*r**2/wz**2) * w0/wz * np.exp(-r**2/wz**2 -1j*phiGouy +1j*(l*phi +k*z +k*r**2/2/(z-1j*zR)))

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

def find_max_Bz(Xv,Yv,E0i,constant_Bz,c,w):
    Ex, Ey = E_input_frame(Xv,Yv,0,take_real=False)
    Bx, By, Bz1 = E2Bv(Xv,Yv,E0i*Ex,E0i*Ey,0,c,w)
    Ex, Ey = E_input_frame(Xv,Yv,2*np.pi/4/w,take_real=False)
    Bx, By, Bz2 = E2Bv(Xv,Yv,E0i*Ex,E0i*Ey,0,c,w)
    Ex, Ey = E_input_frame(Xv,Yv,2*np.pi/2/w,take_real=False)
    Bx, By, Bz3 = E2Bv(Xv,Yv,E0i*Ex,E0i*Ey,0,c,w)
    return max([abs(np.real(Bz1)).max(), abs(np.real(Bz2)).max(), abs(np.real(Bz3)).max()]) + constant_Bz


### ------------------------------------------------------------------------------------------ ###

### Other useful functions

def light_state_contral(keyword_of_state):

    ''' Select the parameters for optical states
    options: 'lg00_l_x','lg00_l_y','lg00_c_l','lg00_c_r','lg01_l_x','lg01_l_y','lg01_c_c','lg01_c_r','lg01_c_a'
    return p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head
    '''

## Linear polatization, Gauss:  for linear polarization, phi0_xy could be any number and s=0
    if keyword_of_state.lower()=='lg00_l_x': p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = [0,0,0,0,0, False,'LG00_linear_x']
    if keyword_of_state.lower()=='lg00_l_y': p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = [0,0,0,0,np.pi/2, False,'LG00_linear_y']
 ## Circular polatization, Gauss: for Circular polarization, phi0_xy=pi/4 and s=+-1
    if keyword_of_state.lower()=='lg00_c_l': p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = [0,0,1,0,np.pi/4, False,'LG00_circular_l']
    if keyword_of_state.lower()=='lg00_c_r': p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = [0,0,-1,0,np.pi/4, False,'LG00_circular_r']
## Linear polatization, LG01: for linear polarization, phi0_xy could be any number and s=0
    if keyword_of_state.lower()=='lg01_l_x': p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = [0,1,0,0,0, False,'LG01_linear_x']
    if keyword_of_state.lower()=='lg01_l_y': p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = [0,1,0,0,np.pi/2, False,'LG01_linear_y']
## Circular polatization, LG01 (Radial + Azimuthal): for Circular polarization, phi0_xy=pi/4 and s=+-1
    if keyword_of_state.lower()=='lg01_c_c': p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = [0,1,1,0,np.pi/4, False,'LG01_circular_c']
    if keyword_of_state.lower()=='lg01_c_l': p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = [0,1,1,0,np.pi/4, False,'LG01_circular_l']
    if keyword_of_state.lower()=='lg01_c_r': p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = [0,1,-1,0,np.pi/4, False,'LG01_circular_r']
## Radial polatization, LG01: for Radial polarization, s=+1, phi0_xy=0, polarization_modulation = True
    if keyword_of_state.lower()=='lg01_c_r': p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = [0,1,1,0,0, True,'LG01_circular_r']
## Azimuthal polatization, LG01: for Azimuthal polarization, s=+1, phi0_xy=np.pi/2, polarization_modulation = True
    if keyword_of_state.lower()=='lg01_c_a': p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = [0,1,1,0,np.pi/2, True,'LG01_circular_a']

    return p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head

def plot_polarization(X,Y,E_x,E_y,*,E0i:float=1.0,title:str='',figsize:(3, 3),scale:float=12,dpi:float=100):
    fig = plt.figure(figsize=figsize,constrained_layout=True,dpi=dpi)
    plt.title(title)
    plt.quiver(X,Y,E_x/E0i,E_y/E0i,scale=12, scale_units='x',width=0.1*abs(X[2]-X[1]))
    plt.xlabel('x ($\mu$m)')
    plt.ylabel('y ($\mu$m)')

def plot_EM(X,Y,E_x,E_y,B_z,*,E0i:float=1.0,title:str='',figsize:(6, 3),scale:float=12,dpi:float=100,take_Bz_range:bool=False,width_quiver:float=0.01):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize,dpi=dpi) #constrained_layout=True,
    fig.suptitle(title)
    ax1.quiver(X,Y,E_x/E0i,E_y/E0i,scale=quiver_scale, scale_units='x',width=width_quiver*abs(X[2]-X[1]))
    ax1.set_xlabel('x ($\mu$m)')
    ax1.set_ylabel('y ($\mu$m)')
    # ax1.text(min(X)*.95, max(Y)*.85, '$|E_{0}|$: '+str(E0i), horizontalalignment='left', fontsize='large')
    ax1.set_aspect("equal")
    if take_Bz_range: Bzmax, Bzmin = [find_max_Bz(Xv,Yv,E0i,constant_Bz,c,w_input), -find_max_Bz(Xv,Yv,E0i,constant_Bz,c,w_input)]
    else: Bzmax, Bzmin = [B_z.max(), B_z.min()]
    contour_Bz = ax2.contourf(Xv, Yv, B_z, levels=50, linewidths=0.0, cmap="PRGn",vmin=Bzmin,vmax=Bzmax)
    cbar = plt.colorbar(contour_Bz)
    cbar.set_label('$B_{z}$ ['+field_units+']')
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

### Unit transfer and check

from datetime import datetime
import pytz

def Unit_check_save(E0i,w0,w,light_source_state,length_units,E_field_units,angular_freq_units,subtitle,solve_time,screenSet,folder_name,Bz_max):

    ## Unit transfer
    w0 = (w0 * ureg(length_units)).to(length_units).magnitude   # Beam waist (suppose w0_x =w0_y)
    tau = tau_0.to(time_units).magnitude   # Unit time
    c = (cons.speed_of_light * ureg('m/s').to(length_units+'/'+time_units) * tau).magnitude   # speed of light (3e8 * tau  m/s)
    p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head = light_state_contral(light_source_state)
    B0 = (1/c) * ureg(field_units)
    A0 = (1/w) * ureg(f"{field_units} * {length_units}")
    E0 = B0 * speed_of_light

    ## Save to data
    file_name = folder_name+'/'+subtitle+output_file_head+"_rec.txt"
    f = open(file_name, "w")
    f.write("File name "+file_name+"\n")
    f.write("Time now (at Stockholm):"+str(datetime.now(pytz.timezone('Europe/Stockholm')))+"\n")
    f.write(":: Parameter record ::\n")
    f.write("[1] Length scale of sample\n")
    f.write("Coherent length (xi): "+str(xi_coherent)+"\n")
    f.write("London penetration depth (london_lambda): "+str(lambdaL)+"\n")
    f.write("Thickness (thickness): "+str(d_thickness)+"\n")
    f.write("Screen length (lambdaL**2/thickness): "+str(screen_length)+"\n")
    f.write("Ratio of length kapa (lambdaL/xi): "+str(lambdaL/xi_coherent)+"\n")
    f.write("[2] Condictivity and  of sample\n")
    f.write("Resistivity (resistivity): "+str(resistivity)+"\n")
    f.write("Condictivity (1/resistivity): "+str(condictivity)+"\n")
    f.write("[3] Time scale of sample\n")
    f.write("Characteristic timescale (tau_0): "+str(tau_0.to(time_units))+"\n") # characteristic timescale for this TDGL model
    f.write("Characteristic rate (1/tau_0): "+str((1/tau_0).to('THz'))+"\n")
    f.write("Speed of light (unit of tau_0): "+str(revised_speed_of_light)+"\n")
    f.write("c input into A_LG (c): "+str(c)+"\n")
    f.write("tau input into A_LG (tau_0.to(s)): "+str(tau)+"\n")
    f.write("k calculated inside A_LG (w/c): "+str(w/c)+"\n")
    f.write("[4] Gap of sample\n")
    f.write("Given inelastic coupling rate (rate_eph): "+str(E_eph)+"\n")
    f.write("Given SC gap energy (gap_0): "+str(E_gap_0)+"\n")
    f.write("Frequency of SC gap (gap_0*0.242): "+str(F_gap_0)+"\n")
    f.write("Strength of inelastic scattering (gamma=2/rate_eph*gap_0): "+str(gamma)+"\n")
    f.write("[5] Parameters of light source\n")
    f.write("Beam size (2w0): "+str(2*w0 * ureg(length_units))+"\n")
    wavelength_of_light = (cons.speed_of_light * ureg('m/s') / ((w/2/np.pi/tau_0).to('THz'))).to(length_units)
    f.write("Wavelength: "+str(wavelength_of_light)+"\n")
    f.write("Beam size / wavelength: "+str((2*w0 * ureg(length_units)) / wavelength_of_light)+"\n")
    if (2*w0 * ureg(length_units) / wavelength_of_light).magnitude<1: print('[Warning] !!! ---Beam size is smaller than wavelength.--- !!!')
    if time_evolute: f.write("Angular frequency of light (w, unit of tau): "+str(w)+"\n")
    if time_evolute: f.write("Angular frequency of light (w, real value): "+str((w/tau_0).to(angular_freq_units))+"\n")
    if time_evolute: f.write("Frequency of light (w/2pi, unit of tau): "+str(w/2/np.pi)+"\n")
    if time_evolute: f.write("Frequency of light (w/2pi, real value): "+str((w/2/np.pi/tau_0).to(angular_freq_units))+"\n")
    if time_evolute!=True: f.write("Time evolution setting is False (no frequency is used).")
    f.write("EM field switch-on at (unit of tau): "+str(t_on)+"\n")
    f.write("EM field switch-off at (unit of tau): "+str(t_off)+"\n")
    f.write("|E0| of light (case of E0=1): "+str(E0.to(E_field_units))+"\n")
    f.write("In-plane |B0| of light (|E0|/c): "+str(B0)+"\n")
    f.write("In-plane |A0| of light (|E0|/w): "+str(A0)+"\n")
    f.write("Check ratio of 2pi|A0|/|B0|: "+str(2*np.pi*A0/B0)+"\n")
    f.write("|Bz| of light: "+str(Bz_max * ureg(field_units))+"\n")
    f.write("Range of input |E0| (E0i): "+str(E0i)+"\n")
    f.write("Degree of LG mode (p): "+str(p)+"\n")
    f.write("Order of LG mode (l): "+str(l)+"\n")
    f.write("Spin number (s): "+str(s)+"\n")
    f.write("Initial phase of time (phi0_t): "+str(phi0_t)+"\n")
    f.write("Initial azimuthal angle (phi0_xy): "+str(phi0_xy)+"\n")
    f.write("Polarization modulation (T/F): "+str(polarization_modulation)+"\n")
    f.write("output_file_head (string): "+output_file_head+"\n")
    f.write("[6] Critical E and B of SC\n")
    f.write("Quantum magnetic flux (Phi0 = h/2e): "+str(cons.h/2/cons.e * ureg('J/A').to('W/A*s'))+"\n")
    f.write("Hc1 of SC (Phi0/4pi/lambda^2*ln[lambda/xi]): "+str(Hc1)+"\n") # [ref] Gennes.P.D., pp.66, Eq (3-56)
    f.write("Hc2 of SC (Phi0/2pi xi^2): "+str(Hc2)+"\n")
    f.write("Hc of SC: "+str(Hc)+"\n")
    f.write("Bc1 of SC (mu0Hc1): "+str((mu_0*Hc1).to('T'))+"\n") # lower critical field
    f.write("Bc2 of SC (mu0Hc2): "+str((mu_0*Hc2).to('T'))+"\n") # upper critical field
    f.write("Bc of SC: "+str(Bc.to('T'))+"\n")
    f.write("Bc*2pi*w0^2/Phi0: "+str(mu_0*Hc*2*np.pi*(w0*ureg(length_units))**2/Phi0)+"\n")
    f.write("Bc2*2pi*w0^2/Phi0: "+str(mu_0*Hc2*2*np.pi*(w0*ureg(length_units))**2/Phi0)+"\n")
    f.write("Critical vector potential A0 (xi*Bc2): "+str(xi_coherent*Bc2)+"\n")
    f.write("Unit current density (J0 = 4*xi*Bc2/mu_0/lambdaL**2): "+str(4*xi_coherent*Bc2/mu_0/lambdaL**2)+"\n")
    f.write("|B0|/Bc: "+str(B0/((mu_0*Hc).to(field_units)))+"\n") # upper critical field
    f.write("|Bz|/Bc: "+str(Bz_max*ureg(field_units)/((mu_0*Hc).to(field_units)))+"\n") # upper critical field
    f.write("Fermi velocity: "+str(vF.to('m/s'))+"\n")
    f.write("Fermi velocity (unit of c): "+str(vF.to('m/s')/(speed_of_light))+"\n")
    f.write("Condensation energy (Bc**2/2/mu0): "+str(Econd.to('eV/'+length_units+'**3'))+"\n")
    f.write("[7] Others\n")
    f.write("Applied constant Bz: "+str(constant_Bz * ureg(field_units))+"\n")
    f.write("Solve time (unit of tau_0): "+str(solve_time)+"\n")
    f.write("Screen Set (T/F): "+str(screenSet)+"\n")
    f.close()

    ## Print and check
    f = open(file_name, "r")
    print(f.read())
    print('... Save to file: '+file_name)
    f.close()

    return tau, c


""" Example of EMwave_G_cir:
applied_vector_potential = EMwave_G_cir(w=1e8, E0=1, phi0=0)
"""

""" Example of A_LG:
E0i = 10
w_input = 2*np.pi
w0 = 0.1
# p, l, s, phi0_t, phi0_xy, polarization_modulation = [0,0,0,0,0, False]       ## Linear polatization, Gauss:  for linear polarization, phi0_xy could be any number and s=0
# p, l, s, phi0_t, phi0_xy, polarization_modulation = [0,0,1,0,np.pi/4, False] ## Circular polatization, Gauss: for Circular polarization, phi0_xy=pi/4 and s=+-1
# p, l, s, phi0_t, phi0_xy, polarization_modulation = [0,1,0,0,0, False]       ## Linear polatization, LG01: for linear polarization, phi0_xy could be any number and s=0
# p, l, s, phi0_t, phi0_xy, polarization_modulation = [0,1,1,0,np.pi/4, False] ## Circular polatization, LG01 (Radial + Azimuthal): for Circular polarization, phi0_xy=pi/4 and s=+-1
# p, l, s, phi0_t, phi0_xy, polarization_modulation = [0,1,1,0,0, True]        ## Radial polatization, LG01: for Radial polarization, s=+1, phi0_xy=0, polarization_modulation = True
# p, l, s, phi0_t, phi0_xy, polarization_modulation = [0,1,1,0,np.pi/2, True]  ## Azimuthal platization, LG01: for Azimuthal polarization, s=+1, phi0_xy=np.pi/2, polarization_modulation = True
applied_vector_potential = A_LG(w=w_input, w0=w0, E0=E0i,
                             phi0_t=phi0_t, phi0_xy=phi0_xy, p=p, l=l, s=s, c=c,
                             tau=tau_0.to(time_units).magnitude, polarization_modulation=polarization_modulation,
                             angular_freq_units=angular_freq_units, length_units=length_units, E_field_units=E_field_units, time_units=time_units,)
"""

""" Example of A_LG_step:
E0i = 10
w_input = 2*np.pi
w0 = 0.1
Bz = 1
t_on, t_off = [0.0, 2.0]
applied_vector_potential = A_LG(w=w_input, w0=w0, E0=E0i,
                             phi0_t=phi0_t, phi0_xy=phi0_xy, p=p, l=l, s=s, c=c,
                             t_on=t_on, t_off=t_off, Bz=Bz, time_evolute=True,
                             tau=tau_0.to(time_units).magnitude, polarization_modulation=polarization_modulation, field_units=field_units,
                             angular_freq_units=angular_freq_units, length_units=length_units, E_field_units=E_field_units, time_units=time_units,)
"""

""" Example of A2B(x, y, z, A)

applied_vector_potential_t = A_LG_t(x=x, y=y, z=z, t=t, w=w_input, w0=w0, E0=E0i,
                             phi0_t=phi0_t, phi0_xy=phi0_xy, p=p, l=l, s=s, c=c,
                             t_on=t_on, t_off=t_off, Bz=Bz, time_evolute=True,
                             tau=tau_0.to(time_units).magnitude, polarization_modulation=polarization_modulation, field_units=field_units,
                             angular_freq_units=angular_freq_units, length_units=length_units, E_field_units=E_field_units, time_units=time_units,)
B = A2B(x, y, z, applied_vector_potential_t)
# Bx = B.x
# By = B.y
# Bz = B.z

E_LG_beam_t = E_LG_t(x=xv, y=yv, z=zv, t=t, w=w_input, w0=w0, E0=E0i,
                             phi0_t=phi0_t, phi0_xy=phi0_xy, p=p, l=l, s=s, c=c,
                             t_on=t_on, t_off=t_off, Bz=Bz_constant, time_evolute=True,
                             tau=tau_0.to(time_units).magnitude, polarization_modulation=polarization_modulation, field_units=field_units,
                             angular_freq_units=angular_freq_units, length_units=length_units, E_field_units=E_field_units, time_units=time_units,)
applied_magnetic_field_t = E2Bv(xv,yv,E,Bz_constant,c,w)
# Bx = applied_magnetic_field_t[:,0]
# By = applied_magnetic_field_t[:,1]
# Bz = applied_magnetic_field_t[:,2]
"""

""" Test of avalability of function: (for the error of saving data)
# pickle.dumps(applied_vector_potential)
# np.dtype(applied_vector_potential)
# print(applied_vector_potential.dtype)
import pickle
import cloudpickle

applied_vector_potential2 = A_LG(w=1, w0=2, E0=3, phi0_t=4, phi0_xy=5, p=1, l=2, s=3, tau=1)

# np.void(cloudpickle.dumps(applied_vector_potential2))
applied_vector_potential2 = A_LG(w=w_input, w0=w0, E0=E0i,
                             phi0_t=phi0_t, phi0_xy=phi0_xy, p=p, l=l, s=s, c=c,
                             tau=tau_0.to(time_units).magnitude, polarization_modulation=polarization_modulation,
                             angular_freq_units=angular_freq_units, length_units=length_units, E_field_units=E_field_units, time_units=time_units,)
# pickle.dumps(applied_vector_potential2)
cloudpickle.dumps(applied_vector_potential2)
"""

""" Select the parameters for optical states via "light_state_contral(keyword_of_state)":
options: 'lg00_l_x','lg00_l_y','lg00_c_l','lg00_c_r','lg01_l_x','lg01_l_y','lg01_c_c','lg01_c_r','lg01_c_a'
return p, l, s, phi0_t, phi0_xy, polarization_modulation, output_file_head
"""

"""Example of Unit_check_save:
Unit_check_save(E0i,w0,w,light_source_type,length_units,E_field_units,angular_freq_units,subtitle)
"""

""" Example of MAKE_ANIMATIONS:
if MAKE_ANIMATIONS:
    test_video = make_video_from_solution(
        test_solution,
        quantities=["order_parameter", "phase", "scalar_potential","vorticity"],
        figsize=(22, 4),
    )
    display(test_video)
"""

# pip install cupy-cuda11x
# pip install cupy-cuda12x

""" Done."""
