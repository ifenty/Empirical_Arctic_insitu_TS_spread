import numpy as np


def calculate_adiabatic_T_gradient(S,T,P):
    """
    % DESCRIPTION:
    %    Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
    """
    
    a0 =  3.5803e-5
    a1 = +8.5258e-6
    a2 = -6.836e-8
    a3 =  6.6228e-10

    b0 = +1.8932e-6
    b1 = -4.2393e-8

    c0 = +1.8741e-8
    c1 = -6.7795e-10
    c2 = +8.733e-12
    c3 = -5.4481e-14

    d0 = -1.1351e-10
    d1 =  2.7759e-12

    e0 = -4.6206e-13
    e1 = +1.8676e-14
    e2 = -2.1687e-16

    adiabatic_T_gradient = (
            a0 
            + (a1 + (a2 + a3 * T) * T) * T 
            + (b0 + b1 * T) * (S-35) 
            + ((c0 + (c1 + (c2 + c3 * T) * T) * T) + (d0 + d1 * T) * (S-35)) * P 
            + (e0 + (e1 + e2 * T) * T ) * P * P
            )

    return adiabatic_T_gradient


def depth_to_pressure(depth, lat):
    """Pressure (dbar) from depth (m) and latitude (deg) as per UNESCO 1983 (Saunders and Fofonoff); broadcasts."""
    magic_number = 5.92e-3 + np.sin(np.abs(lat) * np.pi / 180) ** 2 * 5.25e-3
    return ((1 - magic_number) - np.sqrt(((1 - magic_number) ** 2) - (8.84e-6 * depth))) / 4.42e-6


def calculate_potential_T(T, S, depth, lat):
    """
    % DESCRIPTION:
    %    Potential temperature (referenced to the surface) from in situ T (degC), practical salinity S and
    %    depth (m) / latitude (deg), as per UNESCO 1983 (4th-order Runge-Kutta on the adiabatic gradient).

    Works on any broadcast-compatible numpy / xarray inputs, e.g. T, S (depth, n), depth (depth, 1), lat (n,).
    Missing T or S gives NaN: the caller decides how to fill missing S (never persist a filled S).
    """
    profile_pressures = depth_to_pressure(depth, lat)
    del_profile_pressures = -profile_pressures  # reference pressure is 0

    # theta1
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(S, T, profile_pressures)
    theta = T + 0.5 * del_theta
    q_factor = del_theta

    # theta2
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(S, theta, profile_pressures + 0.5 * del_profile_pressures)
    theta = theta + (1 - 1/np.sqrt(2)) * (del_theta - q_factor)
    q_factor = (2 - np.sqrt(2)) * del_theta + (-2 + 3/np.sqrt(2)) * q_factor

    # theta3
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(S, theta, profile_pressures + 0.5 * del_profile_pressures)
    theta = theta + (1 + 1/np.sqrt(2)) * (del_theta - q_factor)
    q_factor = (2 + np.sqrt(2)) * del_theta + (-2 - 3/np.sqrt(2)) * q_factor

    # theta4
    del_theta = del_profile_pressures * calculate_adiabatic_T_gradient(S, theta, profile_pressures + del_profile_pressures)
    return theta + (del_theta - 2 * q_factor)/6
