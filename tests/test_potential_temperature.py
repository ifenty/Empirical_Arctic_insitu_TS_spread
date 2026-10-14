import numpy as np
import pytest

from aotd.potential_temperature import calculate_adiabatic_T_gradient, calculate_potential_T, depth_to_pressure


def test_adiabatic_gradient_unesco_check_value():
    # UNESCO (1983) technical paper 44 check value: S=40, T=40 degC, P=10000 dbar -> 3.255976e-4 degC/dbar
    assert calculate_adiabatic_T_gradient(40.0, 40.0, 10000.0) == pytest.approx(3.255976e-4, rel=1e-6)


def test_surface_is_identity():
    assert calculate_potential_T(np.array([-1.2, 3.0]), np.array([33.0, 35.0]), 0.0, 75.0) == pytest.approx([-1.2, 3.0])


def test_matches_gsw():
    gsw = pytest.importorskip("gsw")
    depth = np.array([100.0, 1000.0, 3000.0])
    T, S = np.array([0.5, 1.0, 0.0]), np.array([34.5, 34.9, 34.95])
    P = depth_to_pressure(depth, 80.0)
    assert calculate_potential_T(T, S, depth, 80.0) == pytest.approx(gsw.pt_from_t(S, T, P, 0.0), abs=2e-3)


def test_nan_propagates():
    out = calculate_potential_T(np.array([np.nan, 1.0]), np.array([34.0, np.nan]), 500.0, 80.0)
    assert np.isnan(out).all()
