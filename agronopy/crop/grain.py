import math


def modified_henderson_moisture(t, rh, a, b, c):
    """Equilibrium Moisture Equation Constants.

    Empirical equilibrium moisture content equation: modified Henderson
    equation. This equation can be used to calculate equilibrium moisture
    content or equilibrium relative humidity with the constants below for
    typical products.
    """
    return 1 / 100 * (-math.log(1 - rh) / (a * (t + c)))**(1 / b)


def modified_henderson_relative_humidity(t, m, a, b, c):
    """Equilibrium Moisture Equation Constants.

    Empirical equilibrium moisture content equations are the modified Henderson
    equation and the Chung-Pfost equation. These equations can be used to
    calculate equilibrium moisture content or equilibrium relative humidity
    with the constants below for typical products.
    """
    return 1 - math.exp(-a * (t + c) * (100 * m)**b)
