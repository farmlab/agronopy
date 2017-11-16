from agronopy.crop import grain

# Constant for the Hendersen function for the maize
H_A = 0.000086541
H_B = 1.8634
H_C = 49.810


def equilibrium_moisture(t, rh, a=H_A, b=H_B, c=H_C):
    """Equilibrium Moisture Equation Constants for corn yellow dent."""
    return grain.modified_henderson_moisture(t, rh, a, b, c)


def equilibrium_relative_humidity(t, m, a=H_A, b=H_B, c=H_C):
    """Equilibrium Moisture Equation Constants for corn yellow dent."""
    return grain.modified_henderson_relative_humidity(t, m, a, b, c)
