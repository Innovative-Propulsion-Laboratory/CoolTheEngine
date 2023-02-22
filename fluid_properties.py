from CoolProp.CoolProp import PropsSI


def density(P, T, fluid):
    return PropsSI("D", "T", T, "P", P, fluid)


def cp(P, T, fluid):
    return PropsSI("C", "P", P, "T", T, fluid)


def conductivity(P, T, fluid):
    return PropsSI("L", "T", T, "P", P, fluid)


def viscosity(P, T, fluid):
    return PropsSI("V", "T", T, "P", P, fluid)


def sound_speed(P, T, fluid):
    return PropsSI("A", "T", T, "P", P, fluid)


def entropy(P, T, fluid):
    return PropsSI("S", "T", T, "P", P, fluid)


def pressure(S, T, fluid):
    return PropsSI('P', 'T', T, 'S', S, fluid)


def DeltaT(P, T, fluid):
    return PropsSI("T", "P", P, "Q", 0, fluid) - T
