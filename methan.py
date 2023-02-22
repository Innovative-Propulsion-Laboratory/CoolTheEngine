from CoolProp.CoolProp import PropsSI


def densityCH4(P, T, fluid):
    return PropsSI("D", "T", T, "P", P, fluid)


def cpCH4(P, T, fluid):
    return PropsSI("C", "P", P, "T", T, fluid)


def condCH4(P, T, fluid):
    return PropsSI("L", "T", T, "P", P, fluid)


def viscCH4(P, T, fluid):
    return PropsSI("V", "T", T, "P", P, fluid)


def sound_speed_CH4(P, T, fluid):
    return PropsSI("A", "T", T, "P", P, fluid)


def entropyCH4(P, T, fluid):
    return PropsSI("S", "T", T, "P", P, fluid)


def pressureCH4(S, T, fluid):
    return PropsSI('P', 'T', T, 'S', S, fluid)


def DeltaT(P, T, fluid):
    return PropsSI("T", "P", P, "Q", 0, fluid) - T
