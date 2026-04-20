from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from config_enums import CoolantName, FuelName, OxidizerName, WallMaterial


@dataclass
class FieldSpec:
    key: str
    label: str
    group: str
    expected: str  # "float", "int", "enum"
    unit_hint: str = ""
    subgroup: str | None = None
    to_si_factor: float = 1.0
    editable: bool = True
    allow_list: bool = False
    minimum: float | None = None
    maximum: float | None = None
    enum_cls: type | None = None
    default: Any = None


FIELD_SPECS: list[FieldSpec] = [
    # Engine parameters
    FieldSpec("chamber_pressure", "Chamber pressure Pc [bar]", "Engine parameters", "float", minimum=5e5, maximum=400e5, to_si_factor=1e5),
    FieldSpec("oxidizer_name", "Oxidizer name", "Engine parameters", "enum", enum_cls=OxidizerName),
    FieldSpec("fuel_name", "Fuel name", "Engine parameters", "enum", enum_cls=FuelName),
    FieldSpec("ox_mfr", "Oxidizer mass flow rate [kg/s]", "Engine parameters", "float"),
    FieldSpec("fuel_mfr", "Fuel mass flow rate [kg/s]", "Engine parameters", "float"),

    # Coolant parameters
    FieldSpec("coolant_name", "Coolant name", "Coolant parameters", "enum", unit_hint="Dropdown or list", enum_cls=CoolantName, allow_list=True),
    FieldSpec("P_coolant", "Coolant inlet pressure [bar]", "Coolant parameters", "float", unit_hint="Scalar or JSON list", allow_list=True, to_si_factor=1e5),
    FieldSpec("T_coolant", "Coolant inlet temperature [K]", "Coolant parameters", "float", unit_hint="Scalar or JSON list", allow_list=True),
    FieldSpec("coolant_mfr", "Coolant mass flow rate [kg/s]", "Coolant parameters", "float", unit_hint="Scalar or JSON list", allow_list=True),

    # Channel sizing
    FieldSpec("nb_channels", "Number of cooling channels", "Channel sizing", "int", unit_hint="Scalar or JSON list", allow_list=True),
    FieldSpec("channel_widths_inj", "At injector end [mm]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel width", allow_list=True, to_si_factor=1e-3),
    FieldSpec("channel_widths_conv", "At the end of the cylindrical section [mm]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel width", allow_list=True, to_si_factor=1e-3),
    FieldSpec("channel_widths_throat", "At the throat [mm]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel width", allow_list=True, to_si_factor=1e-3),
    FieldSpec("channel_widths_exit", "At the nozzle exit [mm]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel width", allow_list=True, to_si_factor=1e-3),
    FieldSpec("channel_heights_inj", "At injector end [mm]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel height", allow_list=True, to_si_factor=1e-3),
    FieldSpec("channel_heights_conv", "At the end of the cylindrical section [mm]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel height", allow_list=True, to_si_factor=1e-3),
    FieldSpec("channel_heights_throat", "At the throat [mm]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel height", allow_list=True, to_si_factor=1e-3),
    FieldSpec("channel_heights_exit", "At the nozzle exit [mm]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel height", allow_list=True, to_si_factor=1e-3),
    FieldSpec("channel_angles_inj", "At injector end [°]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel angle", allow_list=True),
    FieldSpec("channel_angles_conv", "At the end of the cylindrical section [°]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel angle", allow_list=True),
    FieldSpec("channel_angles_throat", "At the throat [°]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel angle", allow_list=True),
    FieldSpec("channel_angles_exit", "At the nozzle exit [°]", "Channel sizing", "float", unit_hint="Scalar or JSON list", subgroup="Channel angle", allow_list=True),

    # Wall parameters
    FieldSpec("wall_material", "Wall material", "Wall parameters", "enum", unit_hint="Dropdown or list", enum_cls=WallMaterial, allow_list=True),
    FieldSpec("tbc_thickness", "Thermal barrier coating thickness (PDMS/SiO2) [um]", "Wall parameters", "float", unit_hint="Scalar or JSON list", allow_list=True, to_si_factor=1e-6),
    FieldSpec("channel_roughness", "Channel wall roughness (equivalent sand grain roughness ~6x Ra) [um]", "Wall parameters", "float", unit_hint="Scalar or JSON list", allow_list=True, to_si_factor=1e-6),
    FieldSpec("wall_thickness", "Hot wall thickness [mm]", "Wall parameters", "float", unit_hint="Scalar or JSON list", allow_list=True, to_si_factor=1e-3),
    FieldSpec("total_wall_thickness", "Total wall thickness [mm]", "Wall parameters", "float", unit_hint="Scalar or JSON list", allow_list=True, to_si_factor=1e-3),

    # Included in JSON but hidden from UI
    FieldSpec("nb_points", "nb_points", "Hidden", "int", editable=False, default=500),
    FieldSpec("figure_dpi", "figure_dpi", "Hidden", "int", editable=False, default=150),
]
