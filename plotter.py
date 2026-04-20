import matplotlib.pyplot as plt
import cte_tools as t
import matplotlib.backends.backend_pdf
from pathlib import Path
import textwrap
from solver2D import carto2D
from tqdm import tqdm
import numpy as np


def append_run_log_figure(figures: list, log_text: str) -> list:
    """Return a new list with data-summary page(s) inserted as the second figure(s)."""
    if not log_text.strip():
        return list(figures)

    raw_lines = [line.strip() for line in log_text.strip().splitlines() if line.strip()]
    cleaned_lines: list[str] = []
    for line in raw_lines:
        cleaned = line.replace('█', ' ').strip()
        cleaned = ' '.join(cleaned.split())
        if cleaned:
            cleaned_lines.append(cleaned)

    section_patterns = [
        ("2d results at the injection plate", "2D Results at the injection plate:"),
        ("2d results at the throat", "2D Results at the throat:"),
        ("2d results at the nozzle exit", "2D Results at the nozzle exit:"),
    ]
    detail_prefixes = (
        "mean wall temperature at hot gas side",
        "mean wall temperature at coolant side",
        "maximum temperature in the wall",
    )
    global_prefixes = (
        "maximum wall temperature",
        "average wall temperature",
        "maximum wall stress",
        "coolant pressure drop",
        "coolant temperature increase",
    )

    section_data: dict[str, list[str]] = {display: [] for _, display in section_patterns}
    global_data: list[str] = []
    current_section: str | None = None

    for line in cleaned_lines:
        lower_line = line.lower()

        matched_section = None
        for pattern, display in section_patterns:
            if pattern in lower_line:
                matched_section = display
                break
        if matched_section is not None:
            current_section = matched_section
            continue

        if current_section is not None and any(lower_line.startswith(prefix) for prefix in detail_prefixes):
            if "=" in line:
                key, value = line.split("=", 1)
                section_data[current_section].append(f"{key.strip()} = {value.strip()}")
            else:
                section_data[current_section].append(line)
            continue

        if any(lower_line.startswith(prefix) for prefix in global_prefixes):
            if ":" in line:
                key, value = line.split(":", 1)
                global_data.append(f"{key.strip()}: {value.strip()}")
            else:
                global_data.append(line)

    render_rows: list[tuple[str, str]] = []
    for _, display in section_patterns:
        details = section_data[display]
        if len(details) == 0:
            continue
        render_rows.append((display, "header"))
        for detail_line in details:
            render_rows.append((detail_line, "detail"))
        render_rows.append(("", "spacer"))

    if len(global_data) > 0:
        render_rows.append(("Global metrics:", "header"))
        for line in global_data:
            render_rows.append((line, "detail"))

    if len(render_rows) == 0:
        return list(figures)

    lines_per_page = 22
    log_pages: list = []
    for page_start in range(0, len(render_rows), lines_per_page):
        page_rows = render_rows[page_start:page_start + lines_per_page]

        summary_fig = plt.figure(figsize=(6.4, 4.8), dpi=150)
        summary_fig.text(0.5, 0.96, 'CoolTheEngine - Run Data Summary',
                         ha='center', va='top', fontsize=14, weight='bold')

        if len(render_rows) > lines_per_page:
            page_idx = page_start // lines_per_page + 1
            page_total = (len(render_rows) + lines_per_page - 1) // lines_per_page
            summary_fig.text(0.92, 0.96, f'Page {page_idx}/{page_total}', ha='right', va='top', fontsize=9)

        y_pos = 0.88
        line_height = 0.03
        left_margin = 0.08

        for line, row_kind in page_rows:
            if row_kind == "header":
                summary_fig.text(left_margin, y_pos, line, fontsize=10, weight='bold')
            elif row_kind == "detail":
                summary_fig.text(left_margin + 0.04, y_pos, line, fontsize=9)
            else:
                # spacer row
                pass
            y_pos -= line_height

        plt.axis('off')
        log_pages.append(summary_fig)

    merged = list(figures)
    insert_at = 1 if len(merged) >= 1 else 0
    for i, log_fig in enumerate(log_pages):
        merged.insert(insert_at + i, log_fig)
    return merged


def plotter(parameters, data, output_dir: str | Path = "output"):
    """Plotting function for the engine cooling analysis."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Disable warning about too many figure objects in memory.
    plt.rcParams.update({'figure.max_open_warning': 0})

    # List containing all figure objects
    figs = []

    # Unpack parameters
    show_1D, show_2D, figure_dpi, save_plots, cte_config = parameters
    # Create configuration summary page (first page of PDF)
    if save_plots:
        config_fig = plt.figure(figsize=(6.4, 4.8), dpi=figure_dpi)
        config_fig.text(0.5, 0.96, 'CoolTheEngine - Configuration Summary',
                        ha='center', va='top', fontsize=14, weight='bold')

        y_pos = 0.88
        line_height = 0.035
        section_gap = 0.05
        left_margin = 0.08
        fontsize_section = 10
        fontsize_item = 9

        # Engine Parameters
        config_fig.text(left_margin, y_pos, 'Engine Parameters',
                        fontsize=fontsize_section, weight='bold')
        y_pos -= line_height * 1.2
        config_fig.text(left_margin + 0.02, y_pos, f'Chamber Pressure: {cte_config.chamber_pressure/1e5:.2f} bar',
                        fontsize=fontsize_item)
        y_pos -= line_height
        config_fig.text(left_margin + 0.02, y_pos, f'Oxidizer: {cte_config.ox_name}',
                        fontsize=fontsize_item)
        y_pos -= line_height
        config_fig.text(left_margin + 0.02, y_pos, f'Oxidizer Mass Flow Rate: {cte_config.ox_mfr:.3f} kg/s',
                        fontsize=fontsize_item)
        y_pos -= line_height
        config_fig.text(left_margin + 0.02, y_pos, f'Fuel: {cte_config.fuel_name}',
                        fontsize=fontsize_item)
        y_pos -= line_height
        config_fig.text(left_margin + 0.02, y_pos, f'Fuel Mass Flow Rate: {cte_config.fuel_mfr:.3f} kg/s',
                        fontsize=fontsize_item)
        y_pos -= section_gap

        # Coolant Properties
        config_fig.text(left_margin, y_pos, 'Coolant Properties',
                        fontsize=fontsize_section, weight='bold')
        y_pos -= line_height * 1.2
        config_fig.text(left_margin + 0.02, y_pos, f'Coolant: {cte_config.coolant_name}',
                        fontsize=fontsize_item)
        y_pos -= line_height
        config_fig.text(left_margin + 0.02, y_pos, f'Coolant Inlet Temperature: {cte_config.coolant_inlet_temp:.1f} K',
                        fontsize=fontsize_item)
        y_pos -= line_height
        config_fig.text(left_margin + 0.02, y_pos, f'Coolant Inlet Pressure: {cte_config.coolant_inlet_pressure/1e5:.2f} bar',
                        fontsize=fontsize_item)
        y_pos -= line_height
        config_fig.text(left_margin + 0.02, y_pos, f'Coolant Mass Flow Rate: {cte_config.coolant_mfr:.3f} kg/s',
                        fontsize=fontsize_item)
        y_pos -= section_gap

        # Wall Material Properties
        config_fig.text(left_margin, y_pos, 'Wall Material Properties',
                        fontsize=fontsize_section, weight='bold')
        y_pos -= line_height * 1.2
        config_fig.text(left_margin + 0.02, y_pos, f'Wall Material: {cte_config.wall_material}',
                        fontsize=fontsize_item)
        y_pos -= line_height
        config_fig.text(left_margin + 0.02, y_pos, f'TBC Thickness: {cte_config.TBC_thickness*1e6:.1f} µm',
                        fontsize=fontsize_item)
        y_pos -= line_height
        config_fig.text(left_margin + 0.02, y_pos, f'Channel Roughness: {cte_config.channel_roughness*1e6:.1f} µm',
                        fontsize=fontsize_item)
        y_pos -= line_height
        config_fig.text(left_margin + 0.02, y_pos, f'Wall Thickness: {cte_config.wall_thickness*1000:.2f} mm',
                        fontsize=fontsize_item)
        y_pos -= line_height
        config_fig.text(left_margin + 0.02, y_pos, f'Total Wall Thickness: {cte_config.total_wall_thickness*1000:.2f} mm',
                        fontsize=fontsize_item)

        # Channel Configuration - Second Column (Right side)
        right_margin = 0.54
        y_pos_right = 0.88

        config_fig.text(right_margin, y_pos_right, 'Channel Configuration',
                        fontsize=fontsize_section, weight='bold')
        y_pos_right -= line_height * 1.2
        config_fig.text(right_margin + 0.02, y_pos_right, f'Number of Channels: {cte_config.nb_channels}',
                        fontsize=fontsize_item)
        y_pos_right -= section_gap

        # Channel Widths
        config_fig.text(right_margin + 0.02, y_pos_right, 'Channel Widths:', weight='bold',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Injection: {cte_config.width_inj*1000:.2f} mm',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Converging: {cte_config.width_conv*1000:.2f} mm',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Throat: {cte_config.width_throat*1000:.2f} mm',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Exit: {cte_config.width_exit*1000:.2f} mm',
                        fontsize=fontsize_item)
        y_pos_right -= section_gap

        # Channel Heights
        config_fig.text(right_margin + 0.02, y_pos_right, 'Channel Heights:', weight='bold',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Injection: {cte_config.ht_inj*1000:.2f} mm',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Converging: {cte_config.ht_conv*1000:.2f} mm',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Throat: {cte_config.ht_throat*1000:.2f} mm',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Exit: {cte_config.ht_exit*1000:.2f} mm',
                        fontsize=fontsize_item)
        y_pos_right -= section_gap

        # Channel Angles
        config_fig.text(right_margin + 0.02, y_pos_right, 'Channel Angles:', weight='bold',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Injection: {cte_config.beta_inj:.1f}°',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Converging: {cte_config.beta_conv:.1f}°',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Throat: {cte_config.beta_throat:.1f}°',
                        fontsize=fontsize_item)
        y_pos_right -= line_height
        config_fig.text(right_margin + 0.04, y_pos_right, f'Exit: {cte_config.beta_exit:.1f}°',
                        fontsize=fontsize_item)

        plt.axis('off')
        figs.insert(0, config_fig)  # Insert as first page
    # Unpack data (new variable names)
    z_coord_list_mm, r_coord_list_mm, \
        z_coord_list, r_coord_list, \
        cross_section_area_list, \
        initial_fin_thickness, \
        effective_fin_thickness_list, \
        channel_ar_list, \
        channel_width_list, \
        channel_height_list, \
        hydraulic_diameter, \
        initial_channel_cross_section, \
        effective_channel_cross_section, \
        nb_channels, \
        alpha_list, \
        beta_list, \
        x_center_list, \
        y_center_list, \
        gamma_list, \
        mach_list, \
        static_pressure_list, \
        hotgas_total_temp_list, \
        hotgas_recovery_temp_list, \
        hotgas_static_temp_list, \
        hotgas_visc_list, \
        hotgas_cp_list, \
        hotgas_cond_list, \
        hotgas_pr_list, \
        hg_list, \
        sigma_list, \
        hl_normal_list, \
        hl_corrected_list, \
        h_tp_list, \
        molFracH2O, molFracCO2, \
        P_H2O_list, P_CO2_list, \
        q_rad_list_CO2, q_rad_list_H2O, \
        q_rad_list, q_tot_list, \
        CHF_Meyer_list, \
        CHF_Tong_list, \
        coolant_velocity_list, \
        coolant_reynolds_list, \
        coolant_density_list, \
        coolant_cond_list, \
        coolant_cp_list, \
        coolant_visc_list, \
        coolant_temp_list, \
        coolant_pressure_list, \
        coolant_Tsat_list, \
        tbc_temp_list, \
        hotwall_temp_list, \
        coldwall_temp_list, \
        wall_cond_list, \
        wall_material, \
        wall_thickness, \
        hoop_stress_list, \
        thermal_stress_list, \
        max_wall_stress_list = data

    if show_1D:
        material_temp_limit = np.ones_like(z_coord_list_mm) * t.wall_temp_limit(wall_material)
        yield_strength_list = np.array([t.wall_yield_strength(hotwall_temp_list[i], coldwall_temp_list[i], wall_material) for i in range(len(z_coord_list))])

        #  of the profile of the engine
        figs.append(t.one_plot(z_coord_list_mm, r_coord_list_mm,
                               ylabel='Radius [mm]',
                               xlabel='z-coordinate [mm]',
                               title='Profile of the engine',
                               equal_axes=True, ymin=0, dpi=figure_dpi, show=show_1D))

        # Plots of the cross-sectionnal areas
        figs.append(t.one_plot(z_coord_list_mm, cross_section_area_list,
                               title='Cross-sectional area inside the engine',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'Area [$m^2$]', ymin=0, dpi=figure_dpi, show=show_1D))

        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[initial_fin_thickness * 1000,
                                      effective_fin_thickness_list * 1000],
                              y_label_list=['Initial fin thickness',
                                            "Effective fin thickness"],
                              colors_list=['b', 'r'],
                              title='Fin dimensions',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'Length [mm]', ymin=0, dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, channel_ar_list,
                               title=r'Channel aspect ratio',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'Aspect Ratio [-]',
                               ymin=0, dpi=figure_dpi, show=show_1D))

        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[channel_height_list * 1000,
                                      channel_width_list * 1000],
                              y_label_list=['Channel height', "Channel width"],
                              colors_list=['b', 'r'],
                              title='Cooling channels dimensions',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'Length [mm]', ymin=0, dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, hydraulic_diameter * 1000,
                               title=r'Channel hydraulic diameter',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'$D_{hy}$ [$mm$]',
                               ymin=0, dpi=figure_dpi, show=show_1D))

        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[initial_channel_cross_section * 1e6,
                                      effective_channel_cross_section * 1e6],
                              y_label_list=['Initial channel cross-section',
                                            "Effective channel cross-section"],
                              colors_list=['b', 'r'], sci_notation=True,
                              title='Channel cross-section area',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'Area [$mm^2$]', ymin=0, dpi=figure_dpi, show=show_1D))

        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[alpha_list, beta_list],
                              y_label_list=['Alpha angle', "Beta angle"],
                              colors_list=['b', 'r'],
                              title='Cooling channels angles',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'Angle [°]',
                              ymin=min(-1, min(np.min(alpha_list), np.min(beta_list))),
                              ymax=max(1, max(np.max(alpha_list), np.max(beta_list))),
                              dpi=figure_dpi, show=show_1D))

        figs.append(t.channel_around_engine_plot(x_center_list,
                                                 y_center_list,
                                                 z_coord_list,
                                                 r_coord_list,
                                                 dpi=figure_dpi,
                                                 show=show_1D))

        # Plot of the gamma linearisation
        figs.append(t.one_plot(z_coord_list_mm, gamma_list,
                               title=r'Adiabatic constant $\gamma$ of the combustion gases',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'$\gamma$ [-]', dpi=figure_dpi, show=show_1D))

        # Plot of the Mach number in the engine (2D)
        figs.append(t.one_plot(z_coord_list_mm, mach_list,
                               title=r'Mach number',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'$Ma$ [-]', ymin=0, dpi=figure_dpi, show=show_1D))

        # Plot of the static pressure (2D)
        figs.append(t.one_plot(z_coord_list_mm, static_pressure_list,
                               title=r'Static pressure',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'$P$ [Pa]', ymin=0, dpi=figure_dpi, show=show_1D))

        # Plot of the temperatures in the engine
        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[hotgas_static_temp_list, hotgas_recovery_temp_list, hotgas_total_temp_list],
                              y_label_list=[r'Static temperature $T_s$', r'Recovery temperature $T_r$',
                                            r'Total temperature $T_{tot}$'],
                              colors_list=['r', 'b', 'k'],
                              title=r'Combustion gases temperature',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'$T$ [K]', dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, hotgas_recovery_temp_list,
                               title=r'Recovery temperature $T_{aw}$',
                               xlabel=r'Engine axis [$mm$]', fmt='b',
                               ylabel=r'$T$ [K]',
                               ymin=min(hotgas_recovery_temp_list) - 30,
                               ymax=max(hotgas_recovery_temp_list) + 30,
                               dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, hotgas_visc_list,
                               title=r'Hot gas dynamic viscosity',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'$\mu$ [$\mu Pa\cdot s]',
                               dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, hotgas_cp_list,
                               title=r'Hot gas $c_p$',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'$c_p$ [$\frac{J}{K\cdot kg}$]',
                               dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, hotgas_cond_list,
                               title=r'Hot gas conductivity',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'$\lambda$ [$\frac{W}{m \cdot K}$]',
                               dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, hotgas_pr_list,
                               title=r'Hot gas Prandtl number',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'$Pr_g$ [-]',
                               dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, hg_list,
                               title=r'Hot-side convection coefficient $h_g$',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'$h_g$ [$\frac{W}{m^2 \cdot K}$]', dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, sigma_list,
                               title=r'Bartz equation coefficient $\sigma$',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'$\sigma$ [-]',
                               dpi=figure_dpi, show=show_1D))

        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[hl_normal_list, hl_corrected_list, h_tp_list],
                              y_label_list=['Gnielinski', 'Fin effect correction (Popp & Schmidt)', 'Subcooled flow boiling (Lui & Winterton)'],
                              colors_list=['b', 'r', 'g'],
                              title=r'Cold-side convective coefficient $h_l$/$h_{tp}$',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'$h$ [$\frac{W}{m^2 \cdot K}$]',
                              ymin=0, dpi=figure_dpi, show=show_1D))

        # Plot of molar fraction
        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[molFracH2O, molFracCO2],
                              y_label_list=[r'$H_2O$', r'$CO_2$'],
                              colors_list=['r', 'b'],
                              title=r'Molar fraction of combustion products',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'Molar fraction $x_i$ [-]',
                              ymin=0, dpi=figure_dpi, show=show_1D))

        # Plot of partial pressure
        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[P_H2O_list, P_CO2_list],
                              y_label_list=[r'$H_2O$', r'$CO_2$'],
                              colors_list=['r', 'b'],
                              title=r'Partial pressure of combustion products',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'Partial pressure $p_i$ [Pa]',
                              ymin=0, dpi=figure_dpi, show=show_1D, sci_notation=True))

        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[q_rad_list_CO2, q_rad_list_H2O, q_rad_list],
                              y_label_list=['CO2', 'H2O', 'Total'],
                              colors_list=['b', 'r', 'k'],
                              title=r'Radiative heat flux',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'$\dot q_{rad}$ [$\frac{W}{m^2}$]',
                              ymin=0, dpi=figure_dpi, show=show_1D))

        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[q_tot_list, q_rad_list, CHF_Meyer_list, CHF_Tong_list],
                              y_label_list=['Total heat flux', 'Radiative Heat Flux', 'CHF Meyer', 'CHF Tong'],
                              colors_list=['r', 'k', 'g', 'b'],
                              title=r'Heat fluxes',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'$\dot q$ [$\frac{W}{m^2}$]',
                              dpi=figure_dpi, show=show_1D))

        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[coolant_temp_list, coolant_Tsat_list],
                              y_label_list=['Coolant temperature', 'Coolant saturation temperature'],
                              colors_list=['r', 'k'],
                              title=r'Coolant heating',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'$T$ [K]',
                              dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, coolant_pressure_list,
                               title=r'Coolant pressure', fmt='-b',
                               xlabel=r'Engine axis [$mm$]', sci_notation=True,
                               ylabel=r'$P$ [Pa]', dpi=figure_dpi, show=show_1D))

        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[coolant_velocity_list],
                              y_label_list=['Coolant'], fmt='-b',
                              colors_list=['b'],
                              title=r'Coolant velocity',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'$V_l$ [$m/s]', dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, coolant_density_list,
                               title=r'Coolant density $\rho$',
                               xlabel=r'Engine axis [$mm$]', fmt='-b',
                               ylabel=r'$\rho$ [$\frac{kg}{m^3}$]',
                               dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, coolant_reynolds_list,
                               title=r'Coolant Reynolds number',
                               xlabel=r'Engine axis [$mm$]', fmt='-b',
                               ylabel=r'$Re_l$ [-]',
                               dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, coolant_cond_list,
                               title=r'Coolant conductivity',
                               xlabel=r'Engine axis [$mm$]', fmt='-b',
                               ylabel=r'$\lambda_l$ [$\frac{W}{m \cdot K}$]',
                               dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, coolant_cp_list,
                               title=r'Coolant $c_p$',
                               xlabel=r'Engine axis [$mm$]', fmt='-b',
                               ylabel=r'$c_p$ [$\frac{J}{K\cdot kg}$]',
                               dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, coolant_visc_list*1e6,
                               title=r'Coolant dynamic viscosity',
                               xlabel=r'Engine axis [$mm$]', fmt='-b',
                               ylabel=r'$\mu$ [$\mu Pa\cdot s$]',
                               dpi=figure_dpi, show=show_1D))

        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[tbc_temp_list, hotwall_temp_list, coldwall_temp_list, material_temp_limit],
                              y_label_list=['TBC temperature', 'Hot wall', 'Cold wall', "Maximum allowable temperature"],
                              colors_list=['m', 'r', 'b', 'k'],
                              title=r'Wall temperatures $T_{tbc}$, $T_{wg}$ and $T_{wl}$',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'$T$ [K]', dpi=figure_dpi, show=show_1D))

        figs.append(t.n_plots(z_coord_list_mm,
                              y_list=[hoop_stress_list/1e6, thermal_stress_list/1e6, max_wall_stress_list/1e6, yield_strength_list/1e6, -yield_strength_list/1e6],
                              y_label_list=['Hoop stress', 'Thermal stress', "Total wall stress", "Yield strength (tensile)", "Yield strength (compressive)"],
                              colors_list=['g', 'b', "r", 'k', 'k'],
                              title=r'Stress at the inner side of the hotwall (stresses are positive in tension)',
                              xlabel=r'Engine axis [$mm$]',
                              ylabel=r'$\sigma$ [MPa]', dpi=figure_dpi, show=show_1D))

        figs.append(t.one_plot(z_coord_list_mm, wall_cond_list,
                               title=f'Wall conductivity ({wall_material})',
                               xlabel=r'Engine axis [$mm$]',
                               ylabel=r'$\lambda_w$ [$\frac{W}{m \cdot K}$]',
                               dpi=figure_dpi, show=show_1D))

    if show_2D:
        # At the beginning of the chamber
        print("█ 2D Results at the injection plate :                                      █")
        dx = 5e-5  # m
        location = " at the injection plate"
        carto2D(effective_fin_thickness_list[0] + channel_width_list[0],
                channel_width_list[0],
                wall_thickness,
                channel_height_list[0],
                dx,
                hg_list[0],
                wall_cond_list[0],
                hotgas_recovery_temp_list[0],
                hl_corrected_list[0],
                coolant_temp_list[0],
                10, True, 1, location, False)

        # At the throat
        print("█ 2D Results at the throat :                                               █")
        i_throat = np.argmin(r_coord_list)
        dx = 5e-5  # m
        location = " at the throat"
        carto2D(effective_fin_thickness_list[i_throat] + channel_width_list[i_throat],
                channel_width_list[i_throat],
                wall_thickness,
                channel_height_list[i_throat],
                dx,
                hg_list[i_throat],
                wall_cond_list[i_throat],
                hotgas_recovery_temp_list[i_throat],
                hl_corrected_list[i_throat],
                coolant_temp_list[i_throat],
                20, True, 2, location, False)
        # At the end of the divergent
        print("█ 2D Results at the nozzle exit :                                          █")
        dx = 5e-5  # m
        location = " at the nozzle exit"
        carto2D(effective_fin_thickness_list[-1] + channel_width_list[-1],
                channel_width_list[-1],
                wall_thickness,
                channel_height_list[-1],
                dx,
                hg_list[-1],
                wall_cond_list[-1],
                hotgas_recovery_temp_list[-1],
                hl_corrected_list[-1],
                coolant_temp_list[-1],
                10, True, 1, location, False)

    if save_plots:
        pdf = matplotlib.backends.backend_pdf.PdfPages(output_dir / "graphs.pdf")
        for fig in figs:
            fig.savefig(pdf, format='pdf')
        pdf.close()

    return figs


