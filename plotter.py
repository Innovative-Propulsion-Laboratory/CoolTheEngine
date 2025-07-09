import matplotlib.pyplot as plt
import cte_tools as t
import matplotlib.backends.backend_pdf
from solver2D import carto2D
from volume3d import carto3d
from tqdm import tqdm
import numpy as np


def plotter(parameters, data):
    """Plotting function for the engine cooling analysis."""

    # Disable warning about too many figure objects in memory.
    plt.rcParams.update({'figure.max_open_warning': 0})

    # List containing all figure objects
    figs = []

    # Unpack parameters
    plot_detail, show_2D_temperature, do_final_3d_plot, figure_dpi, show, save_plots = parameters

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
        molFracH2O, molFracCO2, \
        P_H2O_list, P_CO2_list, \
        q_rad_list_CO2, q_rad_list_H2O, \
        q_rad_list, q_tot_list, \
        coolant_velocity_list, \
        coolant_reynolds_list, \
        coolant_density_list, \
        coolant_cond_list, \
        coolant_cp_list, \
        coolant_visc_list, \
        coolant_temp_list, \
        coolant_pressure_list, \
        hotwall_temp_list, \
        coldwall_temp_list, \
        wall_cond_list, \
        wall_material, \
        wall_thickness = data

    # Plot of the profile of the engine
    figs.append(t.one_plot(z_coord_list_mm, r_coord_list_mm,
                           ylabel='Radius [mm]',
                           xlabel='z-coordinate [mm]',
                           title='Profile of the engine',
                           equal_axes=True, ymin=0, dpi=figure_dpi, show=show))

    # Plots of the cross-sectionnal areas
    figs.append(t.one_plot(z_coord_list_mm, cross_section_area_list,
                           title='Cross-sectional area inside the engine',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'Area [$m^2$]', ymin=0, dpi=figure_dpi, show=show))

    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[initial_fin_thickness * 1000,
                                  effective_fin_thickness_list * 1000],
                          y_label_list=['Initial fin thickness',
                                        "Effective fin thickness"],
                          colors_list=['b', 'r'],
                          title='Fin dimensions',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'Length [mm]', ymin=0, dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, channel_ar_list,
                           title=r'Channel aspect ratio',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'Aspect Ratio [-]',
                           ymin=0, dpi=figure_dpi, show=show))

    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[channel_height_list * 1000,
                                  channel_width_list * 1000],
                          y_label_list=['Channel height', "Channel width"],
                          colors_list=['b', 'r'],
                          title='Cooling channels dimensions',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'Length [mm]', ymin=0, dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, hydraulic_diameter * 1000,
                           title=r'Channel hydraulic diameter',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'$D_{hy}$ [$mm$]',
                           ymin=0, dpi=figure_dpi, show=show))

    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[initial_channel_cross_section * 1e6,
                                  effective_channel_cross_section * 1e6],
                          y_label_list=['Initial channel cross-section',
                                        "Effective channel cross-section"],
                          colors_list=['b', 'r'], sci_notation=True,
                          title='Channel cross-section area',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'Area [$mm^2$]', ymin=0, dpi=figure_dpi, show=show))

    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[alpha_list, beta_list],
                          y_label_list=['Alpha angle', "Beta angle"],
                          colors_list=['b', 'r'],
                          title='Cooling channels angles',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'Angle [°]',
                          ymin=min(-1, min(np.min(alpha_list), np.min(beta_list))),
                          ymax=max(1, max(np.max(alpha_list), np.max(beta_list))),
                          dpi=figure_dpi, show=show))

    figs.append(t.channel_around_engine_plot(x_center_list,
                                             y_center_list,
                                             z_coord_list,
                                             r_coord_list,
                                             dpi=figure_dpi,
                                             show=show))

    # Plot of the gamma linearisation
    figs.append(t.one_plot(z_coord_list_mm, gamma_list,
                           title=r'Adiabatic constant $\gamma$ of the combustion gases',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'$\gamma$ [-]', dpi=figure_dpi, show=show))

    # Plot of the Mach number in the engine (2D)
    figs.append(t.one_plot(z_coord_list_mm, mach_list,
                           title=r'Mach number',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'$Ma$ [-]', ymin=0, dpi=figure_dpi, show=show))

    # Plot of the static pressure (2D)
    figs.append(t.one_plot(z_coord_list_mm, static_pressure_list,
                           title=r'Static pressure',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'$P$ [Pa]', ymin=0, dpi=figure_dpi, show=show))

    # Plot of the temperatures in the engine
    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[hotgas_static_temp_list, hotgas_recovery_temp_list, hotgas_total_temp_list],
                          y_label_list=[r'Static temperature $T_s$', r'Recovery temperature $T_r$',
                                        r'Total temperature $T_{tot}$'],
                          colors_list=['r', 'b', 'k'],
                          title=r'Combustion gases temperature',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'$T$ [K]', dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, hotgas_recovery_temp_list,
                           title=r'Recovery temperature $T_{aw}$',
                           xlabel=r'Engine axis [$mm$]', fmt='b',
                           ylabel=r'$T$ [K]',
                           ymin=min(hotgas_recovery_temp_list) - 30,
                           ymax=max(hotgas_recovery_temp_list) + 30,
                           dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, hotgas_visc_list,
                           title=r'Hot gas dynamic viscosity',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'$\mu$ [$\mu Pa\cdot s]',
                           dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, hotgas_cp_list,
                           title=r'Hot gas $c_p$',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'$c_p$ [$\frac{J}{K\cdot kg}$]',
                           dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, hotgas_cond_list,
                           title=r'Hot gas conductivity',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'$\lambda$ [$\frac{W}{m \cdot K}$]',
                           dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, hotgas_pr_list,
                           title=r'Hot gas Prandtl number',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'$Pr_g$ [-]',
                           dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, hg_list,
                           title=r'Hot-side convection coefficient $h_g$',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'$P$ [Pa]', dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, sigma_list,
                           title=r'Bartz equation coefficient $\sigma$',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'$\sigma$ [-]',
                           dpi=figure_dpi, show=show))

    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[hl_normal_list, hl_corrected_list],
                          y_label_list=['No correction', 'Correction (Popp & Schmidt)'],
                          colors_list=['b', 'r'],
                          title=r'Cold-side convective coefficient $h_l$',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'$h_l$ [$\frac{W}{m^2 \cdot K}$]',
                          ymin=0, dpi=figure_dpi, show=show))

    # Plot of molar fraction
    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[molFracH2O, molFracCO2],
                          y_label_list=[r'$H_2O$', r'$CO_2$'],
                          colors_list=['r', 'b'],
                          title=r'Molar fraction of combustion products',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'Molar fraction $x_i$ [-]',
                          ymin=0, dpi=figure_dpi, show=show))

    # Plot of partial pressure
    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[P_H2O_list, P_CO2_list],
                          y_label_list=[r'$H_2O$', r'$CO_2$'],
                          colors_list=['r', 'b'],
                          title=r'Partial pressure of combustion products',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'Partial pressure $p_i$ [Pa]',
                          ymin=0, dpi=figure_dpi, show=show, sci_notation=True))

    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[q_rad_list_CO2, q_rad_list_H2O, q_rad_list],
                          y_label_list=['CO2', 'H2O', 'Total'],
                          colors_list=['b', 'r', 'k'],
                          title=r'Radiative heat flux',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'$\dot q_{rad}$ [$\frac{W}{m^2}$]',
                          ymin=0, dpi=figure_dpi, show=show))

    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[q_tot_list, q_rad_list],
                          y_label_list=['Total heat flux', 'Radiative Heat Flux'],
                          colors_list=['r', 'k'],
                          title=r'Total and radiative heat flux',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'$\dot q$ [$\frac{W}{m^2}$]',
                          dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, coolant_temp_list,
                           title=r'Coolant temperature',
                           xlabel=r'Engine axis [$mm$]', fmt='-b',
                           ylabel=r'$T$ [K]', dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, coolant_pressure_list,
                           title=r'Coolant pressure', fmt='-b',
                           xlabel=r'Engine axis [$mm$]', sci_notation=True,
                           ylabel=r'$P$ [Pa]', dpi=figure_dpi, show=show))

    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[coolant_velocity_list],
                          y_label_list=['Coolant'],
                          colors_list=['b'],
                          title=r'Coolant velocity',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'$V_l$ [$m/s]', dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, coolant_density_list,
                           title=r'Coolant density $\rho$',
                           xlabel=r'Engine axis [$mm$]', fmt='-b',
                           ylabel=r'$\rho$ [$\frac{kg}{m^3}$]',
                           dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, coolant_reynolds_list,
                           title=r'Coolant Reynolds number',
                           xlabel=r'Engine axis [$mm$]', fmt='-b',
                           ylabel=r'$Re_l$ [-]',
                           dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, coolant_cond_list,
                           title=r'Coolant conductivity',
                           xlabel=r'Engine axis [$mm$]', fmt='-b',
                           ylabel=r'$\lambda_l$ [$\frac{W}{m \cdot K}$]',
                           dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, coolant_cp_list,
                           title=r'Coolant $c_p$',
                           xlabel=r'Engine axis [$mm$]', fmt='-b',
                           ylabel=r'$c_p$ [$\frac{J}{K\cdot kg}$]',
                           dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, coolant_visc_list*1e6,
                           title=r'Coolant dynamic viscosity',
                           xlabel=r'Engine axis [$mm$]', fmt='-b',
                           ylabel=r'$\mu$ [$\mu Pa\cdot s$]',
                           dpi=figure_dpi, show=show))

    figs.append(t.n_plots(z_coord_list_mm,
                          y_list=[hotwall_temp_list, coldwall_temp_list],
                          y_label_list=['Hot wall', 'Cold wall'],
                          colors_list=['r', 'b'],
                          title=r'Wall temperatures $T_{wg}$ and $T_{wl}$',
                          xlabel=r'Engine axis [$mm$]',
                          ylabel=r'$T$ [K]', dpi=figure_dpi, show=show))

    figs.append(t.one_plot(z_coord_list_mm, wall_cond_list,
                           title=f'Wall conductivity ({wall_material})',
                           xlabel=r'Engine axis [$mm$]',
                           ylabel=r'$\lambda_w$ [$\frac{W}{m \cdot K}$]',
                           dpi=figure_dpi, show=show))

    if show_2D_temperature and show:
        # At the beginning of the chamber
        print("█ 2D Results at the beginning of the chamber :                             █")
        dx = 5e-5  # m
        location = " at the beginning of the chamber"
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
                10, True, 2, location, False)
        # At the end of the divergent
        print("█ 2D Results at the manifold :                                             █")
        dx = 5e-5  # m
        location = " at the manifold"
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

    "Computation for 3D graph"
    if do_final_3d_plot and show:
        nb_points = len(z_coord_list_mm)
        temperature_slice_list = []
        dx = 0.0001

        # Compute a (low-resolution) 2D slice for each point in the engine
        with tqdm(total=nb_points,
                  desc="█ 3D graph computation         ",
                  unit="|   █", bar_format="{l_bar}{bar}{unit}",
                  ncols=76) as progressbar:
            for i in range(0, nb_points):
                temperature_slice = carto2D(effective_fin_thickness_list[i] + channel_width_list[i], channel_width_list[i], wall_thickness,
                                            channel_height_list[i], dx, hg_list[i], wall_cond_list[i], hotgas_recovery_temp_list[i],
                                            hl_corrected_list[i], coolant_temp_list[i], 3, False, 1, "", True)
                temperature_slice_list.append(temperature_slice)
                progressbar.update(1)

        # Stack all these slices in a final 3D plot
        carto3d([0, 0, 0], z_coord_list, r_coord_list, temperature_slice_list, plt.cm.Spectral_r,
                '3D view of wall temperatures (in K)', nb_channels, 0.05)
        print("█                                                                          █")

    if plot_detail >= 1 and save_plots:
        pdf = matplotlib.backends.backend_pdf.PdfPages("output/graphs.pdf")
        for fig in figs:
            fig.savefig(pdf, format='pdf')
        pdf.close()
