"""Prefab analysis protocols for lipid bilayers.

This module defines a set of functions that perform bilayer analysis using
pre-designed (prefab) protocols that automatically setup and run the analyses,
extract and output data, and generate plots. Most of these protocols make use of
the BilayerAnalyzer class, but some use tools from the mda_tools package. These
analyses are meant for quasi-planar lipid bilayer systems.

Example:
    >> from pybilt.bilayer_analyzer import prefab_analysis_protocols as ba_pap
    >> ba_pap.msd_diffusion(structure_file='my_structure',
    ... trajectory_file='my_trajectory', selection_string='resname DOPE POPC',
    ... dump_path='./my_output_directory/')
"""

import numpy as np
import MDAnalysis as mda
try:
    import cPickle as pickle
except ImportError as error:
    import pickle
import pybilt.plot_generation.plot_generation_functions as pgf
from pybilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer
from pybilt.diffusion import diffusion_coefficients as dc
from pybilt.plot_generation.plot_generation_functions import _color_list
from pybilt.common.running_stats import BlockAverager
from pybilt.mda_tools.mda_density_map import position_density_map_2d_leaflet_simple


def msd_diffusion(structure_file, trajectory_file, selection_string,
                  resnames=None, frame_start=0, frame_end=-1, frame_interval=1,
                  dump_path=None):
    """Protocol to compute MSD curves and estimate diffusion coefficents.

    This function uses the BilayerAnalyzer with the analysis 'msd' to compute
    the time series curves for the mean squared displacements of lipids in the
    selection and then uses tools from the diffusion package to estimate the
    diffusion coefficents via three different estimators; i.e., the Einstein
    relation, linear (least squares) line fit, and a non-linear anomalous
    diffusion model fit. The protocol uses two different time ranges for the
    diffusion coefficent estimations: 0 (or first frame) to 10000.0 (assumed to
    be 10000 ps or 10 ns) and 10000.0 to 100000.0 (assumed to 10000 ps or 10 ns
    to 100000 ps or 100 ns). Reports of data are printed to standard out, while
    pickled data is dumped to disk and MSD time series plots are generated and
    also dumped to disk. The output files have the prefix 'msd'.

     Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        resnames (Optional[str, list]): Specify the resnames of the lipid types
            that are to be included in this analysis. Defaults to None, which
            includes all lipid types.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.

    Returns:
        void

    Notes:
        This function summarizes the diffusion coefficent estimates by
        printing to the std out, which could be redirected into a file.
    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # for each lipid resname in the input resnames add an msd analysis for
    # both, upper, and lower leaflets
    if resnames is None:
        resnames = []
    for lipid_type in resnames:
        a_name = "msd_"+lipid_type
        add_string = "msd "+a_name+" resname "+lipid_type
        add_string_upper = "msd " + a_name + "_upper resname " + \
                           lipid_type + " leaflet upper"
        add_string_lower = "msd " + a_name + "_lower resname " + \
                           lipid_type + " leaflet lower"
        analyzer.add_analysis(add_string)
        analyzer.add_analysis(add_string_upper)
        analyzer.add_analysis(add_string_lower)

    # now add plot protocols
    analyzer.add_plot('msd msd_a msd_1 All')
    # plot with each of the lipid types on same plot (leaflet=both)
    l_msd_string = "msd msd_l"
    for lipid_type in resnames:
        l_msd_string += " msd_"+lipid_type+" "+lipid_type
    analyzer.add_plot(l_msd_string)

    # plots for each lipid with upper, lower, both on same plot
    for lipid_type in resnames:
        add_string = "msd " + lipid_type + "_cul msd_" + lipid_type + \
                     " Composite msd_" + lipid_type + "_upper Upper msd_" + \
                     lipid_type + "_lower Lower"
        analyzer.add_plot(add_string)

    analyzer.print_analysis_protocol()
    analyzer.print_plot_protocol()

    # run the analyzer
    analyzer.run_analysis()

    # dump the output
    analyzer.dump_data(path=dump_path)

    # make MSD vs time plots
    analyzer.save_all_plots()

    # compute diffusion coefficients
    # composite

    msd_dat = analyzer.get_analysis_data('msd_1')

    times = msd_dat[:, 0]
    msd_vals = msd_dat[:, 1]
    t0 = times[0]
    t10 = t0 + 10000.0
    t100 = t0 + 100000.0
    # diffusion coeff, whole range
    # Simple application of Einstein relation
    D_e = dc.diffusion_coefficient_Einstein(times, msd_vals)
    # Use linear fit
    D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals)
    # Use anomalous diffusion fit
    D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals)
    print("Composite for all lipids and both leaflets; whole time range")
    print("Values from estimators:")
    print("  Basic Einstein relation:")
    print("    Diffusion coefficient: {}".format(D_e))
    print("  Linear fit:")
    print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
    print("  Anomalous diffusion fit:")
    print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))
    # diffusion coeff, first 10 ns
    # Simple application of Einstein relation
    D_e = dc.diffusion_coefficient_Einstein(times, msd_vals,
                                            time_range=[t0, t10])
    # Use linear fit
    D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals,
                                              time_range=[t0, t10])
    # Use anomalous diffusion fit
    D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals,
                                                 time_range=[t0, t10])
    print("Composite for all lipids and both leaflets; 0-10 ns")
    print("Values from estimators:")
    print("  Basic Einstein relation:")
    print("    Diffusion coefficient: {}".format(D_e))
    print("  Linear fit:")
    print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
    print("  Anomalous diffusion fit:")
    print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))
    # diffusion coeff,  10-100 ns
    if max(times) > t10:
        # Simple application of Einstein relation
        D_e = dc.diffusion_coefficient_Einstein(times, msd_vals,
                                                time_range=[t10, t100])
        # Use linear fit
        D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals,
                                                  time_range=[t10, t100])
        # Use anomalous diffusion fit
        D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals,
                                                     time_range=[t10, t100])
        print("Composite for all lipids and both leaflets; 10-100 ns")
        print("Values from estimators:")
        print("  Basic Einstein relation:")
        print("    Diffusion coefficient: {}".format(D_e))
        print("  Linear fit:")
        print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
        print("  Anomalous diffusion fit:")
        print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))

    #now each individual lipid type
    for lipid_type in resnames:
        msd_dat = analyzer.get_analysis_data("msd_"+lipid_type)

        times = msd_dat[:, 0]
        msd_vals = msd_dat[:, 1]
        t0 = times[0]
        t10 = t0 + 10000.0
        t100 = t0 + 100000.0
        # diffusion coeff, whole range
        # Simple application of Einstein relation
        D_e = dc.diffusion_coefficient_Einstein(times, msd_vals)
        # Use linear fit
        D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals)
        # Use anomalous diffusion fit
        D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals)
        print("Composite for "+lipid_type+" and both leaflets; whole time range")
        print("Values from estimators:")
        print("  Basic Einstein relation:")
        print("    Diffusion coefficient: {}".format(D_e))
        print("  Linear fit:")
        print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
        print("  Anomalous diffusion fit:")
        print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))
        # diffusion coeff, first 10 ns
        # Simple application of Einstein relation
        D_e = dc.diffusion_coefficient_Einstein(times, msd_vals,
                                                time_range=[t0, t10])
        # Use linear fit
        D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals,
                                                  time_range=[t0, t10])
        # Use anomalous diffusion fit
        D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals,
                                                     time_range=[t0, t10])
        print("Composite for "+lipid_type+" and both leaflets; 0-10 ns")
        print("Values from estimators:")
        print("  Basic Einstein relation:")
        print("    Diffusion coefficient: {}".format(D_e))
        print("  Linear fit:")
        print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
        print("  Anomalous diffusion fit:")
        print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))
        if max(times)>t10:
            # diffusion coeff,  10-100 ns
            # Simple application of Einstein relation
            D_e = dc.diffusion_coefficient_Einstein(times, msd_vals,
                                                    time_range=[t10, t100])
            # Use linear fit
            D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals,
                                                      time_range=[t10, t100])
            # Use anomalous diffusion fit
            D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals,
                                                         time_range=[t10, t100])
            print("Composite for "+lipid_type+" and both leaflets; 10-100 ns")
            print("Values from estimators:")
            print("  Basic Einstein relation:")
            print("    Diffusion coefficient: {}".format(D_e))
            print("  Linear fit:")
            print("    Diffusion coefficient: {} Std Error: {}"
                  .format(D_l[0], D_l[1]))
            print("  Anomalous diffusion fit:")
            print("    Diffusion coefficient: {} Alpha value: {}"
                  .format(D_a[0], D_a[1]))
        # now do individual leaflets
        for leaflet in ['upper', 'lower']:
            msd_dat = analyzer.get_analysis_data("msd_" +
                                                 lipid_type+"_"+leaflet)

            times = msd_dat[:, 0]
            msd_vals = msd_dat[:, 1]
            t0 = times[0]
            t10 = t0 + 10000.0
            t100 = t0 + 100000.0
            # diffusion coeff, whole range
            # Simple application of Einstein relation
            D_e = dc.diffusion_coefficient_Einstein(times, msd_vals)
            # Use linear fit
            D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals)
            # Use anomalous diffusion fit
            D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals)
            print("Composite for " + lipid_type + " and "+leaflet+" leaflet; whole time range")
            print("Values from estimators:")
            print("  Basic Einstein relation:")
            print("    Diffusion coefficient: {}".format(D_e))
            print("  Linear fit:")
            print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
            print("  Anomalous diffusion fit:")
            print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))
            # diffusion coeff, first 10 ns
            # Simple application of Einstein relation
            D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[t0, t10])
            # Use linear fit
            D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[t0, t10])
            # Use anomalous diffusion fit
            D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[t0, t10])
            print("Composite for " + lipid_type + " and "+leaflet+" leaflet; 0-10 ns")
            print("Values from estimators:")
            print("  Basic Einstein relation:")
            print("    Diffusion coefficient: {}".format(D_e))
            print("  Linear fit:")
            print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
            print("  Anomalous diffusion fit:")
            print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))
            if max(times) > t10:
                # diffusion coeff,  10-100 ns
                # Simple application of Einstein relation
                D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[t10, t100])
                # Use linear fit
                D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[t10, t100])
                # Use anomalous diffusion fit
                D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[t10, t100])
                print("Composite for " + lipid_type + " and "+leaflet+" leaflet; 10-100 ns")
                print("Values from estimators:")
                print("  Basic Einstein relation:")
                print("    Diffusion coefficient: {}".format(D_e))
                print("  Linear fit:")
                print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
                print("  Anomalous diffusion fit:")
                print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))


    return


def area_per_lipid(structure_file, trajectory_file, selection_string,
                   frame_start=0, frame_end=-1, frame_interval=1,
                   dump_path=None, n_xbins=100, n_ybins=100, name_dict=None):
    """Protocol to compute area per lipid.

    This function uses the BilayerAnalyzer with the analyses 'apl_box' and
    'apl_grid' to estimate the area per lipid using two estimators: the lateral
    area method (composite) and lipid grid approach. Reports of data are
    printed to standard out, while pickled data is dumped to disk and plots
    are generated and also dumped to disk. The generated files have the prefix
    'apl'.

    Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.
        n_xbins (Optional[int]): Specify the number of bins in the 'x' direction
            to use in the lipid grid. Defaults to 100.
        n_ybins (Optional[int]): Specify the number of bins in the 'y' direction
            to use in the lipid grid. Defaults to 100.

    Returns:
        void

    Notes:
        The lipid grid will have dimensions of n_xbins by n_ybins.
    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    # use a subselection of atoms instead of full lipid center of mass, if given
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    # add the apl analyses
    analyzer.add_analysis("apl_box apl_box")
    analyzer.add_analysis("apl_grid apl_grid")

    analyzer.print_analysis_protocol()

    # add the plots
    analyzer.add_plot("apl apl_box apl_box None")
    analyzer.add_plot("apl apl_p apl_box Box apl_grid None")
    analyzer.add_plot("apl apl_grid apl_grid None")

    # adjust the number of bins for gridding
    analyzer.rep_settings['lipid_grid']['n_xbins'] = n_xbins
    analyzer.rep_settings['lipid_grid']['n_ybins'] = n_ybins

    analyzer.print_plot_protocol()

    # run analysis
    analyzer.run_analysis()

    # output data and plots
    analyzer.dump_data(path=dump_path)
    analyzer.save_all_plots()
    # print final ensemble averages to stdout
    apl_box = analyzer.get_analysis_data('apl_box')

    print("Final running average Area per lipid estimates (squared Angstrom): ")
    print("  via the box dimensions: {:0.4f} +- {:0.4f}".format(apl_box[-1][2], apl_box[-1][3]))
    print("  via the gridding procedure: ")
    apl_grid = analyzer.get_analysis_data('apl_grid')
    for item in apl_grid.keys():
        print("    {}: {:0.4f} +- {:0.4f}".format(item, apl_grid[item][-1][2], apl_grid[item][-1][3]))
    ppb = len(apl_box)/10
    if ppb > 1000:
        ppb=1000
    n_b = 9
    while ppb < 10:
        ppb = len(apl_box)/n_b
        n_b-=1
        if n_b == 0:
            ppb=len(apl_box)
            break
    print("Will do block average with {} points per block".format(ppb))
    block_averager = BlockAverager(points_per_block=ppb)
    block_averager.push_container(apl_box[:,1])
    block_average, std_err = block_averager.get()
    print("Block Averaged ({} points with {} blocks) Area per lipid estimates (squared Angstrom): ".format(len(apl_box),block_averager.number_of_blocks()))
    print("  via the box dimensions: {:0.4f} +- {:0.4f}".format(block_average, std_err))
    print("  via the gridding procedure: ")
    apl_grid = analyzer.get_analysis_data('apl_grid')
    for item in apl_grid.keys():
        block_averager = BlockAverager(points_per_block=ppb)
        block_averager.push_container(apl_grid[item][:, 1])
        block_average, std_err = block_averager.get()
        print("    {}: {:0.4f} +- {:0.4f}".format(item, block_average, std_err))
    return


def bilayer_thickness(structure_file, trajectory_file, selection_string,
                      frame_start=0, frame_end=-1, frame_interval=1,
                      dump_path="./", name_dict=None, n_xbins=100, n_ybins=100,
                      plot_grid_map=True):
    """Protocol to compute the bilayer thickness.

    This function uses the BilayerAnalyzer with the analysis
    'bilayer_thickness' to estimate the the thickness of the bilayer via the
    lipid grid approach. Reports of data are printed to standard out, while
    pickled data (numpy array) is dumped to disk and plots are generated and
    also dumped to disk. The generated files have the prefix 'bt' and the grid
    maps have the prefix 'thickness_grid'.

    Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.
        name_dict (Optional[dict]): A dictionary that defines atoms to use
            for each lipid type when computing the center of mass, which is then
            mapped to the lipid grid. The dictionary should have structure
            {'lipid_resname_1':['atom_a', 'atom_b'], 'lipid_resname_2':['atom_c']}.
            Defaults to None, which means the center of mass of the whole lipid
            is used.
        n_xbins (Optional[int]): Specify the number of bins in the 'x' direction
            to use in the lipid grid. Defaults to 100.
        n_ybins (Optional[int]): Specify the number of bins in the 'y' direction
            to use in the lipid grid. Defaults to 100.
        plot_grid_map (Optional[bool]): Specify whether the 2d grid heat map of
            of the thickness should be generated and output for each analyzed
            frame. Defaults to True.

    Returns:
        void

    Notes:
        The lipid grid will have dimensions of n_xbins by n_ybins.

    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    # use a subselection of atoms instead of full lipid center of mass, if given
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    # add the analysis
    analyzer.add_analysis("bilayer_thickness bt")
    # adjust the number of bins for gridding
    analyzer.rep_settings['lipid_grid']['n_xbins'] = n_xbins
    analyzer.rep_settings['lipid_grid']['n_ybins'] = n_ybins

    analyzer.print_analysis_protocol()
    # add the plot
    analyzer.add_plot("bilayer_thickness bt_p bt None")

    analyzer.print_plot_protocol()

    # run analysis
    for dummy_frame in analyzer:
        if plot_grid_map:
            fs = "frame_{:010d}".format(analyzer.reps['com_frame'].number)
            thickgrid = analyzer.reps['lipid_grid'].thickness_grid()
            xyzc = analyzer.reps['lipid_grid'].get_xyzc(leaflet='lower', color_grid=thickgrid)['lower']

            # with sns.color_palette("PuBuGn_d"):
            pgf.plot_lipid_grid_thickness_map_2d(xyzc[0], xyzc[1], thickgrid,
                                                 filename=dump_path+'thickness_grid_'+fs+'.png',
                                                 vmin=30.0, vmax=50.0, interpolation='gaussian')
            pgf.plot_lipid_grid_thickness_map_2d(xyzc[0], xyzc[1], thickgrid,
                                                 filename=dump_path+'thickness_grid_'+fs+'.pdf',
                                                 vmin=30.0, vmax=50.0, interpolation='gaussian')
        # pgf.plot_grid_as_scatter(xyzc, filename=dump_path + 'thickness_grid_' + fs + '.pdf', colorbar=True, vmin=20.0,
        #                           vmax=50.0)

    # output data and plots
    analyzer.dump_data(path=dump_path)
    analyzer.save_all_plots()
    # output final ensemble average to stdout
    bt = analyzer.get_analysis_data('bt')
    print("Bilayer thickness from gridding procedure (Angstrom): {:0.4f} +- {:0.4f}".format(bt[-1][2], bt[-1][3]))
    bt_avg_sq = (bt[:,1].mean())**2
    bt_sq_avg = (bt[:,1]**2).mean()
    fluctuation = np.sqrt(bt_sq_avg - bt_avg_sq)
    print("Bilayer thickness fluctuation (Angstrom): {:0.4f}".format(fluctuation))
    return


def compressibility(structure_file, trajectory_file, selection_string,
                    frame_start=0, frame_end=-1, frame_interval=1,
                    dump_path=None, temperature=298.15):
    """Protocol to compute the compressibilities.

    This function uses the BilayerAnalyzer with the analyses 'acm', 'vcm', and
    'ac' to estimate the area compressibility modulus, volume
    compressibility modulus, and area compressibility of the system. Reports of
    data are printed to standard out, while pickled data is dumped to disk. The
    generated files have the prefixes 'acm', 'vcm', and 'ac' for area
    compressibility modulus, volume compressibility modulus, and area
    compressibility respectively.

    Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.
        temperature (Optional[float]): Set the temperature (in Kelvin) used
            in the simulation. Defaults to 298.15.

    Returns:
        void

    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')

    # add the analyses
    # area compressibility modulus
    analyzer.add_analysis("acm acm temperature "+str(temperature))
    # volume compressibility modulus
    analyzer.add_analysis("vcm vcm temperature " + str(temperature))
    # area compressibility -- (not modulus)
    analyzer.add_analysis("ac ac temperature " + str(temperature))
    analyzer.print_analysis_protocol()
    analyzer.run_analysis()
    analyzer.dump_data(path=dump_path)
    acm = analyzer.get_analysis_data('acm')
    vcm = analyzer.get_analysis_data('vcm')
    ac = analyzer.get_analysis_data('ac')
    print("Area compressibility modulus: {} mN/m".format(acm[-1:,1]))
    print("Area compressibility: {} m/mN".format(ac[-1:,1]))
    print("Volume compressibility modulus: {} J/Angstrom^3".format(vcm[-1:,1]))

    return


def dispvector_correlation(structure_file, trajectory_file, selection_string,
                           frame_start=0, frame_end=-1, frame_interval=1,
                           dump_path="./", name_dict=None):
    """Protocol to compute lipid displacement vector maps and correlations.

    This function uses the BilayerAnalyzer with the analyses 'disp_vec',
    'disp_vec_corr', and 'disp_vec_nncorr' to compute the lipid displacement
    vectors and their correlations (i.e. cos(theta)). The displacement vectors
    are used to  generate the vector map plots. Reports of data are printed to
    standard out, while pickled data is dumped to disk and plots are generated
    and also dumped to disk. The generated pickle data files have the prefix
    'disp_vec' and the vector map plots have prefix 'step_vector_map'.

    Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.

    Returns:
        void

    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    # use a subselection of atoms instead of full lipid
    # center of mass, if given
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    # add the apl analyses
    # compute the displacment vectors for maps
    analyzer.add_analysis("disp_vec disp_vec_upper scale True wrapped True leaflet upper interval "
                          + str(frame_interval))
    analyzer.add_analysis("disp_vec disp_vec_lower scale True wrapped True leaflet lower interval "
                          + str(frame_interval))
    # compute the correlations between a displacement vector and that lipids
    # closest neighbor in the lateral dimensions
    # analyzer.add_analysis("disp_vec_nncorr disp_vec_nncorr_upper leaflet upper interval "
    #                      + str(frame_interval))
    # analyzer.add_analysis("disp_vec_nncorr disp_vec_nncorr_lower leaflet lower interval "
    #                      + str(frame_interval))
    analyzer.add_analysis("disp_vec_corr_avg disp_vec_corr_avg_upper leaflet upper interval " +
                          str(frame_interval))
    analyzer.add_analysis("disp_vec_corr_avg disp_vec_corr_avg_lower leaflet lower interval " +
                          str(frame_interval))
    analyzer.print_analysis_protocol()

    # run analysis
    analyzer.run_analysis()

    # output data and plots
    analyzer.dump_data(path=dump_path)


    # generate the plots/maps for displacement vectors
    disp_vecs = analyzer.get_analysis_data('disp_vec_upper')
    counter = 0
    number = str(len("{}".format(len(disp_vecs))) + 1)
    form = "{:0"+number+"d}"
    for disp_vec in disp_vecs:
        count = form.format(counter)
        filename = "step_vector_map_upper_"+count+".pdf"
        filename_b = "step_vector_map_upper_"+count+".png"
        pgf.plot_step_vectors(disp_vec, filename=dump_path+filename, scaled=True,
                              wrapped=True)
        pgf.plot_step_vectors(disp_vec, filename=dump_path+filename_b, scaled=True,
                              wrapped=True)
        counter += 1
    disp_vecs = analyzer.get_analysis_data('disp_vec_lower')
    counter = 0
    number = str(len("{}".format(len(disp_vecs))) + 1)
    form = "{:0" + number + "d}"
    for disp_vec in disp_vecs:
        count = form.format(counter)
        filename = "step_vector_map_lower_" + count + ".pdf"
        filename_b = "step_vector_map_lower_" + count + ".png"
        pgf.plot_step_vectors(disp_vec, filename=dump_path+filename, scaled=True,
                              wrapped=True)
        pgf.plot_step_vectors(disp_vec, filename=dump_path+filename_b, scaled=True,
                              wrapped=True)
        counter += 1
    disp_vec_corr_avg_upper = analyzer.get_analysis_data('disp_vec_corr_avg_upper')
    disp_vec_corr_avg_lower = analyzer.get_analysis_data('disp_vec_corr_avg_lower')
    ppb = len(disp_vec_corr_avg_upper)/3
    if ppb > 1000:
        ppb=1000
    n_b = 2
    while ppb < 10:
        ppb = len(disp_vec_corr_avg_upper)/n_b
        n_b-=1
        if n_b == 0:
            ppb=len(disp_vec_corr_avg_upper)
            break
    block_averager_upper = BlockAverager(points_per_block=ppb)
    block_averager_upper.push_container(disp_vec_corr_avg_upper[:,1])
    block_average, std_err = block_averager_upper.get()
    print("Block average of disp_vec_corr_avg in the upper leaflet: {} +- {} using {} blocks and {} points per block.".format(block_average, std_err, block_averager_upper.n_blocks, block_averager_upper.points_per_block()))
    block_averager_lower = BlockAverager(points_per_block=ppb)
    block_averager_lower.push_container(disp_vec_corr_avg_lower[:,1])
    block_average, std_err = block_averager_lower.get()
    print("Block average of disp_vec_corr_avg in the lower leaflet: {} +- {} using {} blocks and {} points per block.".format(block_average, std_err, block_averager_lower.n_blocks, block_averager_lower.points_per_block()))
    return

def PN_orientational_angle(structure_file, trajectory_file, selection_string,
                           lipid_resnames, PN_names='default', frame_start=0,
                           frame_end=-1, frame_interval=1, dump_path=None):
    """Protocol to compute PN orientational angles.

    This function uses the BilayerAnalyzer with the analysis 'loa' to compute
    the PN orientational angles of the specified lipids (by resname). Reports
    of data are printed to standard out, while pickled data (numpy array) is
    dumped to disk. The pickle data files have prefix 'loa'.

    Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        lipid_resnames (list): A list of lipid types as specified by their
            resnames.
        PN_names (Optional(dict)): A dictionary specifying the atom names for
            the phosphorous and nitrogen atoms to use for the PN vector. This
            should have the format {'lipid_resname_1':('phosphorous_name',
            nitrogen_name)}. Defaults to 'default', which uses ('P', 'N') for
            all lipid types.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.

    Returns:
        void

    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    # add the loa analyses
    if PN_names is 'default':
        PN_names = dict()
        for resname in lipid_resnames:
            PN_names[resname] = ('P', 'N')
    for resname in lipid_resnames:
        Pn = PN_names[resname][0]
        Nn = PN_names[resname][1]
        analyzer.add_analysis("loa loa_" + resname +
                              "_upper leaflet upper resname " +
                              resname+" ref_atom_1 " + Pn + " ref_atom_2 " +
                              Nn)
        analyzer.add_analysis("loa loa_" + resname +
                              "_lower leaflet lower resname " +
                              resname + " ref_atom_1 " + Pn + " ref_atom_2 " +
                              Nn)
    # compute the correlations between a displacement vector and that lipids
    # closest neighbor in the lateral dimensions
    analyzer.print_analysis_protocol()

    # run analysis
    analyzer.run_analysis()

    # output data and plots
    analyzer.dump_data(path=dump_path)

    for resname in lipid_resnames:
        loa_upper = analyzer.get_analysis_data("loa_" + resname + "_upper")
        loa_lower = analyzer.get_analysis_data("loa_" + resname + "_lower")
        print("Lipid resname {} has average PN orientation anlge of {} in the upper leaflet".format(resname, loa_upper[-1][1]))
        complement = 90.0 - loa_upper[-1][1]
        print("    complement angle: {}".format(complement))
        print("Lipid resname {} has average PN orientation anlge of {} in the lower leaflet".format(resname, np.abs(loa_lower[-1][1])))
        complement = 90.0 - np.abs(loa_lower[-1][1])
        print("    complement angle: {}".format(complement))
        print(" ")
    return


def nearest_neighbor_fraction(structure_file, trajectory_file,
                              selection_string, lipid_resnames, frame_start=0,
                              frame_end=-1, frame_interval=1, dump_path=None,
                              name_dict=None):
    """Protocol to compute the nearest neigbor fraction.

    This function uses the BilayerAnalyzer with analysis 'nnf' to estimate the
    nearest neighbor fractions (also known as fractional interactions) of the
    specified lipids (by resname). Reports of data are printed to standard out,
    while pickled data is dumped to disk and plots are generated and also
    dumped to disk. The generated files have the prefix 'nnf'.

    Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        lipid_resnames (list): A list of lipid types as specified by their
            resnames.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.

    Returns:
        void

    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    # use a subselection of atoms instead of full lipid
    # center of mass, if given
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    nres = len(lipid_resnames)
    pairs = []
    for i in range(nres):
        for j in range(nres):
            pairs.append([lipid_resnames[i], lipid_resnames[j]])
    # add the loa analyses
    for pair in pairs:
        l1 = pair[0]
        l2 = pair[1]
        analyzer.add_analysis("nnf nnf_" + l1 + "_" + l2 + " resname_1 " +
                              l1 + " resname_2 " + l2)
    # compute the correlations between a displacement vector and that lipids
    # closest neighbor in the lateral dimensions
    analyzer.print_analysis_protocol()
    print(" ")

    # run analysis
    analyzer.run_analysis()

    # output data and plots
    analyzer.dump_data(path=dump_path)

    for pair in pairs:
        l1 = pair[0]
        l2 = pair[1]
        t_nnf = analyzer.get_analysis_data("nnf_"+l1+"_"+l2)
        print("Nearest neighbor fraction for lipid pair {} and {} : {:0.4f} +- {:0.4f}".format(l1, l2, t_nnf[-1][2], t_nnf[-1][3]))
    return


def normal_displacement_lipid_type_correlation(structure_file, trajectory_file,
                                               selection_string, frame_start=0,
                                               frame_end=-1, frame_interval=1,
                                               dump_path=None,
                                               com_sub_selection_dict=None):
    """Protocol for correlation between the normal dimension displacement and lipid type.

    This function uses the BilayerAnalyzer with the analysis 'ndcorr' to do a
    local grid analysis to estimate the correlation between between lipid types
    and their deflection along the bilayer normal; this analysis includes all
    lipid types. Reports of data are printed to standard out, while pickled
    data is dumped to disk and plots are generated and also dumped to disk. The
    generated pickle files have the prefix 'norm_disp_correlation' while the
    plots have the prefix 'ndcorr'.

    Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.

    Returns:
        void

    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                             trajectory=trajectory_file,
                             selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')

    if com_sub_selection_dict is not None:
        analyzer.rep_settings['com_frame']['name_dict'] = \
                                            com_sub_selection_dict
    analyzer.add_analysis("ndcorr norm_disp_correlation")
    # compute the correlations between a displacement vector and that lipids
    # closest neighbor in the lateral dimensions
    analyzer.print_analysis_protocol()
    print("--------")

    # run analysis
    analyzer.run_analysis()

    # output data and plots
    analyzer.dump_data(path=dump_path)

    ndcorr = analyzer.get_analysis_data('norm_disp_correlation')
    print(" ")
    print("Normal dimension displacement-lipid type cross correlation results:"
          )
    for leaflet in ndcorr.keys():
        print("  {} leaflet:".format(leaflet))
        for lipid_resname in ndcorr[leaflet].keys():
            mean = ndcorr[leaflet][lipid_resname][-1][2]
            deviation = ndcorr[leaflet][lipid_resname][-1][3]
            print("    Lipid resname {}: {:0.4f} +- {:0.4f}".format(lipid_resname, mean, deviation))

    pgf.plot_displacement_lipid_type_cross_correlation(ndcorr,
                                                       filename="ndcorr.png")
    pgf.plot_displacement_lipid_type_cross_correlation(ndcorr,
                                                       filename="ndcorr.pdf")

    return


def lipid_grid_maps(structure_file, trajectory_file, selection_string,
                    frame_start=0, frame_end=-1, frame_interval=1,
                    dump_path="./", name_dict=None, n_xbins=100, n_ybins=100,
                    type_colors='auto'):
    """Protocol to generate lipid grid map plots.

    This function uses the BilayerAnalyzer with the 'lipid_grid' reprentation
    to estimate the the thickness of the bilayer via the lipid grid approach.
    Grid map lots are generated and dumped to disk. The generated files have
    the prefix 'lipid_grid'.

    Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.
        name_dict (Optional[dict]): A dictionary that defines atoms to use
            for each lipid type when computing the center of mass, which is
            then mapped to the lipid grid. The dictionary should have structure
            {'lipid_resname_1':['atom_a', 'atom_b'],
            'lipid_resname_2':['atom_c']}. Defaults to None, which means the
            center of mass of the whole lipid is used.
        n_xbins (Optional[int]): Specify the number of bins in the 'x'
            direction to use in the lipid grid. Defaults to 100.
        n_ybins (Optional[int]): Specify the number of bins in the 'y'
            direction to use in the lipid grid. Defaults to 100.
        type_colors (Optional[dict]): Specify the colors to use for for each
            lipid type when generating the 2d grid map plot. The dictionary
            should be keyed by the lipid types' resnames. Defaults to 'auto',
            which automatically assigns colors for each lipid type.

    Returns:
        void

    Notes:
        The lipid grid will have dimensions of n_xbins by n_ybins.

    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    # use a subselection of atoms instead of full lipid
    # center of mass, if given
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    # analyzer.settings['print_interval']=1
    # add the analysis
    if type_colors is 'auto':
        type_colors = {}
        for dummy_frame in analyzer:
            unique_resnames = sorted(analyzer.reps['com_frame'].unique_resnames())
            i = 0
            for resname in unique_resnames:
                type_colors[resname] = _color_list[i]
                i+=1
                if i == len(_color_list):
                    i=0
            break
        analyzer.reset()

    analyzer._analysis_protocol.use_objects['lipid_grid'] = True
    # adjust the number of bins for gridding
    analyzer.rep_settings['lipid_grid']['n_xbins'] = n_xbins
    analyzer.rep_settings['lipid_grid']['n_ybins'] = n_ybins
    #run analysis
    for dummy_frame in analyzer:
        fs = "frame_{:010d}".format(analyzer.reps['com_frame'].number)
        xyzc = analyzer.reps['lipid_grid'].get_xyzc(leaflet='lower', color_type_dict=type_colors)['lower']

        pgf.plot_grid_as_scatter(xyzc, filename=dump_path + 'lipid_grid_lower_' + fs + '.png')

        pgf.plot_grid_as_scatter(xyzc, filename=dump_path + 'lipid_grid_lower_' + fs + '.pdf')

        xyzc = analyzer.reps['lipid_grid'].get_xyzc(leaflet='upper', color_type_dict=type_colors)['upper']

        pgf.plot_grid_as_scatter(xyzc, filename=dump_path + 'lipid_grid_upper_' + fs + '.png')

        pgf.plot_grid_as_scatter(xyzc, filename=dump_path + 'lipid_grid_upper_' + fs + '.pdf')


    for resname in sorted(type_colors.keys()):
        print("lipid resname {} is color {}".format(resname,type_colors[resname]))
    return

def distance_cutoff_clustering(structure_file, trajectory_file,
                               selection_string, resnames, frame_start=0,
                               frame_end=-1, frame_interval=1, dump_path="./",
                               name_dict = None):
    """Protocol to compute hiearchical distance cutoff clusters.

    This function uses the BilayerAnalyzer to determine the hiearchical
    distance cutoff clusters of lipids in the lateral dimensions; clustering is
    performed between lipids of the same type using a distance cutoff of
    12.0*distance_units (assumed to be Angstroms). The protocol extracts the
    number of clusters and average cluster size for each lipid type in both the
    'upper' and 'lower' leaflets. Reports of data are printed to standard out,
    while pickled data is dumped to disk and plots are generated and also
    dumped to disk. The generated files have the prefix 'dc_cluster'.

    Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.

    Returns:
        void

    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    # add the analyses
    for resname in resnames:
        add_in = "dc_cluster dc_cluster_{}_upper resname {} leaflet upper".format(resname, resname)
        analyzer.add_analysis(add_in)
        add_in = "dc_cluster dc_cluster_{}_lower resname {} leaflet lower".format(resname, resname)
        analyzer.add_analysis(add_in)
    # use a subselection of atoms instead of full lipid
    # center of mass, if given
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    analyzer.dump_data(path=dump_path)
    print("         ")
    for resname in resnames:
        results = analyzer.get_analysis_data("dc_cluster_{}_upper".format(resname))
        plotname_png = "{}dc_cluster_upper_{}_nclusters.png".format(dump_path, resname)
        plotname_eps = "{}dc_cluster_upper_{}_nclusters.pdf".format(dump_path, resname)
        times = results['nclusters'][:,0]
        means = results['nclusters'][:,2]
        stds = results['nclusters'][:,3]
        pgf.plot_dc_cluster_dat_number([(times, means, stds)], filename=plotname_png)
        pgf.plot_dc_cluster_dat_number([(times, means, stds)], filename=plotname_eps)
        plotname_png = "{}dc_cluster_upper_{}_avg_size.png".format(dump_path, resname)
        plotname_eps = "{}dc_cluster_upper_{}_avg_size.pdf".format(dump_path, resname)
        means = results['avg_size'][:, 2]
        stds = results['avg_size'][:, 3]
        pgf.plot_dc_cluster_dat_number([(times, means, stds)], filename=plotname_png)
        pgf.plot_dc_cluster_dat_number([(times, means, stds)], filename=plotname_eps)
        print("resname {} has distance cutoff cluster results in upper leaflet: ".format(resname))
        for key in results.keys():
            if key is not 'clusters':
                print("  mean {}: {} +- {}".format(key, results[key][-1][2], results[key][-1][3]))

        results = analyzer.get_analysis_data("dc_cluster_{}_lower".format(resname))
        plotname_png = "{}dc_cluster_lower_{}_nclusters.png".format(dump_path,resname)
        plotname_eps = "{}dc_cluster_lower_{}_nclusters.pdf".format(dump_path, resname)
        times = results['nclusters'][:,0]
        means = results['nclusters'][:,2]
        stds = results['nclusters'][:,3]
        pgf.plot_dc_cluster_dat_number([(times, means, stds)], filename=plotname_png)
        pgf.plot_dc_cluster_dat_number([(times, means, stds)], filename=plotname_eps)
        plotname_png = "{}dc_cluster_lower_{}_avg_size.png".format(dump_path, resname)
        plotname_eps = "{}dc_cluster_lower_{}_avg_size.pdf".format(dump_path, resname)
        means = results['avg_size'][:, 2]
        stds = results['avg_size'][:, 3]
        pgf.plot_dc_cluster_dat_number([(times, means, stds)], filename=plotname_png)
        pgf.plot_dc_cluster_dat_number([(times, means, stds)], filename=plotname_eps)
        print("resname {} has distance cutoff cluster results in lower leaflet: ".format(resname))
        for key in results.keys():
            if key is not 'clusters':
                print("  mean {}: {} +- {}".format(key, results[key][-1][2], results[key][-1][3]))
        print("------------------------------------------")

    return


def position_density_maps_2d(structure_file, trajectory_file,
                             bilayer_selection_string, resnames,
                             frame_start=0, frame_end=-1, frame_interval=1,
                             dump_path="./"):
    """Protocol to compute the 2d positional densities and its plots.

    This function uses the
    mda_tools.mda_density_map.position_density_map_2d_leaflet_simple function
    to compute the 2d positional densities of lipids in the lateral dimensions
    and generates their plots. Separate plots are generated for each lipid
    type, denoted by their resnames. The plots are generated and dumped to
    disk. The output files have the prefix 'position_density_2d'.

     Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        resnames (list): Specify the resnames of the lipid types
            that are to be included in this analysis.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.

    Returns:
        void

    """
    u = mda.Universe(structure_file, trajectory_file)
    bilayer_sel = u.select_atoms(bilayer_selection_string)
    x_centers, y_centers, counts \
        = position_density_map_2d_leaflet_simple(u, bilayer_sel, resnames,
                                                 fstart=frame_start,
                                                 fend=frame_end,
                                                 fstep=frame_interval,
                                                 refsel=bilayer_sel,
                                                 scale_to_max=False)
    for key in counts.keys():
        outname_eps = dump_path+'position_density_2d_{}_upper.pdf'.format(key)
        outname_png = dump_path+'position_density_2d_{}_upper.png'.format(key)
        pgf.plot_position_density_map_2d(x_centers,
                                         y_centers,
                                         counts[key]['upper'],
                                         filename=outname_eps,
                                         scaled_to_max=False,
                                         interpolation='gaussian')
        pgf.plot_position_density_map_2d(x_centers,
                                         y_centers,
                                         counts[key]['upper'],
                                         filename=outname_png,
                                         scaled_to_max=False,
                                         interpolation='gaussian')
        outname_eps = dump_path+'position_density_2d_{}_lower.pdf'.format(key)
        outname_png = dump_path+'position_density_2d_{}_lower.png'.format(key)
        pgf.plot_position_density_map_2d(x_centers,
                                         y_centers,
                                         counts[key]['lower'],
                                         filename=outname_eps,
                                         scaled_to_max=False,
                                         interpolation='gaussian')
        pgf.plot_position_density_map_2d(x_centers,
                                         y_centers,
                                         counts[key]['lower'],
                                         filename=outname_png,
                                         scaled_to_max=False,
                                         interpolation='gaussian')
    return

def curvature_grid(structure_file, trajectory_file, selection_string,
                   frame_start=0, frame_end=-1, frame_interval=1,
                   dump_path="./", name_dict=None, n_xbins=100, n_ybins=100,
                   curvature_grid_vmin=None, curvature_grid_vmax=None):
    """Protocol to generate curvature estimates over a grid.

    This function uses the BilayerAnalyzer with the 'lipid_grid' reprentation
    to estimate the curvature of the bilayer via the lipid grid approach.
    Data and plots are generated and dumped to disk. The generated files have
    the prefix 'curvature_grid'.

    Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.
        name_dict (Optional[dict]): A dictionary that defines atoms to use
            for each lipid type when computing the center of mass, which is
            then mapped to the lipid grid. The dictionary should have structure
            {'lipid_resname_1':['atom_a', 'atom_b'],
            'lipid_resname_2':['atom_c']}. Defaults to None, which means the
            center of mass of the whole lipid is used.
        n_xbins (Optional[int]): Specify the number of bins in the 'x'
            direction to use in the lipid grid. Defaults to 100.
        n_ybins (Optional[int]): Specify the number of bins in the 'y'
            direction to use in the lipid grid. Defaults to 100.

    Returns:
        void

    Notes:
        The lipid grid will have dimensions of n_xbins by n_ybins.

    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    # use a subselection of atoms instead of full lipid
    # center of mass, if given
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    # analyzer.settings['print_interval']=1
    # add the analysis
    analyzer._analysis_protocol.use_objects['lipid_grid'] = True
    # adjust the number of bins for gridding
    analyzer.rep_settings['lipid_grid']['n_xbins'] = n_xbins
    analyzer.rep_settings['lipid_grid']['n_ybins'] = n_ybins
    upper_data = []
    lower_data = []
    times = []
    #run analysis
    for dummy_frame in analyzer:
        fs = "frame_{:010d}".format(analyzer.reps['com_frame'].number)
        curvatures = analyzer.reps['lipid_grid'].curvature(filter_sigma=5.0)
        #upper leaflet
        upper_mean = curvatures[0][0]
        x_centers = analyzer.reps['lipid_grid'].leaf_grid['upper'].x_centers
        y_centers = analyzer.reps['lipid_grid'].leaf_grid['upper'].y_centers
        out_tuple = (x_centers, y_centers, upper_mean)
        with open(dump_path+"curvature_grid_xy_mean_curvature_grid_upper_"+fs+".pickle", 'wb') as outfile:
            pickle.dump(out_tuple, outfile)
        pgf.plot_xygrid_as_imshow(x_centers, y_centers, upper_mean,
                                  filename=dump_path+"curvature_grid_upper_"+fs+".pdf",
                                  xlabel='x ($\AA$)', ylabel='y ($\AA$)',
                                  colorbar=True, colorbarlabel='Mean curvature ($\AA^{-1}$)',
                                  vmin=curvature_grid_vmin, vmax=curvature_grid_vmax)
        pgf.plot_xygrid_as_imshow(x_centers, y_centers, upper_mean,
                                  filename=dump_path+"curvature_grid_upper_"+fs+".png",
                                  xlabel='x ($\AA$)', ylabel='y ($\AA$)',
                                  colorbar=True, colorbarlabel='Mean curvature ($\AA^{-1}$)',
                                  vmin=curvature_grid_vmin, vmax=curvature_grid_vmax)
        avg_mean = upper_mean.mean()
        max_mean = upper_mean.max()
        min_mean = upper_mean.min()
        upper_data.append([avg_mean, max_mean, min_mean])
        #lower leaflet
        lower_mean = curvatures[1][0]
        x_centers = analyzer.reps['lipid_grid'].leaf_grid['lower'].x_centers
        y_centers = analyzer.reps['lipid_grid'].leaf_grid['lower'].y_centers
        out_tuple = (x_centers, y_centers, lower_mean)
        with open(dump_path+"curvature_grid_xy_mean_curvature_grid_lower_"+fs+".pickle", 'wb') as outfile:
            pickle.dump(out_tuple, outfile)
        pgf.plot_xygrid_as_imshow(x_centers, y_centers, lower_mean,
                                  filename=dump_path+"curvature_grid_lower_"+fs+".pdf",
                                  xlabel='x ($\AA$)', ylabel='y ($\AA$)',
                                  colorbar=True, colorbarlabel='Mean curvature ($\AA^{-1}$)',
                                  vmin=curvature_grid_vmin, vmax=curvature_grid_vmax)
        pgf.plot_xygrid_as_imshow(x_centers, y_centers, lower_mean,
                                  filename=dump_path+"curvature_grid_lower_"+fs+".png",
                                  xlabel='x ($\AA$)', ylabel='y ($\AA$)',
                                  colorbar=True, colorbarlabel='Mean curvature ($\AA^{-1}$)',
                                  vmin=curvature_grid_vmin, vmax=curvature_grid_vmax)
        avg_mean = lower_mean.mean()
        max_mean = lower_mean.max()
        min_mean = lower_mean.min()
        lower_data.append([avg_mean, max_mean, min_mean])
        times.append(analyzer.reps['com_frame'].time)
        um_max = np.abs(upper_mean).max()
        um_factor = 1.0
        if um_max < 0.10:
            um_ratio = 0.10/um_max
            um_factor = um_ratio*100.0
        lm_max = np.abs(lower_mean).max()
        lm_factor = 1.0
        if lm_max < 0.10:
            lm_ratio = 0.10/lm_max
            lm_factor = 100.0*lm_ratio
        analyzer.reps['lipid_grid'].write_pdb(dump_path+"curvature_grid_lipid_grid_"+fs+".pdb",
                                              beta_grid_upper=upper_mean*um_factor,
                                              beta_grid_lower=lower_mean*lm_factor,
                                              use_gaussian_filter=True,
                                              filter_sigma=5.0)
    upper_data = np.array(upper_data)
    lower_data = np.array(lower_data)
    times = np.array(times)/1000.0
    with open(dump_path+"curvature_grid_times_upper_data.pickle", 'wb') as outfile:
        pickle.dump((times, upper_data), outfile)
    with open(dump_path+"curvature_grid_times_lower_data.pickle", 'wb') as outfile:
        pickle.dump((times, lower_data), outfile)
    pgf.plot([(times, upper_data[:, 0]), (times, lower_data[:, 0])],
             name_list=['upper leaflet', 'lower leaflet'], show=False,
             save=True, xlabel='Time (ns)', ylabel='Average Mean Curvature ($\AA^{-1}$)',
             filename=dump_path+"curvature_grid_average_mean_curvature.pdf")
    pgf.plot([(times, upper_data[:, 0]), (times, lower_data[:, 0])],
             name_list=['upper leaflet', 'lower leaflet'], show=False,
             save=True, xlabel='Time (ns)', ylabel='Average Mean Curvature ($\AA^{-1}$)',
             filename=dump_path+"curvature_grid_average_mean_curvature.png")
    print("Mean Average Mean Curvature in upper leaflet: {}".format(upper_data[:,0].mean()))
    print("Mean Average Mean Curvature in lower leaflet: {}".format(lower_data[:,0].mean()))
    pgf.plot([(times, upper_data[:, 1]), (times, lower_data[:, 1])],
             name_list=['upper leaflet', 'lower leaflet'], show=False,
             save=True, xlabel='Time (ns)', ylabel='Max Mean Curvature ($\AA^{-1}$)',
             filename=dump_path+"curvature_grid_max_mean_curvature.pdf")
    pgf.plot([(times, upper_data[:, 1]), (times, lower_data[:, 1])],
             name_list=['upper leaflet', 'lower leaflet'], show=False,
             save=True, xlabel='Time (ns)', ylabel='Max Mean Curvature ($\AA^{-1}$)',
             filename=dump_path+"curvature_grid_max_mean_curvature.png")
    print("Max Max Mean Curvature in upper leaflet: {}".format(upper_data[:, 1].max()))
    print("Max Max Mean Curvature in lower leaflet: {}".format(lower_data[:, 1].max()))
    pgf.plot([(times, upper_data[:, 2]), (times, lower_data[:, 2])],
             name_list=['upper leaflet', 'lower leaflet'], show=False,
             save=True, xlabel='Time (ns)', ylabel='Min Mean Curvature ($\AA^{-1}$)',
             filename=dump_path+"curvature_grid_min_mean_curvature.pdf")
    pgf.plot([(times, upper_data[:, 2]), (times, lower_data[:, 2])],
             name_list=['upper leaflet', 'lower leaflet'], show=False,
             save=True, xlabel='Time (ns)', ylabel='Min Mean Curvature ($\AA^{-1}$)',
             filename=dump_path+"curvature_grid_min_mean_curvature.png")
    pgf.plot([(times, upper_data[:, 0]), (times, lower_data[:, 0]),
             (times, upper_data[:, 1]), (times, lower_data[:, 1]),
             (times, upper_data[:, 2]), (times, lower_data[:, 2])],
             name_list=['Avg-ul', 'Avg-ll',
             'Max-ul', 'Max-ll', 'Min-ul',
             'Min-ll'], show=False,
             save=True, xlabel='Time (ns)', ylabel='Mean Curvature ($\AA^{-1}$)',
             filename=dump_path+"curvature_grid_all_mean_curvature.pdf")
    pgf.plot([(times, upper_data[:, 0]), (times, lower_data[:, 0]),
             (times, upper_data[:, 1]), (times, lower_data[:, 1]),
             (times, upper_data[:, 2]), (times, lower_data[:, 2])],
             name_list=['Avg-ul', 'Avg-ll',
             'Max-ul', 'Max-ll', 'Min-ul',
             'Min-ll'], show=False,
             save=True, xlabel='Time (ns)', ylabel='Mean Curvature ($\AA^{-1}$)',
             filename=dump_path+"curvature_grid_all_mean_curvature.png")
    print("Min Min Mean Curvature in the upper leaflet: {}".format(upper_data[:, 2].min()))
    print("Min Min Mean Curvature in the lower leaflet: {}".format(lower_data[:, 2].min()))
    return

def com_lateral_rdf(structure_file, trajectory_file,
                    bilayer_selection_string, resnames, name_dict=None,
                    frame_start=0, frame_end=-1, frame_interval=1,
                    dump_path="./"):
    """Protocol to compute the 2d RDFs for lipid types in the bilayer lateral plane.

    This function uses the
    BilayerAnalyzer with the
    bilayer_anlayzer.analysis_protocols.COMLateralRDFProtocol analysis protocol
    to compute the 2d RDFs of lipids in the lateral dimensions of the bilayer
    and generates their plots. Separate plots are generated for each pair of lipid
    types (denoted by their resnames). The plots are generated and dumped to
    disk. The output files have the prefix 'com_lateral_rdf'.

     Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        resnames (list): Specify the resnames of the lipid types
            that are to be included in this analysis.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.

    Returns:
        void

    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=bilayer_selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    # use a subselection of atoms instead of full lipid
    # center of mass, if given
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    res_pairs = []
    for lipid_type_a in resnames:
        for lipid_type_b in resnames:
            if [lipid_type_a, lipid_type_b] not in res_pairs:
                res_pairs.append([lipid_type_a, lipid_type_b])
    # add the analyses
    for pair in res_pairs:
        add_in = "com_lateral_rdf com_lateral_rdf_{}-{}_upper resname_1 {} resname_2 {} leaflet upper range_outer 40.0 n_bins 80".format(pair[0], pair[1], pair[0], pair[1])
        analyzer.add_analysis(add_in)
        add_in = "com_lateral_rdf com_lateral_rdf_{}-{}_lower resname_1 {} resname_2 {} leaflet lower range_outer 40.0 n_bins 80".format(pair[0], pair[1], pair[0], pair[1])
        analyzer.add_analysis(add_in)

    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    analyzer.dump_data(path=dump_path)

    # Generate plots
    all_data = dict()
    all_data['upper'] = []
    all_data['lower'] = []
    all_data['both'] = []
    all_names = dict()
    all_names['upper'] = []
    all_names['lower'] = []
    all_names['both'] = []
    hl_x = np.array([0.0, 40.0])
    hl_y = np.array([1.0, 1.0])
    hl = (hl_x, hl_y)
    for key in all_names.keys():
        all_names[key].append("1.0")
        all_data[key].append(hl)
    for pair in res_pairs:
        item = "com_lateral_rdf_{}-{}_upper".format(pair[0], pair[1])
        rdf, bins = analyzer.get_analysis_data(item)
        rdf_u = rdf
        name = "{}-{}".format(pair[0], pair[1])
        all_data['upper'].append((bins, rdf))
        all_names['upper'].append(name)
        pgf.plot([(bins, rdf)], filename=dump_path+item+".png", xlabel="Radial Distance ($\AA$)", ylabel="Radial Distribution")
        pgf.plot([(bins, rdf)], filename=dump_path+item+".eps", xlabel="Radial Distance ($\AA$)", ylabel="Radial Distribution")

        item = "com_lateral_rdf_{}-{}_lower".format(pair[0], pair[1])
        rdf, bins = analyzer.get_analysis_data(item)
        rdf_l = rdf
        name = "{}-{}".format(pair[0], pair[1])
        all_data['lower'].append((bins, rdf))
        all_names['lower'].append(name)

        pgf.plot([hl, (bins, rdf)], filename=dump_path+item+".png", xlabel="Radial Distance ($\AA$)", ylabel="Radial Distribution")
        pgf.plot([hl, (bins, rdf)], filename=dump_path+item+".eps", xlabel="Radial Distance ($\AA$)", ylabel="Radial Distribution")

        item = "com_lateral_rdf_{}-{}_both".format(pair[0], pair[1])
        rdf = (rdf_u + rdf_l)/2.0
        name = "{}-{}".format(pair[0], pair[1])
        all_data['both'].append((bins, rdf))
        all_names['both'].append(name)
        pgf.plot([hl, (bins, rdf)], filename=dump_path+item+".png", xlabel="Radial Distance ($\AA$)", ylabel="Radial Distribution")
        pgf.plot([hl, (bins, rdf)], filename=dump_path+item+".eps", xlabel="Radial Distance ($\AA$)", ylabel="Radial Distribution")
    for key in all_data.keys():
        pgf.plot(all_data[key], name_list=all_names[key], filename=dump_path+"com_lateral_rdf_all_"+key+".png", xlabel="Radial Distance ($\AA$)", ylabel="Radial Distribution")
        pgf.plot(all_data[key], name_list=all_names[key], filename=dump_path+"com_lateral_rdf_all_"+key+".eps", xlabel="Radial Distance ($\AA$)", ylabel="Radial Distribution")
    return

def _estimate_correlation_length(bins, averages):
    """Estimate the correlation length of a spatial velocity correlation
    function. The correlation length is typically defined as the radial
    distance at which the correlation first becomes zero.
    """
    n_point = len(bins)
    for i in range(0, n_point-1):
        bottom = averages[i]
        top = averages[i+1]
        if bottom > 0.0 and top < 0.0:
            return (bins[i]+bins[i+1])/2.0
    return None

def spatial_velocity_correlation_functions(structure_file, trajectory_file,
                    bilayer_selection_string, resnames, name_dict=None,
                    frame_start=0, frame_end=-1, frame_interval=1,
                    dump_path="./"):
    """Protocol to compute the 2d spatial velocity correlation functions for
    lipid types in the bilayer lateral plane.

    This function uses the
    BilayerAnalyzer with the
    bilayer_anlayzer.analysis_protocols.SpatialVelocityCorrelationFunctionProtocol
    analysis protocol
    to compute the 2d spatial velocity correlation functions lipid-lipid interactions
    in the lateral dimensions of the bilayer
    and generates their plots. Separate plots are generated for each pair of lipid
    types (denoted by their resnames) in each leaflet, as well as the composite
    of all lipids and over both leaflets. The plots are generated and dumped to
    disk. The output files have the prefix 'spatial_velocity_corr'.

    Args:
        structure_file (str): The path and filename of the structure file.
        trajectory_file (str, list): The path and filename of the trajectory
            file. Also accepts a list of path/filenames.
        selection_string (str): The MDAnalysis compatible string used to select
            the bilayer components.
        resnames (list): Specify the resnames of the lipid types
            that are to be included in this analysis.
        frame_start (Optional[int]): Specify the index of the first frame
            of the trajectory to include in the analysis. Defaults to 0,
            which is the first frame of the trajectory.
        frame_end (Optional[int]): Specify the index of the last frame
            of the trajectory to include in the analysis. Defaults to -1,
            which is the last frame of the trajectory.
        frame_interval (Optional[int]): Specify the interval between frames of
            the trajectory to include in the analysis, or the frame frequency.
            Defaults to 1, which includes all frames [frame_start, frame_end].
        dump_path (Optional[str]): Specify a file path for the output files.
            Defaults to None, which outputs in the current directory.

    Returns:
        void

    """
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=bilayer_selection_string)

    analyzer.set_frame_range(frame_start, frame_end, 1)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    # use a subselection of atoms instead of full lipid
    # center of mass, if given
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    res_pairs = []
    for lipid_type_a in resnames:
        #res_pairs.append([lipid_type_a, 'all'])
        for lipid_type_b in resnames:
            if [lipid_type_a, lipid_type_b] not in res_pairs:
                res_pairs.append([lipid_type_a, lipid_type_b])
    # add the analyses
    res_pairs.append(['all', 'all'])
    for pair in res_pairs:
        add_in = "spatial_velocity_corr spatial_velocity_corr_{}-{}_upper resname_1 {} resname_2 {} leaflet upper range_outer 75.0 n_bins 150 interval {}".format(pair[0], pair[1], pair[0], pair[1], frame_interval)
        analyzer.add_analysis(add_in)
        add_in = "spatial_velocity_corr spatial_velocity_corr_{}-{}_lower resname_1 {} resname_2 {} leaflet lower range_outer 75.0 n_bins 150 interval {}".format(pair[0], pair[1], pair[0], pair[1], frame_interval)
        analyzer.add_analysis(add_in)
        add_in = "spatial_velocity_corr spatial_velocity_corr_{}-{}_both resname_1 {} resname_2 {} leaflet both range_outer 75.0 n_bins 150 interval {}".format(pair[0], pair[1], pair[0], pair[1], frame_interval)
        analyzer.add_analysis(add_in)

    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    analyzer.dump_data(path=dump_path)

    # Generate plots
    all_data = dict()
    all_data['upper'] = []
    all_data['lower'] = []
    all_data['both'] = []
    all_names = dict()
    all_names['upper'] = []
    all_names['lower'] = []
    all_names['both'] = []
    for pair in res_pairs:
        item = "spatial_velocity_corr_{}-{}_upper".format(pair[0], pair[1])
        bins, averages = analyzer.get_analysis_data(item)
        name = "{}-{}".format(pair[0], pair[1])
        all_data['upper'].append((bins, averages))
        all_names['upper'].append(name)
        pgf.plot([(bins, averages)], filename=dump_path+item+".png", xlabel="Radial Distance ($\AA$)", ylabel="Velocity Correlation")
        pgf.plot([(bins, averages)], filename=dump_path+item+".eps", xlabel="Radial Distance ($\AA$)", ylabel="Velocity Correlation")
        corr_dist = _estimate_correlation_length(bins, averages)
        print("Correlation distance for pair {}-{} in the upper leaflet: {} Angstroms".format(pair[0], pair[1], corr_dist))
        item = "spatial_velocity_corr_{}-{}_lower".format(pair[0], pair[1])
        bins, averages = analyzer.get_analysis_data(item)
        name = "{}-{}".format(pair[0], pair[1])
        all_data['lower'].append((bins, averages))
        all_names['lower'].append(name)
        pgf.plot([(bins, averages)], filename=dump_path+item+".png", xlabel="Radial Distance ($\AA$)", ylabel="Velocity Correlation")
        pgf.plot([(bins, averages)], filename=dump_path+item+".eps", xlabel="Radial Distance ($\AA$)", ylabel="Velocity Correlation")
        corr_dist = _estimate_correlation_length(bins, averages)
        print("Correlation distance for pair {}-{} in the lower leaflet: {} Angstroms".format(pair[0], pair[1], corr_dist))

        item = "spatial_velocity_corr_{}-{}_both".format(pair[0], pair[1])
        bins, averages = analyzer.get_analysis_data(item)
        name = "{}-{}".format(pair[0], pair[1])
        all_data['both'].append((bins, averages))
        all_names['both'].append(name)
        pgf.plot([(bins, averages)], filename=dump_path+item+".png", xlabel="Radial Distance ($\AA$)", ylabel="Velocity Correlation")
        pgf.plot([(bins, averages)], filename=dump_path+item+".eps", xlabel="Radial Distance ($\AA$)", ylabel="Velocity Correlation")
        corr_dist = _estimate_correlation_length(bins, averages)
        print("Correlation distance for pair {}-{} in both leaflets: {} Angstroms".format(pair[0], pair[1], corr_dist))
    for key in all_data.keys():
        pgf.plot(all_data[key], name_list=all_names[key], filename=dump_path+"spatial_velocity_corr_all_"+key+".png", xlabel="Radial Distance ($\AA$)", ylabel="Velocity Correlation")
        pgf.plot(all_data[key], name_list=all_names[key], filename=dump_path+"spatial_velocity_corr_all_"+key+".eps", xlabel="Radial Distance ($\AA$)", ylabel="Velocity Correlation")
    return
