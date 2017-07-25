import numpy as np

from pybilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer
from pybilt.mda_tools import diffusion_coefficients as dc




def msd_diffusion(structure_file, trajectory_file, selection_string, resnames=[], frame_start=0, frame_end=-1, frame_interval=1, dump_path=None):

    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    #for each lipid resname in the input resnames add an msd analysis for both, upper, and lower leaflets
    for lipid_type in resnames:
        a_name = "msd_"+lipid_type
        add_string = "msd "+a_name+" resname "+lipid_type
        add_string_upper = "msd "+a_name+"_upper resname "+lipid_type+" leaflet upper"
        add_string_lower = "msd "+a_name+"_lower resname "+lipid_type+" leaflet lower"
        analyzer.add_analysis(add_string)
        analyzer.add_analysis(add_string_upper)
        analyzer.add_analysis(add_string_lower)

    #now add plot protocols
    analyzer.add_plot('msd msd_a msd_1 All')
    # plot with each of the lipid types on same plot (leaflet=both)
    l_msd_string = "msd msd_l"
    for lipid_type in resnames:
        l_msd_string+=" msd_"+lipid_type+" "+lipid_type
    analyzer.add_plot(l_msd_string)

    #plots for each lipid with upper, lower, both on same plot
    for lipid_type in resnames:
        add_string = "msd "+lipid_type+"_cul msd_"+lipid_type+" Composite msd_"+lipid_type+"_upper Upper msd_"+lipid_type+"_lower Lower"
        analyzer.add_plot(add_string)

    analyzer.print_analysis_protocol()
    analyzer.print_plot_protocol()

    #run the analyzer
    analyzer.run_analysis()

    #dump the output
    analyzer.dump_data(path=dump_path)

    #make MSD vs time plots
    analyzer.save_all_plots()

    #compute diffusion coefficients
    #composite

    msd_dat = analyzer.get_analysis_data('msd_1')

    times = msd_dat[:,0]
    msd_vals = msd_dat[:,1]
    #diffusion coeff, whole range
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
    #diffusion coeff, first 10 ns
    # Simple application of Einstein relation
    D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[0.0, 10000.0])
    # Use linear fit
    D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[0.0, 10000.0])
    # Use anomalous diffusion fit
    D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[0.0, 10000.0])
    print("Composite for all lipids and both leaflets; 0-10 ns")
    print("Values from estimators:")
    print("  Basic Einstein relation:")
    print("    Diffusion coefficient: {}".format(D_e))
    print("  Linear fit:")
    print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
    print("  Anomalous diffusion fit:")
    print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))
    #diffusion coeff,  10-100 ns
    if max(times)>10000.0:
        # Simple application of Einstein relation
        D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[10000.0, 100000.0])
        # Use linear fit
        D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[10000.0, 100000.0])
        # Use anomalous diffusion fit
        D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[10000.0, 100000.0])
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
        D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[0.0, 10000.0])
        # Use linear fit
        D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[0.0, 10000.0])
        # Use anomalous diffusion fit
        D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[0.0, 10000.0])
        print("Composite for "+lipid_type+" and both leaflets; 0-10 ns")
        print("Values from estimators:")
        print("  Basic Einstein relation:")
        print("    Diffusion coefficient: {}".format(D_e))
        print("  Linear fit:")
        print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
        print("  Anomalous diffusion fit:")
        print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))
        if max(times)>10000.0:
            # diffusion coeff,  10-100 ns
            # Simple application of Einstein relation
            D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[10000.0, 100000.0])
            # Use linear fit
            D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[10000.0, 100000.0])
            # Use anomalous diffusion fit
            D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[10000.0, 100000.0])
            print("Composite for "+lipid_type+" and both leaflets; 10-100 ns")
            print("Values from estimators:")
            print("  Basic Einstein relation:")
            print("    Diffusion coefficient: {}".format(D_e))
            print("  Linear fit:")
            print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
            print("  Anomalous diffusion fit:")
            print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))
        #now do individual leaflets
        for leaflet in ['upper', 'lower']:
            msd_dat = analyzer.get_analysis_data("msd_" + lipid_type+"_"+leaflet)

            times = msd_dat[:, 0]
            msd_vals = msd_dat[:, 1]
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
            D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[0.0, 10000.0])
            # Use linear fit
            D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[0.0, 10000.0])
            # Use anomalous diffusion fit
            D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[0.0, 10000.0])
            print("Composite for " + lipid_type + " and "+leaflet+" leaflet; 0-10 ns")
            print("Values from estimators:")
            print("  Basic Einstein relation:")
            print("    Diffusion coefficient: {}".format(D_e))
            print("  Linear fit:")
            print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
            print("  Anomalous diffusion fit:")
            print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))
            if max(times) > 10000.0:
                # diffusion coeff,  10-100 ns
                # Simple application of Einstein relation
                D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[10000.0, 100000.0])
                # Use linear fit
                D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[10000.0, 100000.0])
                # Use anomalous diffusion fit
                D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[10000.0, 100000.0])
                print("Composite for " + lipid_type + " and "+leaflet+" leaflet; 10-100 ns")
                print("Values from estimators:")
                print("  Basic Einstein relation:")
                print("    Diffusion coefficient: {}".format(D_e))
                print("  Linear fit:")
                print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
                print("  Anomalous diffusion fit:")
                print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))


    return


def area_per_lipid(structure_file, trajectory_file, selection_string, frame_start=0, frame_end=-1,
                  frame_interval=1, dump_path=None):
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    #remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    #add the apl analyses
    analyzer.add_analysis("apl_box apl_box")
    analyzer.add_analysis("apl_grid apl_grid")

    analyzer.print_analysis_protocol()

    #add the plots
    analyzer.add_plot("apl apl_box apl_box None")
    analyzer.add_plot("apl apl_p apl_box Box apl_grid None")
    analyzer.add_plot("apl apl_grid apl_grid None")

    analyzer.print_plot_protocol()

    #run analysis
    analyzer.run_analysis()

    #output data and plots
    analyzer.dump_data(path=dump_path)
    analyzer.save_all_plots()

    return

def bilayer_thickness(structure_file, trajectory_file, selection_string, frame_start=0, frame_end=-1,
                  frame_interval=1, dump_path=None, name_dict=None):
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

    analyzer.print_analysis_protocol()
    #add the plot
    analyzer.add_plot("bilayer_thickness bt_p bt")

    analyzer.print_plot_protocol()

    #run analysis
    analyzer.run_analysis()

    #output data and plots
    analyzer.dump_data(path=dump_path)
    analyzer.save_all_plots()

    return

def compressibility(structure_file, trajectory_file, selection_string, frame_start=0, frame_end=-1,
                  frame_interval=1, dump_path=None, temperature=298.15):

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
    acm = analyzer.get_analysis_data('acm')
    vcm = analyzer.get_analysis_data('vcm')
    ac = analyzer.get_analysis_data('ac')
    print(acm)
    print("Area compressibility modulus: {} mN/m".format(acm[-1:,1]))
    print("Area compressibility: {} m/mN".format(ac[-1:,1]))
    print("Volume compressibility modulus: {} J/Angstrom^3".format(vcm[-1:,1]))

    return

def dispvector_correlation(structure_file, trajectory_file, selection_string, frame_start=0, frame_end=-1,
                  frame_interval=1, dump_path=None):
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    #remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    #add the apl analyses
    analyzer.add_analysis("apl_box apl_box")
    analyzer.add_analysis("apl_grid apl_grid")

    analyzer.print_analysis_protocol()

    #add the plots
    analyzer.add_plot("apl apl_box apl_box None")
    analyzer.add_plot("apl apl_p apl_box Box apl_grid None")
    analyzer.add_plot("apl apl_grid apl_grid None")

    analyzer.print_plot_protocol()

    #run analysis
    analyzer.run_analysis()

    #output data and plots
    analyzer.dump_data(path=dump_path)
    analyzer.save_all_plots()

    return