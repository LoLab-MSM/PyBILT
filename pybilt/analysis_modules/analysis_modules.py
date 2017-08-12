import numpy as np
import MDAnalysis as mda

from pybilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer
from pybilt.diffusion import diffusion_coefficients as dc
import pybilt.plot_generation.plot_generation_functions as pgf
from pybilt.plot_generation.plot_generation_functions import _color_list
from pybilt.common.running_stats import BlockAverager
from pybilt.mda_tools.mda_density_map import position_density_map_2d_leaflet_simple

def msd_diffusion(structure_file, trajectory_file, selection_string, resnames=None, frame_start=0,
                  frame_end=-1, frame_interval=1, dump_path=None):

    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    #for each lipid resname in the input resnames add an msd analysis for both, upper, and lower leaflets
    if resnames is None:
        resnames = []
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
    t0 = times[0]
    t10 = t0+10000.0
    t100 = t0+100000.0
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
    D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[t0, t10])
    # Use linear fit
    D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[t0, t10])
    # Use anomalous diffusion fit
    D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[t0, t10])
    print("Composite for all lipids and both leaflets; 0-10 ns")
    print("Values from estimators:")
    print("  Basic Einstein relation:")
    print("    Diffusion coefficient: {}".format(D_e))
    print("  Linear fit:")
    print("    Diffusion coefficient: {} Std Error: {}".format(D_l[0], D_l[1]))
    print("  Anomalous diffusion fit:")
    print("    Diffusion coefficient: {} Alpha value: {}".format(D_a[0], D_a[1]))
    #diffusion coeff,  10-100 ns
    if max(times)>t10:
        # Simple application of Einstein relation
        D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[t10, t100])
        # Use linear fit
        D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[t10, t100])
        # Use anomalous diffusion fit
        D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[t10, t100])
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
        D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[t0, t10])
        # Use linear fit
        D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[t0, t10])
        # Use anomalous diffusion fit
        D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[t0, t10])
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
            D_e = dc.diffusion_coefficient_Einstein(times, msd_vals, time_range=[t10, t100])
            # Use linear fit
            D_l = dc.diffusion_coefficient_linear_fit(times, msd_vals, time_range=[t10, t100])
            # Use anomalous diffusion fit
            D_a = dc.diffusion_coefficient_anomalous_fit(times, msd_vals, time_range=[t10, t100])
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


def area_per_lipid(structure_file, trajectory_file, selection_string, frame_start=0, frame_end=-1,
                  frame_interval=1, dump_path=None, n_xbins=100, n_ybins=100):
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

    #adjust the number of bins for gridding
    analyzer.rep_settings['lipid_grid']['n_xbins'] = n_xbins
    analyzer.rep_settings['lipid_grid']['n_ybins'] = n_ybins

    analyzer.print_plot_protocol()

    #run analysis
    analyzer.run_analysis()

    #output data and plots
    analyzer.dump_data(path=dump_path)
    analyzer.save_all_plots()
    #print final ensemble averages to stdout
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

def bilayer_thickness(structure_file, trajectory_file, selection_string, frame_start=0, frame_end=-1,
                  frame_interval=1, dump_path="./", name_dict=None, n_xbins=100, n_ybins=100, plot_grid_map=True):
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
    #adjust the number of bins for gridding
    analyzer.rep_settings['lipid_grid']['n_xbins'] = n_xbins
    analyzer.rep_settings['lipid_grid']['n_ybins'] = n_ybins

    analyzer.print_analysis_protocol()
    #add the plot
    analyzer.add_plot("bilayer_thickness bt_p bt None")

    analyzer.print_plot_protocol()

    #run analysis
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

    #output data and plots
    analyzer.dump_data(path=dump_path)
    analyzer.save_all_plots()
    #output final ensemble average to stdout
    bt = analyzer.get_analysis_data('bt')
    print("Bilayer thickness from gridding procedure (Angstrom): {:0.4f} +- {:0.4f}".format(bt[-1][2], bt[-1][3]))

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
    analyzer.dump_data(path=dump_path)
    acm = analyzer.get_analysis_data('acm')
    vcm = analyzer.get_analysis_data('vcm')
    ac = analyzer.get_analysis_data('ac')
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
    #compute the displacment vectors for maps
    analyzer.add_analysis("disp_vec disp_vec_upper scale True wrapped True leaflet upper interval "+str(frame_interval))
    analyzer.add_analysis("disp_vec disp_vec_lower scale True wrapped True leaflet lower interval "+str(frame_interval))
    #compute the full correlation matrix between displacement vectors (i.e. cos(theta))
    analyzer.add_analysis("disp_vec_corr disp_vec_corr interval "+str(frame_interval))
    #comput the correlations between a displacement vector and that lipids closest neighbor in the lateral dimensions
    analyzer.add_analysis("disp_vec_nncorr disp_vec_nncorr_upper leaflet upper interval "+str(frame_interval))
    analyzer.add_analysis("disp_vec_nncorr disp_vec_nncorr_lower leaflet lower interval "+str(frame_interval))
    analyzer.print_analysis_protocol()

    #run analysis
    analyzer.run_analysis()

    #output data and plots
    analyzer.dump_data(path=dump_path)


    #generate the plots/maps for displacement vectors
    disp_vecs = analyzer.get_analysis_data('disp_vec_upper')
    counter=0
    number = str(len("{}".format(len(disp_vecs)) )+1)
    form = "{:0"+number+"d}"
    for disp_vec in disp_vecs:
        count = form.format(counter)
        filename = "step_vector_map_upper_"+count+".pdf"
        filename_b = "step_vector_map_upper_"+count+".png"
        pgf.plot_step_vectors(disp_vec, filename=filename, scaled=True, wrapped=True)
        pgf.plot_step_vectors(disp_vec, filename=filename_b, scaled=True, wrapped=True)
        counter+=1

    disp_vecs = analyzer.get_analysis_data('disp_vec_lower')
    counter=0
    number = str(len("{}".format(len(disp_vecs)) )+1)
    form = "{:0"+number+"d}"
    for disp_vec in disp_vecs:
        count = form.format(counter)
        filename = "step_vector_map_lower_"+count+".pdf"
        filename_b = "step_vector_map_lower_"+count+".png"
        pgf.plot_step_vectors(disp_vec, filename=filename, scaled=True, wrapped=True)
        pgf.plot_step_vectors(disp_vec, filename=filename_b, scaled=True, wrapped=True)
        counter+=1

    disp_vec_corrs = analyzer.get_analysis_data('disp_vec_corr')
    counter = 0
    number = str(len("{}".format(len(disp_vecs))) + 1)
    form = "{:0" + number + "d}"
    for disp_vec_corr in disp_vec_corrs:
        corr_mat = disp_vec_corr[0]
        count = form.format(counter)
        filename = "step_vector_correlation_map_" + count + ".pdf"
        filename_b = "step_vector_correlation_map_" + count + ".png"
        pgf.plot_corr_mat_as_scatter(corr_mat, filename=filename)
        pgf.plot_corr_mat_as_scatter(corr_mat, filename=filename_b)
        counter+=1
    return

def PN_orientational_angle(structure_file, trajectory_file, selection_string, lipid_resnames, PN_names='default', frame_start=0, frame_end=-1,
                  frame_interval=1, dump_path=None):
    analyzer = BilayerAnalyzer(structure=structure_file,
                             trajectory=trajectory_file,
                             selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    #remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    #add the loa analyses
    if PN_names is 'default':
        PN_names = dict()
        for resname in lipid_resnames:
            PN_names[resname] = ('P', 'N')
    for resname in lipid_resnames:
        Pn = PN_names[resname][0]
        Nn = PN_names[resname][1]
        analyzer.add_analysis("loa loa_"+resname+"_upper leaflet upper resname "+resname+" ref_atom_1 "+Pn+" ref_atom_2 "+Nn)
        analyzer.add_analysis("loa loa_"+resname+"_lower leaflet lower resname "+resname+" ref_atom_1 "+Pn+" ref_atom_2 "+Nn)
    #comput the correlations between a displacement vector and that lipids closest neighbor in the lateral dimensions
    analyzer.print_analysis_protocol()

    #run analysis
    analyzer.run_analysis()

    #output data and plots
    analyzer.dump_data(path=dump_path)

    for resname in lipid_resnames:
        loa_upper = analyzer.get_analysis_data("loa_"+resname+"_upper")
        loa_lower = analyzer.get_analysis_data("loa_"+resname+"_lower")
        print("Lipid resname {} has average PN orientation anlge of {} in the upper leaflet".format(resname,loa_upper[-1][1]))
        complement = 90.0 - loa_upper[-1][1]
        print("    complement angle: {}".format(complement))
        print("Lipid resname {} has average PN orientation anlge of {} in the lower leaflet".format(resname,np.abs(loa_lower[-1][1])))
        complement = 90.0 - np.abs(loa_lower[-1][1])
        print("    complement angle: {}".format(complement))
        print(" ")
    return

def nearest_neighbor_fraction(structure_file, trajectory_file, selection_string, lipid_resnames, frame_start=0, frame_end=-1,
                  frame_interval=1, dump_path=None):
    analyzer = BilayerAnalyzer(structure=structure_file,
                             trajectory=trajectory_file,
                             selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    #remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    nres = len(lipid_resnames)
    pairs = []
    #for i in range(nres):
    #    pairs.append([lipid_resnames[i], lipid_resnames[i]])
    for i in range(nres):
        for j in range(nres):
            pairs.append([lipid_resnames[i], lipid_resnames[j]])
    #add the loa analyses
    for pair in pairs:
        l1 = pair[0]
        l2 = pair[1]
        analyzer.add_analysis("nnf nnf_"+l1+"_"+l2+" resname_1 "+l1+" resname_2 "+l2)
    #comput the correlations between a displacement vector and that lipids closest neighbor in the lateral dimensions
    analyzer.print_analysis_protocol()
    print(" ")

    #run analysis
    analyzer.run_analysis()

    #output data and plots
    analyzer.dump_data(path=dump_path)

    for pair in pairs:
        l1 = pair[0]
        l2 = pair[1]
        t_nnf = analyzer.get_analysis_data("nnf_"+l1+"_"+l2)
        print("Nearest neighbor fraction for lipid pair {} and {} : {:0.4f} +- {:0.4f}".format(l1,l2,t_nnf[-1][2],t_nnf[-1][3]))
    return

def normal_displacement_lipid_type_correlation(structure_file, trajectory_file, selection_string, frame_start=0, frame_end=-1,
                  frame_interval=1, dump_path=None, com_sub_selection_dict=None):

    analyzer = BilayerAnalyzer(structure=structure_file,
                             trajectory=trajectory_file,
                             selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    #remove the default msd analysis
    analyzer.remove_analysis('msd_1')

    if com_sub_selection_dict is not None:
        analyzer.rep_settings['com_frame']['name_dict'] = com_sub_selection_dict
    analyzer.add_analysis("ndcorr norm_disp_correlation")
    #comput the correlations between a displacement vector and that lipids closest neighbor in the lateral dimensions
    analyzer.print_analysis_protocol()
    print("--------")

    #run analysis
    analyzer.run_analysis()

    #output data and plots
    analyzer.dump_data(path=dump_path)

    ndcorr = analyzer.get_analysis_data('norm_disp_correlation')
    print(" ")
    print("Normal dimension displacement-lipid type cross correlation results:")
    for leaflet in ndcorr.keys():
        print("  {} leaflet:".format(leaflet))
        for lipid_resname in ndcorr[leaflet].keys():
            mean = ndcorr[leaflet][lipid_resname][-1][2]
            deviation = ndcorr[leaflet][lipid_resname][-1][3]
            print("    Lipid resname {}: {:0.4f} +- {:0.4f}".format(lipid_resname, mean, deviation))

    pgf.plot_displacement_lipid_type_cross_correlation(ndcorr, filename="ndcorr.png")
    pgf.plot_displacement_lipid_type_cross_correlation(ndcorr, filename="ndcorr.pdf")

    return


def lipid_grid_maps(structure_file, trajectory_file, selection_string, frame_start=0, frame_end=-1,
                  frame_interval=1, dump_path="./", name_dict=None, n_xbins=100, n_ybins=100, type_colors='auto'):
    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    # use a subselection of atoms instead of full lipid center of mass, if given
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    #analyzer.settings['print_interval']=1
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

    analyzer.analysis_protocol.use_objects['lipid_grid'] = True
    #adjust the number of bins for gridding
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

def distance_cutoff_clustering(structure_file, trajectory_file, selection_string, resnames, frame_start=0, frame_end=-1,
                  frame_interval=1, dump_path="./"):

    analyzer = BilayerAnalyzer(structure=structure_file,
                               trajectory=trajectory_file,
                               selection=selection_string)

    analyzer.set_frame_range(frame_start, frame_end, frame_interval)
    # remove the default msd analysis
    analyzer.remove_analysis('msd_1')
    #add the analyses
    for resname in resnames:
        add_in = "dc_cluster dc_cluster_{}_upper resname {} leaflet upper".format(resname,resname)
        analyzer.add_analysis(add_in)
        add_in = "dc_cluster dc_cluster_{}_lower resname {} leaflet lower".format(resname,resname)
        analyzer.add_analysis(add_in)

    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    analyzer.dump_data(path=dump_path)
    print("         ")
    for resname in resnames:
        results = analyzer.get_analysis_data("dc_cluster_{}_upper".format(resname))
        plotname_png = "{}dc_cluster_upper_{}_nclusters.png".format(dump_path,resname)
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

def position_density_maps_2d(structure_file, trajectory_file, bilayer_selection_string, resnames,
                             frame_start=0, frame_end=-1,frame_interval=1, dump_path="./"):

    u = mda.Universe(structure_file, trajectory_file)
    bilayer_sel = u.select_atoms(bilayer_selection_string)
    x_centers, y_centers, counts = position_density_map_2d_leaflet_simple(u, bilayer_sel, resnames,
                                                            fstart=frame_start, fend=frame_end, fstep=frame_interval,
                                                                          refsel=bilayer_sel, scale_to_max=False)
    for key in counts.keys():
        outname_eps = 'position_density_2d_{}_upper.pdf'.format(key)
        outname_png = 'position_density_2d_{}_upper.png'.format(key)
        pgf.plot_position_density_map_2d(x_centers, y_centers, counts[key]['upper'], filename=outname_eps, scaled_to_max=False,
                                         interpolation='gaussian')
        pgf.plot_position_density_map_2d(x_centers, y_centers, counts[key]['upper'], filename=outname_png, scaled_to_max=False,
                                     interpolation='gaussian')
        outname_eps = 'position_density_2d_{}_lower.pdf'.format(key)
        outname_png = 'position_density_2d_{}_lower.png'.format(key)
        pgf.plot_position_density_map_2d(x_centers, y_centers, counts[key]['lower'], filename=outname_eps, scaled_to_max=False,
                                         interpolation='gaussian')
        pgf.plot_position_density_map_2d(x_centers, y_centers, counts[key]['lower'], filename=outname_png, scaled_to_max=False,
                                     interpolation='gaussian')
    return
