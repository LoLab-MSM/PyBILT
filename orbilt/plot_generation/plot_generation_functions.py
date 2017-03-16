'''
A set of functions to generate plots/figures from the lipid bilayer analysis outputs.
These functions use matplotlib (http://matplotlib.org/index.html) along with Seaborn (
https://stanford.edu/~mwaskom/software/seaborn/index.html). 

'''
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import seaborn as sns
import numpy as np


# the default savefig params can be different from the display params
# e.g., you may want a higher resolution, or to make the figure
# background white
sfig_params = { 
    'savefig.dpi' : 300,
    'savefig.format' : 'eps'
    }
mpl.rcParams.update(sfig_params)
params = {'figure.figsize': [8.5, 6.0], 'font.size': 14, 'axes.labelsize': 18,
    'legend.fontsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,}
mpl.rcParams.update(params)
#sns.set_style("whitegrid")
#sns.set_style("white")
#sns.set(context="paper", font="monospace")
sns.set_style("ticks")


def plot_step_vectors(vectors_resnames, filename='step_vectors.eps',save=True, show=False):
    '''
    Generates a single plot with the lipid displacement vectors (or step vectors)
    Takes a single frame of the output from:
        MemSys.StepVector
    Corresponding colors (if multiple lipid types are included) can be 
    generated using:
        MemSys.StepVectorColors
    '''
    color_list = ['red','green', 'blue','black','orange','purple', 'yellow']
    sns.set_style("whitegrid")
    x = vectors_resnames[0][:,0]
    y=vectors_resnames[0][:,1]
    vx = vectors_resnames[0][:,2]
    vy = vectors_resnames[0][:,3]
    step_vec_plot = plt.figure()
    resnames = set(vectors_resnames[1])
    i = 0
    color_dict = {}
    for res in resnames:
        color_dict[res] = color_list[i]
        i+=1
        if i > len(color_list)-1:
            i=0
    colors = []
    for residue in vectors_resnames[1]:
        colors.append(color_dict[residue])

    Q = plt.quiver(x,y,vx,vy,color=colors)
    label_string = ""
    for resname in color_dict:
        label_string+=resname+":"+color_dict[resname]+" "
    qk = plt.quiverkey(Q, 0.25, 0.95, 2, label_string, labelpos='E',
                   coordinates='figure')
    #else:
    #    plt.quiver(x,y,vx,vy)
    #plt.title('Lateral Displacement Vectors')
    if save:
        plt.savefig(filename)
    if show:
        return plt.show()
    plt.close()
    return

def plot_msd(msd_dat_list,name_list=None,filename='msd.eps',time_in='ps',time_out='ns',show=False, interval=1,save=True):
    '''
    Generates a single plot with Mean Squared Displacement curves
    Takes outputs from:
        MemSys.CalcMSD 
        MemSys.CalcMSD_parallel
    The outputs are passed to function in a list input: apl_dat_list
    '''
#    params = {
#    'axes.labelsize': 20,
#    'text.fontsize': 20,
#    'legend.fontsize': 20,
#    'xtick.labelsize': 16,
#    'ytick.labelsize': 16,
#    'text.usetex': False,
#    'figure.figsize': [8.0, 6.0]
#    }
#    params = {'figure.figsize': [10.0, 8.0]}
#    mpl.rcParams.update(params)
#        
    i = 0
    for msd_dat in msd_dat_list:
        msd_d = msd_dat.copy()
        t = msd_d[::interval,0]
        if time_in == 'ps' and time_out == 'ns':
            t/=1000.0
        elif time_in == 'ns' and time_out == 'ps':
            t*=1000.0
        msd = msd_d[::interval,1]
        if name_list is not None:
            plt.plot(t, msd, linewidth=2.0,label=name_list[i])
        else:
            plt.plot(t, msd, linewidth=2.0)
        i+=1
        #plt.title("Mean Sqared Displacement vs. Time")
    xlabel = "Time ("+time_out+")"
    plt.xlabel(xlabel)
    plt.ylabel("Distance in the lateral plane ($\AA^2$)")
    if name_list is not None:
        plt.legend(loc=0)
    if save:    
        plt.savefig(filename)
    if show:
        return plt.show()
    plt.close()
    return


def plot_area_per_lipid(apl_dat_list,name_list=None,filename='apl.eps',time_in='ps',time_out='ns',save=True,show=False, interval=1):
    '''
    Generates a single plot with area per lipid (apl) curves
    Takes outputs from:
        MemSys.CalcAreaPerLipid_Box 
        MemSys.CalcAreaPerLipid_ClosestNeighborCircle
    The outputs are passed to function in a list input: apl_dat_list
    '''
    #print "filename: ", filename
    i = 0
    for apl_dat in apl_dat_list:
        apl_d = apl_dat.copy()
      #  print "apl_d"
      #  print apl_d
        t = apl_d[::interval,0]
        if time_in == 'ps' and time_out == 'ns':
            #print "switching time units from ps to ns"
            t/=1000.0
        elif time_in == 'ns' and time_out == 'ps':
            t*=1000.0
        apl = apl_d[::interval,2]
        apl_dev = apl_d[::interval,3]
        if name_list is not None:
            #print "plotting",name_list[i]," with errorbars"
            #print t
            #print apl
            plt.errorbar(t, apl, yerr=apl_dev,linewidth=2.0,label=name_list[i])
        else:
            plt.errorbar(t, apl, yerr=apl_dev,linewidth=2.0)
        i+=1
        #plt.title("Mean Sqared Displacement vs. Time")
    xlabel = "Time ("+time_out+")"
    plt.xlabel(xlabel)
    plt.ylabel("Area per lipid ($\AA^2$)")
    if name_list is not None:
        plt.legend(loc=0)
    if save:    
        plt.savefig(filename)
    if show:
        return plt.show()
    plt.close()
    return


def plot_cluster_dat_number(clust_dat_list,name_list=None,filename='clust_number.eps',time_in='ps',time_out='ns',show=False):
    '''
    Generates a single of the average number of clusters (vs. time)
    using output data from:
        MemSys.CheckClustering     
    The outputs are passed to function in a list input: clust_dat_list
    '''
    i = 0
    for cl_dat in clust_dat_list:
        cl_loc = cl_dat.copy()
        t = cl_loc[:,0]
        if time_in == 'ps' and time_out == 'ns':
            #print "switching time units from ps to ns"
            t/=1000.0
        elif time_in == 'ns' and time_out == 'ps':
            t*=1000.0
        cl = cl_loc[:,5]
        cl_dev = cl_loc[:,6]
        if name_list is not None:
            #print "plotting",name_list[i]," with errorbars"
            #print t
            #print apl
            plt.errorbar(t, cl, yerr=cl_dev,linewidth=2.0,label=name_list[i])
        else:
            plt.errorbar(t, cl, yerr=cl_dev,linewidth=2.0)
        i+=1
        #plt.title("Mean Sqared Displacement vs. Time")
    xlabel = "Time ("+time_out+")"
    plt.xlabel(xlabel)
    plt.ylabel("Average Number of Clusters")
    if name_list is not None:
        plt.legend(loc=0)
        
    plt.savefig(filename)
    if show:
        return plt.show()
    plt.close()
    return

def plot_cluster_dat_size(clust_dat_list,name_list=None,filename='clust_size.eps',time_in='ps',time_out='ns',show=False):
    '''
    Generates a single plot of the average cluster size (vs time)
    using output data from:
        MemSys.CheckClustering     
    The outputs are passed to function in a list input: clust_dat_list
    '''
    i = 0
    for cl_dat in clust_dat_list:
        cl_loc = cl_dat.copy()
        t = cl_loc[:,0]
        if time_in == 'ps' and time_out == 'ns':
            #print "switching time units from ps to ns"
            t/=1000.0
        elif time_in == 'ns' and time_out == 'ps':
            t*=1000.0
        cl = cl_loc[:,7]
        cl_dev = cl_loc[:,8]
        if name_list is not None:
            #print "plotting",name_list[i]," with errorbars"
            #print t
            #print apl
            plt.errorbar(t, cl, yerr=cl_dev,linewidth=2.0,label=name_list[i])
        else:
            plt.errorbar(t, cl, yerr=cl_dev,linewidth=2.0)
        i+=1
        #plt.title("Mean Sqared Displacement vs. Time")
    xlabel = "Time ("+time_out+")"
    plt.xlabel(xlabel)
    plt.ylabel("Average Size of Cluster (lipids per cluster)")
    if name_list is not None:
        plt.legend(loc=0)
        
    plt.savefig(filename)
    if show:
        return plt.show()
    plt.close()
    return


def plot_cluster_maps(clusters, filename='cluster_map.eps',show=False):
    '''
    Generates a single plot of the lipid cluster map 
    Takes a single frame of the output from:
        MemSys.ExportClustersForPlotting
    '''
    sns.set_style("whitegrid")
    x = clusters[0]
    y=clusters[1]
    c = clusters[2]
    plt.scatter(x,y,c=c,s=800)
    #plt.title('Lateral Displacement Vectors')
    plt.savefig(filename)
    if show:
        return plt.show()
    plt.close()
    return

def plot_density_profile(dp_out_list, save=True, filename='density_profile.eps', show=False, label_list=None, ylabel='Density'):
    """ Plot density profiles
    This function can be used to plot the results of density profiles functions 
    in the mda_density_profile module.
    
    Args:
        dp_out_list (list of tuples): A list of the tuple outputs of the profile calculation functions
        save (bool, optional): Default is True. Saves the plot output as an image file if True.
        filename (str, optional): The name out the image file that will be created if save=True.
        show (bool, optional): Default is False. Display the plot (plt.show) if True.
        label_list (list of str : None, optional): Default is None. Allows a list of strings used to 
            label the plot lines.
    """
    i = 0
    for item in dp_out_list:
        if label_list is not None:

            plt.plot(item[0], item[1], linewidth=2.0, label=label_list[i])
        else:
            plt.plot(item[0], item[1], linewidth=2.0)
        i+=1
    if label_list is not None:
        plt.legend(loc=0)
    plt.ylabel(ylabel)
    plt.xlabel('Position Along the Normal')   
    #plt.xlabel(ylabel)
    #plt.ylabel('Position Along the Normal')  
    if save:
        plt.savefig(filename)
    if show:
        return plt.show()
    plt.close()
    return


def plot_grid_as_scatter(in_xyzc, save=True, filename='lipid_grid.eps', show=False, colorbar=False, cmap=None, vmin=None, vmax=None):
    cma = plt.cm.get_cmap('viridis')
    if cmap is not None:
        cma = cmap
   # fig = plt.figure()
   # ax = fig.add_subplot(111)
   ## dx = np.abs(in_xyzc[0][1]-in_xyzc[0][0])
    #dy = np.abs(in_xyzc[1][1]-in_xyzc[1][0])
    #ds = [dx]
    #if dy > dx:
    #    ds = [dy]  
     
    #ds_in_points = np.diff(ax.transData.transform(zip([0]*1, ds)))[0] 
    #print ds_in_points
    if vmin is not None and vmax is None:
        plt.scatter(in_xyzc[0], in_xyzc[1], c=in_xyzc[3], marker='s',s=50, cmap=cma, vmin=vmin)
    elif vmax is not None and vmin is None:
        plt.scatter(in_xyzc[0], in_xyzc[1], c=in_xyzc[3], marker='s',s=50, cmap=cma, vmax=vmax)
    elif vmin is not None and vmax is not None:
        plt.scatter(in_xyzc[0], in_xyzc[1], c=in_xyzc[3], marker='s', s=50, cmap=cma, vmin=vmin, vmax=vmax)
    else:
        plt.scatter(in_xyzc[0], in_xyzc[1], c=in_xyzc[3], marker='s',s=50, cmap=cma)
    #print in_xyzc[3]
    #plt.scatter(in_xyzc[0], in_xyzc[1], c=in_xyzc[3], marker='s',s=50, cmap=cma)
    #cax, kw = mpl.colorbar.make_axes(plt.gca())
    #norm = mpl.colors.Normalize(vmin = min(in_xyzc[3]), vmax = max(in_xyzc[3]), clip = False)

    #c = mpl.colorbar.ColorbarBase(cax, cmap=cma, norm=norm)
    if colorbar:
        plt.colorbar()
    if save:
        plt.savefig(filename)
    if show:
        return plt.show()
    plt.close()
    return

def plot_average_deuterium_op(dop_dat_list,name_list=None,filename='dop.eps',time_in='ps',time_out='ns',show=False, interval=1):
    '''
    Generates a single plot of the average deuterium order parameter vs. time

    The outputs are passed to function in a list input: dop_dat_list
    '''
    #print "filename: ", filename
    i = 0
    for dop_dat in dop_dat_list:
        dop_d = dop_dat.copy()
        t = dop_d[::interval,0]
        if time_in == 'ps' and time_out == 'ns':
            #print "switching time units from ps to ns"
            t/=1000.0
        elif time_in == 'ns' and time_out == 'ps':
            t*=1000.0
        dop = dop_d[::interval,1]
        dop_dev = dop_d[::interval,2]
        if name_list is not None:
            #print "plotting",name_list[i]," with errorbars"
            #print t
            #print dop
            plt.errorbar(t, dop, yerr=dop_dev,linewidth=2.0,label=name_list[i])
        else:
            plt.errorbar(t, dop, yerr=dop_dev,linewidth=2.0)
        i+=1
        #plt.title("Mean Sqared Displacement vs. Time")
    xlabel = "Time ("+time_out+")"
    plt.xlabel(xlabel)
    plt.ylabel("Average Deuterium Order Parameter")
    if name_list is not None:
        plt.legend(loc=0)
        
    plt.savefig(filename)
    if show:
        return plt.show()
    plt.close()
    return

def plot_bilayer_thickness(bt_dat_list,name_list=None,filename='bilayer_thickness.eps',time_in='ps',time_out='ns',show=False, interval=1,save=True):
    '''
    Generates a single plot with bilayer thickness curves
    Takes outputs from:
        
    The outputs are passed to function in a list input: bt_dat_list
    '''
#    params = {
#    'axes.labelsize': 20,
#    'text.fontsize': 20,
#    'legend.fontsize': 20,
#    'xtick.labelsize': 16,
#    'ytick.labelsize': 16,
#    'text.usetex': False,
#    'figure.figsize': [8.0, 6.0]
#    }
#    params = {'figure.figsize': [10.0, 8.0]}
#    mpl.rcParams.update(params)
#        
    i = 0
    for bt_dat in bt_dat_list:
        bt_d = bt_dat.copy()
        t = bt_d[::interval,0]
        if time_in == 'ps' and time_out == 'ns':
            t/=1000.0
        elif time_in == 'ns' and time_out == 'ps':
            t*=1000.0
        bt = bt_d[::interval,1]
        error = bt_d[::interval,2]
        if name_list is not None:
            plt.errorbar(t, bt, yerr=error,linewidth=2.0,label=name_list[i])
        else:
            plt.errorbar(t, bt, yerr=error,linewidth=2.0)
        i+=1
        #plt.title("Mean Sqared Displacement vs. Time")
    xlabel = "Time ("+time_out+")"
    plt.xlabel(xlabel)
    plt.ylabel("Bilayer thickness ($\AA$)")
    if name_list is not None:
        plt.legend(loc=0)
    if save:    
        plt.savefig(filename)
    if show:
        return plt.show()
    plt.close()
    return