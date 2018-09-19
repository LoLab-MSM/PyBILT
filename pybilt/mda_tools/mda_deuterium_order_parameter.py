from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import MDAnalysis as mda
#import the running stats class
from pybilt.common.running_stats import RunningStats
from six.moves import range

def vector_to_unit(v):
    #returns the unit vector version of v: v_u = v/|v|
    # magnitude of v
    v_mag = np.sqrt(np.dot(v, v))
    return v / v_mag

def angle_between_vectors(v1, v2):
    # computes the angle (radians) between vectors v1 and v2
    # dot(v1, v2) = |v1||v2|cos(theta), so theta = acos[ dot(v1, v2)/(|v1||v2|) ]
    #convert input vectors to unit vectors -- takes care of magintudes
    v1_u = vector_to_unit(v1)
    v2_u = vector_to_unit(v2)
    # dot product
    dot = np.dot(v1_u,v2_u)
    # clip the dot product to ensure it is in bounds for arccos
    dot_clip = np.clip(dot, -1.0, 1.0)
    # get the angle
    theta = np.arccos(dot_clip)
    return theta

def cos_angle_between_vectors(v1, v2):
    # computes the cosine of the angle between vectors v1 and v2
    # dot(v1, v2) = |v1||v2|cos(theta), so cos(theta) = dot(v1, v2)/(|v1||v2|)
    #convert input vectors to unit vectors -- takes care of magintudes
    v1_u = vector_to_unit(v1)
    v2_u = vector_to_unit(v2)
    # dot product
    dot = np.dot(v1_u,v2_u)
    # clip the dot product to ensure it is in bounds
    dot_clip = np.clip(dot, -1.0, 1.0)
    return dot_clip


def build_acyl_index_lists(membrane_lipid_sel):
    acyl_carbons = []
    acyl_hydrogens = []
    for atom in membrane_lipid_sel.atoms:
        #check if carbon
        if atom.name[0] is 'C' or atom.type[0] is 'C':
            nb = len(atom.bonds)
            if nb == 4:
                hbond = []
                cbond = []
                nch_flag = True
                for bond in atom.bonds:
                    for batom in bond.atoms:
                        btn=batom.name[0]
                        btt=batom.type[0]
                        nch_h = btn == 'H' or btt == 'H'
                        nch_c = btn == 'C'or btt == 'C'
                        #print "btn ",btn," btt ",btt
                        #print "nch_h ",nch_h," nch_c ",nch_c
                        #make sure the carbon is only bonded to hydrogens other carbon
                        if nch_flag:
                            nch_flag = nch_h or nch_c
                            #print "nch_flag ",nch_flag
                        #if the current bond is to hydrogen add it to list
                        if nch_h:
                            hbond.append(batom.index)
                #passes if acyl carbon
                #print "nch_flag ",nch_flag
                if nch_flag and len(hbond) == 2:
                    acyl_carbons.append(atom.index)
                    acyl_hydrogens.append(hbond)
    return acyl_carbons,acyl_hydrogens

#helper functions

def _get_bilayer_norm_vector(norm_axis):
    if norm_axis is 'z':
        return np.array([0.0,0.0,1.0])
    elif norm_axis is 'x':
        return np.array([1.0,0.0,0.0])
    elif norm_axis is 'y':
        return np.array([0.0,1.0,0.0])

def _adjust_frame_range_for_slicing(fstart, fend, nframes):
    if fend != -1:
        fend+=1
    if fend == (nframes-1) and fstart == (nframes-1):
        fend+=1
    if fend == fstart:
        fend+=1
    if fstart<0:
        fstart+=nframes
    if fend < 0:
        fend+=nframes+1
    return fstart, fend

#average over all acyl groups of all the lipids based on description in Moore et al. 2001 Biophysical Journal 81(5) 2484-2494
def average_deuterium_order_Moore(trajectory,membrane_sel, fstart=0,fend=-1,fstep=1, norm_axis='z', make_whole=False):


    bilayer_norm = _get_bilayer_norm_vector(norm_axis)
    nframes = len(trajectory)

    #adjust the frame end points for slicing
    fstart, fend = _adjust_frame_range_for_slicing(fstart, fend, len(trajectory))
    nframes = (fend-fstart)/fstep + 1
    print(("doing frame slice points {} to {} with step/interval {}".format(fstart, fend, fstep)))
    print(("total of {} frames".format(nframes)))
    #build the index lists of acyl components
    print("building index lists for acyl groups")
    acyl_carbons,acyl_hydrogens = build_acyl_index_lists(membrane_sel)
    print(("there are {} acyl groups".format(len(acyl_carbons))))
    #configuration and time average for Scd = < 0.5 ( 3 cos**2(beta) - 1) >
    Scd = RunningStats()
    Scd_out = np.zeros((nframes,6))
    Scd_all = RunningStats()
    f=0
    for frame in trajectory[fstart:fend:fstep]:
        if make_whole:
            for res in membrane_sel.residues:
                #print(membrane_sel)
                #print(res)
                mda.lib.mdamath.make_whole(res.atoms)
        Scd_i = RunningStats()
        for i in range(len(acyl_carbons)):

            c_pos = frame.positions[acyl_carbons[i]]
            h1_pos = frame.positions[acyl_hydrogens[i][0]]
            h2_pos = frame.positions[acyl_hydrogens[i][1]]
            acyl_norm = np.cross(h1_pos - c_pos, h2_pos - c_pos)
            #patch for broken bonds due to pbc wrapping
            dh1 = h1_pos - c_pos
            dh2 = h2_pos - c_pos
            rh1 = np.sqrt(np.dot(dh1, dh1))
            rh2 = np.sqrt(np.dot(dh1, dh1))
            if (rh1 < 5.0) and (rh2 < 5.0):
                Scd_i.push( 0.50*(3.0* cos_angle_between_vectors(acyl_norm, bilayer_norm)**2 -1.0) )
                Scd_all.push( 0.50*(3.0* cos_angle_between_vectors(acyl_norm, bilayer_norm)**2 -1.0) )
        Scd_f = Scd_i.mean()
        Scd.push(Scd_f)
        Scd_out[f,0]=frame.time
        Scd_out[f, 1] = Scd_f
        Scd_out[f,2]=Scd.mean()
        Scd_out[f,3]=Scd.deviation()
        Scd_out[f,4]=Scd_all.mean()
        Scd_out[f,5]=Scd_all.deviation()
        f+=1
    return Scd_out

#based on description in
def deuterium_order_parameter(mda_universe, lipid_segids, lipid_acyl_carbons,
                              lipid_acyl_hydrogens, first_frame=0,
                              last_frame=-1, frame_interval=1, norm_axis='z'):
    """Compute the deteurium (or acyl chain) order parameter
    This function computes the average deuterium (or acyl chain) order
    parameter for acyl chain carbons from the trajectory data in a MDAnalysis
    universe object. The lipid are assumed to be in a planar bilayer or
    monolayer and aligned such that the normal to the lipid layer surface
    points along either the x, y, or z axis (z by default). This function
    does not directly account for split molecules due to wrapping. As a patch,
    the function tests the C-H bond distances each frame and if they large,
    suggesting wrapped molecule splitting, that carbon is ignored that frame
    and not added to the average.
    Args:
        mda_universe (MDAnalysis.Universe): The MDAnalysis.Universe object
            containing the structure and trajectory of the lipid system to
            be analyzed.
        lipid_segids (dict): A dictionary keyed by the lipid RESNAMEs and
            and valued with the SEGID for that lipid type.
        lipid_acyl_carbons (dict): A dictionary keyed by the lipid RESNAMEs.
            The values are dictionaries keyed by integer acyl chain carbon
            positions that you want to include in the analysis for that lipid
            type and has values of the acyl chain carbon name. This variable is
            used to make appropriate atom selections of the acyl chain carbons.
            e.g.: lipid_acyl_carbons = {'POPC':{1:'CA1', 2:'CA2', ...}, ...}
        lipid_acyl_hydrogens (dict): A dictionary keyed by the lipid RESNAMEs.
            The values are dictionaries keyed by integer acyl chain carbon
            positions that you want to include in the analysis for that lipid
            type and has values of a list containing the names of the two
            hydrogens for that acyl chain carbon position. This variable is
            used to make appropriate atom selections of the acyl chain
            hydrogens.
            e.g.: lipid_acyl_hydrogens = {'POPC':{1:['H1A', 'H1B'], ...}, ...}
        first_frame (int, optional): Define the first frame of the trajectory
            from which to start the analysis. Default: 0, the first frame in
            mda_universe.trajectory
        last_frame (int, optional): Define the last frame of the trajectory
            to include in the analysis. Default: -1, the last frame in
            mda_universe.trajectory
        frame_interval (int, optional): Define an interval for looping over
            the frames in the trajectory. Default: 1, every frame
        norm_axis (str): Define the axis that is normal to the lipid layer
            surface. Options: 'x', 'y', or 'z'. Default: 'z'.

        Returns:
            (tuple): 2 elements - element 0: (dict): A dictionary with same key
                structure as lipid_acyl_carbons but with interior values
                containing the ensemble+time averaged order parameter values.
                element 1: (dict): Same as first element but containing the
                standard deviations of the averages.

        References:
            1. Moore et al. 2001 Biophysical Journal 81(5) 2484-2494
            2. Vermeer Eur Biophys J (2007) 36:919-931
    """
    trajectory = mda_universe.trajectory
    bilayer_norm = _get_bilayer_norm_vector(norm_axis)
    nframes = len(trajectory)
    #setup the output container
    averagers = {}
    for lipid_type in lipid_acyl_carbons.keys():
        averagers[lipid_type] = {}
        for number in lipid_acyl_carbons[lipid_type]:
            averagers[lipid_type][number] = RunningStats()

    # adjust the frame end points for slicing
    first_frame, last_frame = _adjust_frame_range_for_slicing(first_frame, last_frame, nframes)

    nframes = (last_frame - first_frame) / frame_interval
    print(("doing frame slice points {} to {} with step/interval {}".format(first_frame, last_frame, frame_interval)))
    print(("total of {} frames".format(nframes)))

    for frame in trajectory[first_frame:last_frame:frame_interval]:
        for lipid_type in lipid_acyl_carbons.keys():
            segid = lipid_segids[lipid_type]
            for position in lipid_acyl_carbons[lipid_type].keys():
                carbon_name = lipid_acyl_carbons[lipid_type][position]
                hydrogen_names = lipid_acyl_hydrogens[lipid_type][position]
                Cs_str = "segid {} and resname {} and name {}".format(segid, lipid_type, carbon_name)
                Cs = mda_universe.select_atoms(Cs_str)
                H1s_str = "segid {} and resname {} and name {}".format(segid, lipid_type, hydrogen_names[0])
                H1s = mda_universe.select_atoms(H1s_str)
                H2s_str = "segid {} and resname {} and name {}".format(segid, lipid_type, hydrogen_names[1])
                H2s = mda_universe.select_atoms(H2s_str)
                n_lipid = len(Cs)
                for i in range(n_lipid):
                    v_ch1 = H1s[i].position - Cs[i].position
                    v_ch2 = H2s[i].position - Cs[i].position
                    cos_beta_ch1 = cos_angle_between_vectors(v_ch1,
                                                             bilayer_norm)
                    cos_beta_ch2 = cos_angle_between_vectors(v_ch2,
                                                             bilayer_norm)
                    scd_ch1 = 0.50*(3.0* cos_beta_ch1**2 -1.0)
                    scd_ch2 = 0.50*(3.0* cos_beta_ch1**2 -1.0)
                    # Length condition is a patch for broken bonds due to pbc wrapping
                    rh1 = np.sqrt(np.dot(v_ch1, v_ch1))
                    rh2 = np.sqrt(np.dot(v_ch2, v_ch2))
                    if (rh1 < 5.0) and (rh2 < 5.0):
                        # Take the average of the two hydrogens to get
                        # the instantaneous value for this carbon
                        scd_c = (scd_ch1 + scd_ch2)/2.0
                        averagers[lipid_type][position].push(scd_c)
    averages = dict()
    stds = dict()
    for lipid_type in lipid_acyl_carbons.keys():
        averages[lipid_type] = dict()
        stds[lipid_type] = dict()
        for position in lipid_acyl_carbons[lipid_type].keys():
            averages[lipid_type][position] = averagers[lipid_type][position].mean()
            stds[lipid_type][position] = averagers[lipid_type][position].deviation()
    return averages, stds

#average over all acyl groups of all the lipids based on description in Vermeer Eur Biophys J (2007) 36:919-931
def average_deuterium_order_Vermeer(trajectory,membrane_sel, fstart=0,fend=-1,fstep=1, norm_axis='z'):

    bilayer_norm = _get_bilayer_norm_vector(norm_axis)
    nframes = len(trajectory)

    #adjust the frame end points for slicing
    fstart, fend = _adjust_frame_range_for_slicing(fstart, fend, nframes)

    nframes = (fend - fstart)/fstep
    print(("doing frame slice points {} to {} with step/interval {}".format(fstart, fend, fstep)))
    print(("total of {} frames".format(nframes)))

    #build the index lists of acyl components
    print("building index lists for acyl groups")
    acyl_carbons,acyl_hydrogens = build_acyl_index_lists(membrane_sel)
    print("there are ",len(acyl_carbons)," acyl groups")
    #configuration and time average for Scd = < 0.5 ( 3 cos**2(beta) - 1) >
    Scd = RunningStats()
    # Scd_out = np.zeros((nframes,3))
    Scd_out = np.zeros((nframes,6))
    Scd_all = RunningStats()
    f=0
    for frame in trajectory[fstart:fend:fstep]:
        curr_time = frame.time
        Scd_i = RunningStats()
        for i in range(len(acyl_carbons)):

            c_i = acyl_carbons[i]
            h1_i = acyl_hydrogens[i][0]
            h2_i = acyl_hydrogens[i][1]
            c_pos = frame.positions[c_i]
            h1_pos = frame.positions[h1_i]
            h2_pos = frame.positions[h2_i]
            v_ch1 = h1_pos - c_pos
            v_ch2 = h2_pos - c_pos
            cos_beta_ch1 = cos_angle_between_vectors(v_ch1,bilayer_norm)
            cos_beta_ch2 = cos_angle_between_vectors(v_ch2,bilayer_norm)
            scd_ch1 = 0.50*(3.0* cos_beta_ch1**2 -1.0)
            scd_ch2 = 0.50*(3.0* cos_beta_ch1**2 -1.0)
            # Length condition is a patch for broken bonds due to pbc wrapping
            rh1 = np.sqrt(np.dot(v_ch1, v_ch1))
            rh2 = np.sqrt(np.dot(v_ch2, v_ch2))
            if (rh1 < 5.0) and (rh2 < 5.0):
                Scd_i.push(scd_ch1)
                Scd_i.push(scd_ch2)
                Scd_all.push(scd_ch1)
                Scd_all.push(scd_ch2)
        Scd_f = Scd_i.mean()
        Scd.push(Scd_f)
        Scd_out[f,0]=frame.time
        Scd_out[f, 1] = Scd_f
        Scd_out[f,2]=Scd.mean()
        Scd_out[f,3]=Scd.deviation()
        Scd_out[f,4]=Scd_all.mean()
        Scd_out[f,5]=Scd_all.deviation()
        f+=1
    return Scd_out
