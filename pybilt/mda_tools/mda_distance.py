import numpy as np
from MDAnalysis.analysis import align



def com_com_distances(universe, mda_selection_pairs, fstart=0, fend=-1, fstep=1):
    """Center of mass to Center of mass distance.
    This function computes the distance between the centers of mass between pairs of MDAnalysis atoms selections
    across the the MD trajectory.

    Args:
        universe (MDAnalysis.Universe): The MDAnalysis universe object to run the analysis on.
        mda_selection_pairs (list): A list of 2 element lists or tuples containing pairs of MDAnalsysis
           atom selection objects to compute the distance between.
        fstart (int): Optional, the first frame to include in the analysis. Default: 0 (or the first frame)
        fend (int): Optional, the last frame to include in the analysis. Default: -1 (or the last frame)
        fstep (int): Optional, the interval between frames in the analysis when looping from fstart to fend.
            Default: 1 (or every frame)

    Returns:
        (np.array), (list): Returns two outputs. The first is an Numpy array with the timeseries simulation times
            corresponding to the frames in the analysis. The second is list of Numpy arrays with the distances; the
            order in the list corresponds to the atom selection pairs in the mda_selection_pairs input.
    """

    nframes = len(universe.trajectory)
    # adjust the end point for slicing
    if fend != -1:
        fend += 1
    if fend == (nframes - 1) and fstart == (nframes - 1):
        fend += 1
    if fend == fstart:
        fend += 1
    if fstart < 0:
        fstart += nframes
    if fend < 0:
        fend += nframes + 1
    times = []
    pair_dists = []
    for pair in mda_selection_pairs:
        pair_dists.append([])
    for frame in universe.trajectory[fstart:fend:fstep]:
        times.append(frame.time)
        i = 0
        for pair in mda_selection_pairs:
            sel_1 = pair[0]
            sel_2 = pair[1]
            com_1 = sel_1.atoms.center_of_mass()
            com_2 = sel_2.atoms.center_of_mass()
            d_com = com_2 - com_1
            dist = np.sqrt(np.dot(d_com, d_com))
            pair_dists[i].append(dist)
            i+=1
    times = np.array(times)
    i=0
    for vals in pair_dists:
        pair_dists[i] = np.array(vals)
        i+=1
    return times, pair_dists

def com_com_distances_plane(universe, mda_selection_pairs, fstart=0, fend=-1, fstep=1, plane='xy'):
    """Center of mass to Center of mass distance in two dimensions (a plane).
    This function computes the distance between the centers of mass between pairs of MDAnalysis atoms selections
    across the the MD trajectory, but only uses the 2d coordinates of the specified axis plane.

    Args:
        universe (MDAnalysis.Universe): The MDAnalysis universe object to run the analysis on.
        mda_selection_pairs (list): A list of 2 element lists or tuples containing pairs of MDAnalsysis
           atom selection objects to compute the distance between.
        fstart (int): Optional, the first frame to include in the analysis. Default: 0 (or the first frame)
        fend (int): Optional, the last frame to include in the analysis. Default: -1 (or the last frame)
        fstep (int): Optional, the interval between frames in the analysis when looping from fstart to fend.
            Default: 1 (or every frame)
        plane (str): Optional, the 2d axis plane to compute the distance in. Default: 'xy' (or the xy plane)

    Returns:
        (np.array), (list): Returns two outputs. The first is an Numpy array with the timeseries simulation times
            corresponding to the frames in the analysis. The second is list of Numpy arrays with the distances; the
            order in the list corresponds to the atom selection pairs in the mda_selection_pairs input.
    """
    lat_ind = [0, 1]
    dir_ind = 2
    if plane is 'yx':
        dir_ind = 2
        lat_ind = [0, 1]
    elif plane is 'xz' or plane is 'zx':
        dir_ind = 1
        lat_ind = [0, 2]
    elif plane is 'zy' or plane is 'yz':
        dir_ind = 0
        lat_ind = [1, 2]
    #indices = mda_selection.indices

    nframes = len(universe.trajectory)
    # adjust the end point for slicing
    if fend != -1:
        fend += 1
    if fend == (nframes - 1) and fstart == (nframes - 1):
        fend += 1
    if fend == fstart:
        fend += 1
    if fstart < 0:
        fstart += nframes
    if fend < 0:
        fend += nframes + 1
    times = []
    pair_dists = []
    for pair in mda_selection_pairs:
        pair_dists.append([])
    for frame in universe.trajectory[fstart:fend:fstep]:
        times.append(frame.time)
        i = 0
        for pair in mda_selection_pairs:
            sel_1 = pair[0]
            sel_2 = pair[1]
            com_1 = sel_1.atoms.center_of_mass()
            com_2 = sel_2.atoms.center_of_mass()
            plane_val_1 = com_1[lat_ind]
            plane_val_2 = com_2[lat_ind]
            d_com = plane_val_2 - plane_val_1
            dist = np.sqrt(np.dot(d_com, d_com))
            pair_dists[i].append(dist)
            i+=1
    times = np.array(times)
    i=0
    for vals in pair_dists:
        pair_dists[i] = np.array(vals)
        i+=1
    return times, pair_dists


def com_com_distances_axis(universe, mda_selection_pairs, fstart=0, fend=-1, fstep=1, axis='z'):
    """Center of mass to Center of mass distance in one dimension (along an axis).
    This function computes the distance between the centers of mass between pairs of MDAnalysis atoms selections
    across the the MD trajectory, but only uses the 1d coordinate of the specified axis.

    Args:
        universe (MDAnalysis.Universe): The MDAnalysis universe object to run the analysis on.
        mda_selection_pairs (list): A list of 2 element lists or tuples containing pairs of MDAnalsysis
           atom selection objects to compute the distance between.
        fstart (int): Optional, the first frame to include in the analysis. Default: 0 (or the first frame)
        fend (int): Optional, the last frame to include in the analysis. Default: -1 (or the last frame)
        fstep (int): Optional, the interval between frames in the analysis when looping from fstart to fend.
            Default: 1 (or every frame)
        axis (str): Optional, the 1d axis to compute the distance in. Default: 'z' (or the z axis)

    Returns:
        (np.array), (list): Returns two outputs. The first is an Numpy array with the timeseries simulation times
            corresponding to the frames in the analysis. The second is list of Numpy arrays with the distances; the
            order in the list corresponds to the atom selection pairs in the mda_selection_pairs input.
    """
    lat_ind = [0, 1]
    dir_ind = 2
    if axis is 'x':
        dir_ind = 0
        lat_ind = [1, 2]
    elif axis is 'y':
        dir_ind = 1
        lat_ind = [0, 2]
    #indices = mda_selection.indices

    nframes = len(universe.trajectory)
    # adjust the end point for slicing
    if fend != -1:
        fend += 1
    if fend == (nframes - 1) and fstart == (nframes - 1):
        fend += 1
    if fend == fstart:
        fend += 1
    if fstart < 0:
        fstart += nframes
    if fend < 0:
        fend += nframes + 1
    times = []
    pair_dists = []
    for pair in mda_selection_pairs:
        pair_dists.append([])
    for frame in universe.trajectory[fstart:fend:fstep]:
        times.append(frame.time)
        i = 0
        for pair in mda_selection_pairs:
            sel_1 = pair[0]
            sel_2 = pair[1]
            com_1 = sel_1.atoms.center_of_mass()
            com_2 = sel_2.atoms.center_of_mass()
            norm_val_1 = com_1[dir_ind]
            norm_val_2 = com_2[dir_ind]
            dist = np.abs(norm_val_2 - norm_val_1)
            pair_dists[i].append(dist)
            i+=1
    times = np.array(times)
    i=0
    for vals in pair_dists:
        pair_dists[i] = np.array(vals)
        i+=1
    return times, pair_dists

def com_com_distances_axis_align(universe, mda_selection_pairs, align_struct_universe, align_sel_string, fstart=0,
                                 fend=-1, fstep=1, axis='z'):
    """Center of mass to Center of mass distance in one dimension (along an axis) after structure alignment.
    This function computes the distance between the centers of mass between pairs of MDAnalysis atoms selections
    across the the MD trajectory, but only uses the 1d coordinate of the specified axis, and aligns the structure
    to some selection of atoms of a reference structure.

    Args:
        universe (MDAnalysis.Universe): The MDAnalysis universe object to run the analysis on.
        mda_selection_pairs (list): A list of 2 element lists or tuples containing pairs of MDAnalsysis
           atom selection objects to compute the distance between.
        align_struct_universe (MDAnalsysi.Universe): The MDAnalsysis universe object of the reference structure to
           align the system to.
        align_sel_string (str): A MDAnalysis selection string to use for the structure alignment.
        fstart (int): Optional, the first frame to include in the analysis. Default: 0 (or the first frame)
        fend (int): Optional, the last frame to include in the analysis. Default: -1 (or the last frame)
        fstep (int): Optional, the interval between frames in the analysis when looping from fstart to fend.
            Default: 1 (or every frame)
        axis (str): Optional, the 1d axis to compute the distance in. Default: 'z' (or the z axis)

    Returns:
        (np.array), (list): Returns two outputs. The first is an Numpy array with the timeseries simulation times
            corresponding to the frames in the analysis. The second is list of Numpy arrays with the distances; the
            order in the list corresponds to the atom selection pairs in the mda_selection_pairs input.
    """

    lat_ind = [0, 1]
    dir_ind = 2
    if axis is 'x':
        dir_ind = 0
        lat_ind = [1, 2]
    elif axis is 'y':
        dir_ind = 1
        lat_ind = [0, 2]
    #indices = mda_selection.indices

    nframes = len(universe.trajectory)
    # adjust the end point for slicing
    if fend != -1:
        fend += 1
    if fend == (nframes - 1) and fstart == (nframes - 1):
        fend += 1
    if fend == fstart:
        fend += 1
    if fstart < 0:
        fstart += nframes
    if fend < 0:
        fend += nframes + 1
    times = []
    pair_dists = []
    for pair in mda_selection_pairs:
        pair_dists.append([])
    for frame in universe.trajectory[fstart:fend:fstep]:
        times.append(frame.time)
        # now do the alignment
        align.alignto(universe, align_struct_universe, select=align_sel_string, mass_weighted=True)
        i = 0
        for pair in mda_selection_pairs:
            sel_1 = pair[0]
            sel_2 = pair[1]
            com_1 = sel_1.atoms.center_of_mass()
            com_2 = sel_2.atoms.center_of_mass()
            norm_val_1 = com_1[dir_ind]
            norm_val_2 = com_2[dir_ind]
            dist = np.abs(norm_val_2 - norm_val_1)
            pair_dists[i].append(dist)
            i+=1
    times = np.array(times)
    i=0
    for vals in pair_dists:
        pair_dists[i] = np.array(vals)
        i+=1
    return times, pair_dists

