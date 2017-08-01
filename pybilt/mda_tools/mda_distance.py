import numpy as np
from MDAnalysis.analysis import align

def com_com_distance_axis_multi_align(universe, mda_selection_pairs, align_struct_universe, align_sel_string, fstart=0, fend=-1,
                                     fstep=1, axis='z'):
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
