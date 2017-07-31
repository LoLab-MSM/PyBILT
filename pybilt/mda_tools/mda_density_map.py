import numpy as np
from MDAnalysis.analysis import align


def position_density_map_2d_multi_align(universe, mda_selections, align_struct_universe, align_sel_string, fstart=0, fend=-1,
                                     fstep=1, norm_axis='z', nbins=100, reference=(0.0, 0.0), refsel=None, normalize=False, scale_to_max=False):
    lat_ind = [0, 1]
    dir_ind = 2
    if norm_axis is 'x':
        dir_ind = 0
        lat_ind = [1, 2]
    elif norm_axis is 'y':
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

    # get the maximum box dimensions in the lateral plane
    bxm = 0.0
    bym = 0.0
    nframes = 0
    sel_x = []
    sel_y = []
    system_sel = universe.select_atoms('all')
    for frame in universe.trajectory[fstart:fend:fstep]:
        #get the unaligned system com and axis coordinate
        system_com = system_sel.atoms.center_of_mass()
        system_x = system_com[lat_ind][0]
        system_y = system_com[lat_ind][1]
        #now do the alignment and get new com and axis coordinates
        align.alignto(universe, align_struct_universe, select=align_sel_string, mass_weighted=True)
        system_com_a = system_sel.atoms.center_of_mass()
        system_x_a = system_com[lat_ind][0]
        system_y_a = system_com[lat_ind][1]
        dx_a = system_x_a - system_x
        dy_a = system_y_a - system_y
        bxc = frame.dimensions[lat_ind][0]
        byc = frame.dimensions[lat_ind][0]

        if bxc > bxm:
            bxm = bxc
        if byc > bym:
            bym = byc
        ref_sel_x = 0.0
        ref_sel_y = 0.0
        if refsel is not None:
            ref_com = refsel.atoms.center_of_mass()
            ref_sel_x = ref_com[lat_ind][0]
            ref_sel_y = ref_com[lat_ind][1]
            sel_x.append(-ref_sel_x)
            sel_y.append(-ref_sel_y)
        nframes += 1
    shiftxmax = 0.0
    shiftxmin = 0.0
    shiftymax = 0.0
    shiftymin = 0.0
    if refsel is not None:
        # reference=sel_z_avg/nframes
        shiftxmax = min(sel_x)
        shiftxmin = max(sel_x)
        shiftymax = max(sel_y)
        shiftymin = min(sel_y)
    # build the 2d density map axes
    minx = 0.0 + dx_a + shiftxmin - 0.15*(bxm)
    maxx = bxm + dx_a + shiftxmax + 0.15*(bxm)
    miny = 0.0 + dy_a + shiftymin - 0.15*(bym)
    maxy = bym + dy_a + shiftymax + 0.15*(bym)
    x_edges = np.linspace(minx, maxx, (nbins + 1), endpoint=True)
    x_incr = x_edges[1] - x_edges[0]
    x_incr_h = x_incr / 2.0
    y_edges = np.linspace(miny, maxy, (nbins + 1), endpoint=True)
    y_incr = y_edges[1] - y_edges[0]
    y_incr_h = y_incr / 2.0

    x_centers = np.zeros(nbins)
    x_nedges = len(x_edges)
    for i in xrange(1, x_nedges):
        j = i - 1
        x_centers[j] = x_edges[j] + x_incr_h
    y_centers = np.zeros(nbins)
    y_nedges = len(x_edges)
    for i in xrange(1, y_nedges):
        j = i - 1
        y_centers[j] = y_edges[j] + y_incr_h

    #counts = np.zeros(nbins)
    if refsel is None:
        sel_x = np.zeros(nframes)
        sel_y = np.zeros(nframes)
    else:
        sel_x = np.array(sel_x)
        sel_y = np.array(sel_y)
    f = 0
    out_counts = {}
    for key in mda_selections.keys():
        out_counts[key] = np.zeros((nbins,nbins))
    for frame in universe.trajectory[fstart:fend:fstep]:

        # now do the alignment
        align.alignto(universe, align_struct_universe, select=align_sel_string, mass_weighted=True)

        for key in mda_selections.keys():
            indices = mda_selections[key].indices
            counts_f = np.zeros((nbins, nbins))
            sel_pos = frame._pos[indices]
            xpos = sel_pos[:, lat_ind[0]]
            ypos = sel_pos[:, lat_ind[1]]
            sel_x_curr = sel_x[f]
            sel_y_curr = sel_y[f]
            xpos += sel_x_curr
            ypos += sel_y_curr
            x_push_index = (xpos - minx) / x_incr
            y_push_index = (ypos - miny) / y_incr
            j = 0
            for i in range(len(x_push_index)):
                ii = int(np.floor(x_push_index[i]))
                jj = int(np.floor(y_push_index[i]))
                if ii >= nbins:
                    ii = nbins - 1
                elif ii < 0:
                    ii = 0
                if jj >= nbins:
                    jj = nbins -1
                elif jj < 0:
                    jj = 0
                counts_f[ii][jj] += 1.0
                j += 1
            #counts_f /= binvolume
            out_counts[key] += counts_f
        f += 1
    nframes = float(nframes)
    #print(out_counts)
    if normalize:
        for key in mda_selections.keys():
            out_counts[key] /= out_counts[key].sum()
            #out_counts[key] /= np.max(out_counts[key])
    elif scale_to_max:
        for key in mda_selections.keys():
            #out_counts[key] /= out_counts[key].sum()
            out_counts[key] /= np.max(out_counts[key])
    x_centers -= reference[0]
    y_centers -= reference[1]

    return x_centers, y_centers, out_counts

