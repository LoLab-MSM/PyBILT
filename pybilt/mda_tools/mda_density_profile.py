import numpy as np
import pybilt.common.gaussian as gauss

from MDAnalysis.analysis import align

#dictionary of elements and their valence electron counts - used for electron profile density
valence_dict = {'H':1,'C':4,'N':5,'O':6,'P':5}
#dictionary of elements and their atomic numbers - used for electron profile density
Z_dict = {'H':1,'C':6,'N':7,'O':8,'P':15}
#dictionary of elements and their LJ sizes -- CHARMM 36 forcefield
LJ_dict = {'H':1.25,'C':2.00,'N':1.85,'O':1.70,'P':2.15}

#Preston B. Moore, Carlos F. Lopez, Michael L. Klein, Dynamical Properties of a Hydrated Lipid Bilayer
#                from a Multinanosecond Molecular Dynamics Simulation, Biophysical Journal, Volume 81, Issue 5, 2001,
#                Pages 2484-2494, ISSN 0006-3495, http://dx.doi.org/10.1016/S0006-3495(01)75894-8.
#                (http://www.sciencedirect.com/science/article/pii/S0006349501758948)
def electron_density_profile(trajectory,mda_selection, fstart=0,fend=-1,fstep=1, axis='z',nbins=100,reference=0.0,refsel=None,valence=True):
    lat_ind = [0,1]
    dir_ind = 2
    ec_dict = valence_dict
    if not valence:
        ec_dict = Z_dict
    if    axis is 'x':
        dir_ind=0
        lat_ind=[1,2]
    elif axis is 'y':
        dir_ind=1
        lat_ind=[0,2]
    indices = mda_selection.indices
    natoms = len(mda_selection)
    #build the charge array
    charges = np.zeros(natoms)
    a=0
    for atom in mda_selection:
        element = atom.name[0]
        electrons = ec_dict[element]
        partial = atom.charge
        total = electrons-partial
        charges[a]=total
        a+=1

    nframes = len(trajectory)
    #adjust the end point for slicing
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

    #get the maximum box dimension along axis
    bzm=0.0
    nframes = 0
    sel_z = []
    for frame in trajectory[fstart:fend:fstep]:
        bzc = frame.dimensions[dir_ind]
        if bzc>bzm:
            bzm=bzc
        ref_sel_z = 0.0
        if refsel is not None:
            ref_com = refsel.center_of_mass()
            ref_sel_z = ref_com[dir_ind]
            sel_z.append(-ref_sel_z)
        nframes+=1
    shiftzmax = 0.0
    shiftzmin = 0.0
    if refsel is not None:
        #reference=sel_z_avg/nframes
        shiftzmax = min(sel_z)
        shiftzmin = max(sel_z)
    #build the profile axis
    minz = 0.0+shiftzmin
    maxz = bzm+shiftzmax
    #print "minz ",minz," maxz ",maxz
    edges = np.linspace(minz,maxz,(nbins+1),endpoint=True)
    incr = edges[1]-edges[0]
    incr_h = incr/2.0
    centers = np.zeros(nbins)
    nedges = len(edges)
    for i in xrange(1,nedges):
        j=i-1
        centers[j]=edges[j]+incr_h

    counts = np.zeros(nbins)
    #print sel_z
    if refsel is None:
        sel_z = np.zeros(nframes)
    else:
        sel_z = np.array(sel_z)
    f=0
    for frame in trajectory[fstart:fend:fstep]:

        bx = frame.dimensions[lat_ind[0]]
        by = frame.dimensions[lat_ind[1]]
        binvolume = incr*bx*by
        counts_f = np.zeros(nbins)
        sel_pos = frame.positions[indices]
        zpos = sel_pos[:,dir_ind]
        sel_z_curr = sel_z[f]
        zpos+= sel_z_curr
        push_index = (zpos-minz)/incr
        j=0
        for i in push_index:
            ii = int(np.floor(i))
            if ii >= nbins:
                ii = nbins-1
            elif ii < 0:
                ii = 0
            counts_f[ii]+=charges[j]
            j+=1
        counts_f/=binvolume
        counts+=counts_f
        f+=1
    counts/=nframes
    centers-=reference
    return (centers, counts)

#assigns the charge of each atom as a gaussian along the profile direction
#Preston B. Moore, Carlos F. Lopez, Michael L. Klein, Dynamical Properties of a Hydrated Lipid Bilayer
#                from a Multinanosecond Molecular Dynamics Simulation, Biophysical Journal, Volume 81, Issue 5, 2001,
#                Pages 2484-2494, ISSN 0006-3495, http://dx.doi.org/10.1016/S0006-3495(01)75894-8.
#                (http://www.sciencedirect.com/science/article/pii/S0006349501758948)
def electron_density_profile_gaussians(trajectory,mda_selection, fstart=0,fend=-1,fstep=1, axis='z',nbins=100,reference=0.0,refsel=None,valence=True,size_to_sigma=2.0):
    lat_ind = [0,1]
    dir_ind = 2
    ec_dict = valence_dict
    if not valence:
        ec_dict = Z_dict
    if    axis is 'x':
        dir_ind=0
        lat_ind=[1,2]
    elif axis is 'y':
        dir_ind=1
        lat_ind=[0,2]
    indices = mda_selection.indices
    natoms = len(mda_selection)
    #build the charge array
    charges = np.zeros(natoms)
    sizes = np.zeros(natoms)
    a=0
    for atom in mda_selection:
        element = atom.name[0]
        electrons = ec_dict[element]
        partial = atom.charge
        total = electrons-partial
        charges[a]=total
        sizes[a] = LJ_dict[element]
        a+=1

    nframes = len(trajectory)
    #adjust the end point for slicing
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

    #get the maximum box dimension along axis
    bzm=0.0
    nframes = 0
    sel_z = []
    for frame in trajectory[fstart:fend:fstep]:
        bzc = frame.dimensions[dir_ind]
        if bzc>bzm:
            bzm=bzc
        ref_sel_z = 0.0
        if refsel is not None:
            ref_com = refsel.center_of_mass()
            ref_sel_z = ref_com[dir_ind]
            sel_z.append(-ref_sel_z)
        nframes+=1
    shiftzmax = 0.0
    shiftzmin = 0.0
    if refsel is not None:
        #reference=sel_z_avg/nframes
        shiftzmax = min(sel_z)
        shiftzmin = max(sel_z)
    #build the profile axis
    minz = 0.0+shiftzmin
    maxz = bzm+shiftzmax
    #print "minz ",minz," maxz ",maxz
    edges = np.linspace(minz,maxz,(nbins+1),endpoint=True)
    incr = edges[1]-edges[0]
    incr_h = incr/2.0
    centers = np.zeros(nbins)
    nedges = len(edges)
    for i in xrange(1,nedges):
        j=i-1
        centers[j]=edges[j]+incr_h

    counts = np.zeros(nbins)
    #print sel_z
    if refsel is None:
        sel_z = np.zeros(nframes)
    else:
        sel_z = np.array(sel_z)
    f=0
    for frame in trajectory[fstart:fend:fstep]:
        print "doing frame ",f
        bx = frame.dimensions[lat_ind[0]]
        by = frame.dimensions[lat_ind[1]]
        binvolume = incr*bx*by
        counts_f = np.zeros(nbins)
        sel_pos = frame.positions[indices]
        zpos = sel_pos[:,dir_ind]
        sel_z_curr = sel_z[f]
        zpos+= sel_z_curr
        #loop over atoms - build gaussian and get the charge it contributes to each profile bin
        j=0
        for z in zpos:
            sigma = sizes[j]/size_to_sigma
            gc = gauss.GaussianRange((minz,maxz),z,sigma,npoints=(2*nbins))
            #now loop over the edges of the profile axis
            for i in xrange(1,nbins+1):
                zl = edges[i-1]
                zu = edges[i]

                counts_f[i-1]+= (charges[j]*gc.integrate_range(zl,zu))
            j+=1


        counts_f/=binvolume
        counts+=counts_f
        f+=1
    counts/=nframes
    centers-=reference
    return (centers, counts)


def mass_density_profile(trajectory,mda_selection, fstart=0,fend=-1,fstep=1, axis='z',nbins=100,reference=0.0,refsel=None):
    lat_ind = [0,1]
    dir_ind = 2
    if    axis is 'x':
        dir_ind=0
        lat_ind=[1,2]
    elif axis is 'y':
        dir_ind=1
        lat_ind=[0,2]
    indices = mda_selection.indices
    natoms = len(mda_selection)
    #build the mass array
    masses = np.zeros(natoms)
    a=0
    for atom in mda_selection:
        masses[a]=atom.mass
        a+=1

    nframes = len(trajectory)
    #adjust the end point for slicing
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

    #get the maximum box dimension along axis
    bzm=0.0
    nframes = 0
    sel_z = []
    for frame in trajectory[fstart:fend:fstep]:
        bzc = frame.dimensions[dir_ind]
        if bzc>bzm:
            bzm=bzc
        ref_sel_z = 0.0
        if refsel is not None:
            ref_com = refsel.center_of_mass()
            ref_sel_z = ref_com[dir_ind]
            sel_z.append(-ref_sel_z)
        nframes+=1
    shiftzmax = 0.0
    shiftzmin = 0.0
    if refsel is not None:
        #reference=sel_z_avg/nframes
        shiftzmax = min(sel_z)
        shiftzmin = max(sel_z)
    #build the profile axis
    minz = 0.0+shiftzmin
    maxz = bzm+shiftzmax
    edges = np.linspace(minz,maxz,(nbins+1),endpoint=True)
    incr = edges[1]-edges[0]
    incr_h = incr/2.0
    centers = np.zeros(nbins)
    nedges = len(edges)
    for i in xrange(1,nedges):
        j=i-1
        centers[j]=edges[j]+incr_h

    counts = np.zeros(nbins)
    if refsel is None:
        sel_z = np.zeros(nframes)
    else:
        sel_z = np.array(sel_z)
    f=0
    for frame in trajectory[fstart:fend:fstep]:

        bx = frame.dimensions[lat_ind[0]]
        by = frame.dimensions[lat_ind[1]]
        binvolume = incr*bx*by
        counts_f = np.zeros(nbins)
        sel_pos = frame.positions[indices]
        zpos = sel_pos[:,dir_ind]
        sel_z_curr = sel_z[f]
        zpos+= sel_z_curr
        push_index = (zpos-minz)/incr
        j=0
        for i in push_index:
            ii = int(np.floor(i))
            if ii >= nbins:
                ii = nbins-1
            elif ii < 0:
                ii = 0
            counts_f[ii]+=masses[j]
            j+=1
        counts_f/=binvolume
        counts+=counts_f
        f+=1
    counts/=nframes
    centers-=reference
    return (centers, counts)

class SizeError(Exception):
    def __init__(self):
        self.value = "Sizes of input arrays are not the same!"
    def __str__(self):
        return repr(self.value)


def get_intersections(x, y1, y2):

    nx = len(x)
    ny1 = len(y1)
    ny2 = len(y2)
    if nx != ny1 or nx != ny2 or ny1 != ny2:
        raise SizeError()
    yh = np.zeros(nx)
    yl = np.zeros(nx)
    if y1[0]<y2[0]:
        yl = y1
        yh = y2
    else:
        yl = y2
        yh = y1
    output = []
    for i in xrange(1,nx):
        ylc = yl[i]
        ylp = yl[i-1]
        yhp = yh[i-1]
        yhc = yh[i]
        if ylc > yhc and ylp < yhp:
            #get the intersection
            xc = (x[i]+x[i-1])*0.50
            yc = (ylc+yhp+yl[i-1]+yh[i])*0.250
            #now flip the low and high arrays
            temp = yl
            yl = yh
            yh = temp

            output.append([xc,yc])

    return output


def mass_density_profile_multi_align(universe, mda_selections, align_struct_universe, align_sel_string, fstart=0, fend=-1,
                                     fstep=1, axis='z', nbins=100, reference=0.0, refsel=None):
    lat_ind = [0, 1]
    dir_ind = 2
    if axis is 'x':
        dir_ind = 0
        lat_ind = [1, 2]
    elif axis is 'y':
        dir_ind = 1
        lat_ind = [0, 2]
    #indices = mda_selection.indices

    # build the mass array
    masses = {}
    for key in mda_selections.keys():
        natoms = len(mda_selections[key])
        masses[key] = np.zeros(natoms)
        a = 0
        for atom in mda_selections[key]:
            masses[key][a] = atom.mass
            a  += 1

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

    # get the maximum box dimension along axis
    bzm = 0.0
    nframes = 0
    sel_z = []
    system_sel = universe.select_atoms('all')
    for frame in universe.trajectory[fstart:fend:fstep]:
        #get the unaligned system com and axis coordinate
        system_com = system_sel.atoms.center_of_mass()
        system_z = system_com[dir_ind]
        #now do the alignment and get new com and axis coordinates
        align.alignto(universe, align_struct_universe, select=align_sel_string, mass_weighted=True)
        system_com_a = system_sel.atoms.center_of_mass()
        system_z_a = system_com_a[dir_ind]
        dz_a = system_z_a - system_z
        bzc = frame.dimensions[dir_ind]
        if bzc > bzm:
            bzm = bzc
        ref_sel_z = 0.0
        if refsel is not None:
            ref_com = refsel.center_of_mass()
            ref_sel_z = ref_com[dir_ind]
            sel_z.append(-ref_sel_z)
        nframes += 1
    shiftzmax = 0.0
    shiftzmin = 0.0
    if refsel is not None:
        # reference=sel_z_avg/nframes
        shiftzmax = min(sel_z)
        shiftzmin = max(sel_z)
    # build the profile axis
    minz = 0.0 + dz_a + shiftzmin - 0.15*(bzm)
    maxz = bzm + dz_a + shiftzmax + 0.15*(bzm)
    edges = np.linspace(minz, maxz, (nbins + 1), endpoint=True)
    incr = edges[1] - edges[0]
    incr_h = incr / 2.0
    centers = np.zeros(nbins)
    nedges = len(edges)
    for i in xrange(1, nedges):
        j = i - 1
        centers[j] = edges[j] + incr_h

    if refsel is None:
        sel_z = np.zeros(nframes)
    else:
        sel_z = np.array(sel_z)
    f = 0
    out_counts = {}
    for key in mda_selections.keys():
        out_counts[key] = np.zeros(nbins)
    for frame in universe.trajectory[fstart:fend:fstep]:

        bx = frame.dimensions[lat_ind[0]]
        by = frame.dimensions[lat_ind[1]]
        # now do the alignment
        align.alignto(universe, align_struct_universe, select=align_sel_string, mass_weighted=True)
        binvolume = incr * bx * by
        for key in mda_selections.keys():
            indices = mda_selections[key].indices
            counts_f = np.zeros(nbins)
            sel_pos = frame.positions[indices]
            zpos = sel_pos[:, dir_ind]
            sel_z_curr = sel_z[f]
            zpos += sel_z_curr
            push_index = (zpos - minz) / incr
            j = 0
            for i in push_index:
                ii = int(np.floor(i))
                if ii >= nbins:
                    ii = nbins - 1
                elif ii < 0:
                    ii = 0
                counts_f[ii] += masses[key][j]
                j += 1
            counts_f /= binvolume
            out_counts[key] += counts_f
        f += 1
    nframes = float(nframes)
    for key in mda_selections.keys():

        out_counts[key] /= nframes
    centers -= reference

    return centers, out_counts
