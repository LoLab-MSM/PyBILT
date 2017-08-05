'''
    Classes and functions to implement lipid COM gridding and analysis for lipid bilayers. Acts on MemSys objects.
    The gridding and anlaysis procedures are based on
    the decriptions given in Gapsys et al. J Comput Aided Mol Des (2013) 27:845-858,
    which is itself a modified version of the GridMAT-MD method by Allen et al. Vol. 30, No. 12 Journal of Computational Chemistry.
    However, I have currently left out several bits of the extra functionality, e.g. the handling of an embedded protein.
'''

import numpy as np
# import my running stats class
from pybilt.common.running_stats import *


class LipidGrid_2d(object):
    def __init__(self, mda_frame, mda_universe, mda_frame_resids, plane, nxbins=50, nybins=50, embedded_protein=None):
        # store the frame and leaflet
        self.frame = mda_frame
        # self.leaflet = ms_leaflet
        # get the x and y indices
        ix = plane[0]
        iy = plane[1]
        iz = [i for i in [0, 1, 2] if i not in plane][0]
        # get the box dimemsions
        box = mda_frame.box[plane]
        boxx = box[ix]
        boxy = box[iy]
        # save the numbers of bins
        self.x_nbins = nxbins
        self.y_nbins = nybins
        # initialize the edges of the and centers of the gridpoints
        # x
        self.x_min = 0.0
        self.x_max = boxx
        self.x_edges = np.linspace(self.x_min, self.x_max, (nxbins + 1), endpoint=True)
        self.x_incr = self.x_edges[1] - self.x_edges[0]
        x_incr_h = self.x_incr / 2.0
        self.x_centers = np.zeros(nxbins)
        self.x_nedges = len(self.x_edges)
        for i in xrange(1, self.x_nedges):
            j = i - 1
            self.x_centers[j] = self.x_edges[j] + x_incr_h

        # y
        self.y_min = 0.0
        self.y_max = boxy
        self.y_edges = np.linspace(self.y_min, self.y_max, (nybins + 1), endpoint=True)
        self.y_incr = self.y_edges[1] - self.y_edges[0]
        y_incr_h = self.y_incr / 2.0
        self.y_centers = np.zeros(nybins)
        self.y_nedges = len(self.y_edges)
        for i in xrange(1, self.x_nedges):
            j = i - 1
            self.y_centers[j] = self.y_edges[j] + y_incr_h
        self.x_length = self.x_max - self.x_min
        self.y_length = self.y_max - self.y_min
        # get the lipid indices for this leaflet
        resids = mda_frame_resids
        #void_ind = max(indices) + 1
        # now assign lipids to the gridpoints
        self.lipid_grid = np.zeros((nxbins, nybins), dtype=np.int)
        self.lipid_grid_z = np.zeros((nxbins, nybins))
        self.lipid_grid_resnames = []
        bxh = boxx / 2.0
        byh = boxy / 2.0
        cx = 0
        for x in self.x_centers:
            self.lipid_grid_resnames.append([])
            cy = 0
            for y in self.y_centers:
                r_min = 1.0e10
                i_min = 0
                z_min = 0.0
                resname_min = "UNK"
                # check lipid COMs
                for i in resids:
                    res_sel_string = "resid "+str(i)
                    res_sel = mda_universe.select_atoms(res_sel_string)
                    resname = res_sel.resname
                    res_indices = res_sel.indices
                    for index in res_indices:
                        pos = mda_frame._pos[index]
                        xi = pos[ix]
                        yi = pos[iy]
                        zi = pos[iz]
                        # print "iz ",iz," zi ",zi
                        dx = x - xi
                        dy = y - yi
                        # Minimum image -- coordinates must be pre-wrapped
                        if np.absolute(dx) > bxh:
                            dx = boxx - np.absolute(x - bxh) - np.absolute(xi - bxh)
                        if np.absolute(dy) > bxh:
                            dy = boxy - np.absolute(y - byh) - np.absolute(yi - byh)
                        rxy = np.sqrt(dx ** 2 + dy ** 2)
                        if rxy < r_min:
                            r_min = rxy
                            i_min = i
                            z_min = zi
                            resname_min = resname
                        # if embedded_protein is not None:

                        # print "i_min ",i_min," z_min ",z_min
                        # if cutoff is not None:

                        # else:
                self.lipid_grid[cx, cy] = i_min
                self.lipid_grid_z[cx, cy] = z_min
                self.lipid_grid_resnames[cx].append(resname_min)
                cy += 1
            cx += 1

    def get_index_at(self, ix, iy):
        return self.lipid_grid[ix, iy]

    def get_z_at(self, ix, iy):
        return self.lipid_grid_z[ix, iy]

    # Outputs the grid as an xyz coordinate file
    def write_xyz(self, xyz_name):
        # Open up the file to write to
        xyz_out = open(xyz_name, "w")
        npoints = self.x_nbins * self.y_nbins
        comment = "Leaflet Grid " + self.leaflet.name
        xyz_out.write(str(npoints))
        xyz_out.write("\n")
        xyz_out.write(comment)
        xyz_out.write("\n")

        cx = 0
        for x in self.x_centers:
            cy = 0
            for y in self.y_centers:
                # get the z coordinate
                z = self.lipid_grid_z[cx, cy]
                # get the lipid resname
                ic = self.lipid_grid[cx, cy]
                oname = self.lipid_grid_resnames[cx][cy]
                # write to file


                line = str(oname) + " " + str(x) + " " + str(y) + " " + str(z)

                xyz_out.write(line)
                xyz_out.write("\n")
                cy += 1
            cx += 1
        xyz_out.close()
        return


class LipidGrids(object):
    def __init__(self, mda_frame, mda_universe, leaflets, plane, nxbins=50, nybins=50, embedded_protein=None):
        # store the frame and leaflet
        self.frame = mda_frame
        self.leaflets = leaflets
        self.plane = plane
        self.norm = [i for i in [0, 1, 2] if i not in plane][0]
        self.nbins_x = nxbins
        self.nbins_y = nybins
        self.leaf_grid = {}
        self.myframe = mda_frame.frame
        # initialize the grids
        # upper
        upper_resids = leaflets['upper'].get_member_resids()
        self.leaf_grid['upper'] = LipidGrid_2d(mda_frame, mda_universe, upper_resids, plane, nxbins=nxbins, nybins=nybins)
        # lower
        lower_resids = leaflets['lower'].get_member_resids()
        self.leaf_grid['lower'] = LipidGrid_2d(mda_frame, mda_universe, lower_resids, plane, nxbins=nxbins, nybins=nybins)
        return

    def thickness_grid(self):
        tgrid = np.zeros((self.nbins_x, self.nbins_y))
        for ix in xrange(self.nbins_x):
            for iy in xrange(self.nbins_y):
                zu = self.leaf_grid['upper'].get_z_at(ix, iy)
                zl = self.leaf_grid['lower'].get_z_at(ix, iy)
                dz = zu - zl
                tgrid[ix, iy] = dz
                if dz < 0.0:
                    print "Warning!!--MD frame number ", self.myframe, " --Value thickness less than zero (", dz, ") at grid point ", ix, " ", iy
        return tgrid

    def average_thickness(self, return_grid=False):
        trun = RunningStats()
        tgrid = self.thickness_grid()
        for ix in xrange(self.nbins_x):
            for iy in xrange(self.nbins_y):
                tc = tgrid[ix, iy]
                trun.push(tc)
        avg_out = (trun.mean(), trun.deviation())
        if return_grid:
            return avg_out, tgrid
        else:
            return avg_out

    def map_to_grid(self, com_values_dict, leaflet='both'):
        do_leaflet = []
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        elif leaflet == "upper" or leaflet == "lower":
            do_leaflet.append(leaflet)
        else:
            # unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'", leaflet, "\' from the LeafletGrids of frame ", self.myframe
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            do_leaflet.append('upper')
            do_leaflet.append('lower')

        out_dict = {}
        for leaf in do_leaflet:
            out_dict[leaf] = np.zeros((self.nbins_x, self.nbins_y))
            for ix in xrange(self.nbins_x):
                for iy in xrange(self.nbins_y):
                    com_ind = self.leaf_grid[leaf].get_index_at(ix, iy)
                    value = com_values_dict[com_ind]
                    out_dict[leaf][ix, iy] = value

        return out_dict

    def area_per_lipid(self):
        do_leaflet = []
        do_leaflet.append('upper')
        do_leaflet.append('lower')
        # get the unique type/resnames in the system
        resnames = []
        for leaf in do_leaflet:
            for group in self.leaflets[leaf].groups:
                gname = group.name()
                if gname not in resnames:
                    resnames.append(gname)
        # initialize counters for each residue/type
        area_run_per_res_type = {}

        for name in resnames:
            area_run_per_res_type[name] = RunningStats()

        area_per_lipid = {}

        area_run = RunningStats()
        for leaf in do_leaflet:
            area_per_bin = self.leaf_grid[leaf].x_incr * self.leaf_grid[leaf].y_incr
            lip_ind = self.leaflets[leaf].get_member_indices()
            for i in lip_ind:
                rname = self.frame.lipidcom[i].type
                locations = np.where(self.leaf_grid[leaf].lipid_grid == i)
                nlocs = len(locations[0])
                # print locations
                # print 'nlocs ',nlocs
                area = area_per_bin * nlocs
                area_per_lipid[i] = area
                area_run_per_res_type[rname].push(area)

                area_run.push(area)

        average_per_res = {}
        for name in resnames:
            average = area_run_per_res_type[name].mean()
            std = area_run_per_res_type[name].deviation()
            average_per_res[name] = (average, std)
        system_average = area_run.mean()
        system_dev = area_run.deviation()

        output = (system_average, average_per_res, area_per_lipid)
        return output

    def curvature(self):
        nxb = self.nbins_x
        nyb = self.nbins_y
        # first order derivtives
        sx_u = np.zeros((nxb, nyb))
        sy_u = np.zeros((nxb, nyb))
        sx_l = np.zeros((nxb, nyb))
        sy_l = np.zeros((nxb, nyb))

        for ix in xrange(nxb):
            for iy in xrange(nyb):
                ixp = ix - 1
                if ixp < 0:
                    ixp += nxb
                ixn = ix + 1
                if ixn >= nxb:
                    ixn -= nxb
                iyp = ix - 1
                if iyp < 0:
                    iyp += nyb
                iyn = iy + 1
                if iyn >= nyb:
                    iyn -= nyb
                    # upper
                    ## using central difference for numerical first derivative
                sx = self.leaf_grid['upper'].lipid_grid_z[ixn, iy] - self.leaf_grid['upper'].lipid_grid_z[ixp, iy]
                sx /= (self.leaf_grid['upper'].x_incr) ** 2
                sy = self.leaf_grid['upper'].lipid_grid_z[ix, iyn] - self.leaf_grid['upper'].lipid_grid_z[ix, iyp]
                sy /= (self.leaf_grid['upper'].y_incr) ** 2
                sx_u[ix, iy] = sx
                sy_u[ix, iy] = sy
                # lower
                sx = self.leaf_grid['lower'].lipid_grid_z[ixn, iy] - self.leaf_grid['lower'].lipid_grid_z[ixp, iy]
                sx /= (self.leaf_grid['lower'].x_incr) ** 2
                sy = self.leaf_grid['lower'].lipid_grid_z[ix, iyn] - self.leaf_grid['lower'].lipid_grid_z[ix, iyp]
                sy /= (self.leaf_grid['lower'].y_incr) ** 2
                sx_l[ix, iy] = sx
                sy_l[ix, iy] = sy
        # now do second order derivatives - central difference numerical derivative of the first derivative
        ssx_u = np.zeros((nxb, nyb))
        ssy_u = np.zeros((nxb, nyb))
        ssxy_u = np.zeros((nxb, nyb))
        ssx_l = np.zeros((nxb, nyb))
        ssy_l = np.zeros((nxb, nyb))
        ssxy_l = np.zeros((nxb, nyb))
        for ix in xrange(nxb):
            for iy in xrange(nyb):
                ixp = ix - 1
                if ixp < 0:
                    ixp += nxb
                ixn = ix + 1
                if ixn >= nxb:
                    ixn -= nxb
                iyp = ix - 1
                if iyp < 0:
                    iyp += nyb
                iyn = iy + 1
                if iyn >= nyb:
                    iyn -= nyb
                    # upper
                    ## using central difference for numerical first derivative
                ssx = sx_u[ixn, iy] - sx_u[ixp, iy]
                ssx /= (self.leaf_grid['upper'].x_incr) ** 2
                ssy = sy_u[ix, iyn] - sy_u[ix, iyp]
                ssy /= (self.leaf_grid['upper'].y_incr) ** 2
                ssxy = sx_u[ix, iyn] - sx_u[ix, iyp]
                ssxy /= (self.leaf_grid['upper'].y_incr) ** 2
                ssx_u[ix, iy] = ssx
                ssy_u[ix, iy] = ssy
                ssxy_u[ix, iy] = ssxy

                # lower
                ssx = sx_l[ixn, iy] - sx_l[ixp, iy]
                ssx /= (self.leaf_grid['lower'].x_incr) ** 2
                ssy = sy_l[ix, iyn] - sy_l[ix, iyp]
                ssy /= (self.leaf_grid['lower'].y_incr) ** 2
                ssxy = sx_l[ix, iyn] - sx_l[ix, iyp]
                ssxy /= (self.leaf_grid['upper'].y_incr) ** 2
                ssx_l[ix, iy] = ssx
                ssy_l[ix, iy] = ssy
                ssxy_l[ix, iy] = ssxy
        # now get curvatures
        curv_mean_u = np.zeros((nxb, nyb))
        curv_gauss_u = np.zeros((nxb, nyb))
        curv_mean_l = np.zeros((nxb, nyb))
        curv_gauss_l = np.zeros((nxb, nyb))
        dx_u = self.leaf_grid['upper'].x_incr
        dy_u = self.leaf_grid['upper'].y_incr
        dx_l = self.leaf_grid['lower'].x_incr
        dy_l = self.leaf_grid['lower'].y_incr
        for ix in xrange(nxb):
            for iy in xrange(nyb):
                # upper
                sx = sx_u[ix, iy]
                sy = sy_u[ix, iy]
                ssx = ssx_u[ix, iy]
                ssy = ssy_u[ix, iy]
                ssxy = ssxy_u[ix, iy]
                sx_v = np.array([self.leaf_grid['upper'].x_centers[ix] + dx_u, 0.0, sx])
                sy_v = np.array([0.0, self.leaf_grid['upper'].y_centers[iy] + dy_u, sy])
                ssx_v = np.array([self.leaf_grid['upper'].x_centers[ix] + dx_u, 0.0, ssx])
                ssy_v = np.array([0.0, self.leaf_grid['upper'].y_centers[iy] + dy_u, ssy])
                ssxy_v = np.array([0.0, self.leaf_grid['upper'].y_centers[iy] + dy_u, ssxy])
                E = np.dot(sx_v, sx_v)
                F = np.dot(sx_v, sy_v)
                G = np.dot(sy_v, sy_v)
                n = np.cross(sx_v, sy_v)
                n /= np.linalg.norm(n)
                L = np.dot(ssx_v, n)
                M = np.dot(ssxy_v, n)
                N = np.dot(ssy_v, n)
                # mean curvature
                J = (E * N + G * L - 2.0 * F * M) / (2.0 * (E * G - F) ** 2)
                # Gaussian curvature
                K = (L * N - M ** 2) / (E * G - F ** 2)
                curv_mean_u[ix, iy] = J
                curv_gauss_u[ix, iy] = K
                # lower
                sx = sx_l[ix, iy]
                sy = sy_l[ix, iy]
                ssx = ssx_l[ix, iy]
                ssy = ssy_l[ix, iy]
                ssxy = ssxy_l[ix, iy]
                sx_v = np.array([self.leaf_grid['lower'].x_centers[ix] + dx_u, 0.0, sx])
                sy_v = np.array([0.0, self.leaf_grid['lower'].y_centers[iy] + dy_u, sy])
                ssx_v = np.array([self.leaf_grid['lower'].x_centers[ix] + dx_u, 0.0, ssx])
                ssy_v = np.array([0.0, self.leaf_grid['lower'].y_centers[iy] + dy_u, ssy])
                ssxy_v = np.array([0.0, self.leaf_grid['lower'].y_centers[iy] + dy_u, ssxy])
                E = np.dot(sx_v, sx_v)
                F = np.dot(sx_v, sy_v)
                G = np.dot(sy_v, sy_v)
                n = np.cross(sx_v, sy_v)
                n /= np.linalg.norm(n)
                L = np.dot(ssx_v, n)
                M = np.dot(ssxy_v, n)
                N = np.dot(ssy_v, n)
                # mean curvature
                J = (E * N + G * L - 2.0 * F * M) / (2.0 * (E * G - F) ** 2)
                # Gaussian curvature
                K = (L * N - M ** 2) / (E * G - F ** 2)
                curv_mean_l[ix, iy] = J
                curv_gauss_l[ix, iy] = K

        return ((curv_mean_u, curv_gauss_u), (curv_mean_l, curv_gauss_l))

    def grid_to_dict(self, in_grid, leaflet='upper'):
        out_dict = {}
        for ix in xrange(self.nbins_x):
            for iy in xrange(self.nbins_y):
                l_i = self.leaf_grid[leaflet].get_index_at(ix, iy)
                grid_val = in_grid[ix, iy]
                out_dict[l_i] = grid_val
        return out_dict

    def get_xyzc(self, leaflet='both', zvalue_dict=None, color_dict=None, color_grid=None, color_type_dict=None):
        do_leaflet = []
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        elif leaflet == "upper" or leaflet == "lower":
            do_leaflet.append(leaflet)
        else:
            # unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'", leaflet, "\' from the LeafletGrids of frame ", self.myframe
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        out_dict = {}
        npoints = (self.nbins_x * self.nbins_y)
        X = np.zeros(npoints)
        Y = np.zeros(npoints)
        Z = np.zeros(npoints)
        C = np.zeros(npoints)
        if color_dict is not None:
            if len(color_dict.shape) == 2:
                C = np.zeros((npoints, color_dict.shape[1]))
        if color_type_dict is not None:
            dict_type = type(color_type_dict[color_type_dict.keys()[0]])
            if dict_type is str:
                C = np.zeros(npoints, dtype=np.str)
        for leaf in do_leaflet:
            npt = 0
            cx = 0
            for x in self.leaf_grid[leaf].x_centers:
                cy = 0
                for y in self.leaf_grid[leaf].y_centers:
                    # get the z coordinate
                    z = self.leaf_grid[leaf].get_z_at(cx, cy)
                    ic = self.leaf_grid[leaf].get_index_at(cx, cy)
                    # optionally pull z value from lipid index dictionary
                    if zvalue_dict is not None:
                        z = zvalue_dict[ic]
                    X[npt] = x
                    Y[npt] = y
                    Z[npt] = z
                    if color_dict is not None:
                        C[npt] = color_dict[ic]
                    if color_grid is not None:
                        C[npt] = color_grid[cx, cy]
                    if color_type_dict is not None:
                        ltype = self.frame.lipidcom[ic].type
                        color_curr = color_type_dict[ltype]
                        C[npt] = color_curr

                    npt += 1
                    cy += 1
                cx += 1
                # if color_dict is not None and len(color_dict.shape)==1:
                #     col_min = min(C)
                #     C-=col_min
                #     col_max = max(C)
                #     C/=col_max

                # elif color_grid is not None:
                #     col_min = min(C)
                #     C-=col_min
                #     col_max = max(C)
                #     C/=col_max

            out_dict[leaf] = (X, Y, Z, C)
        return out_dict

    def write_xyz(self, leaflet='both', zvalue_dict='Default', out_path="./"):
        do_leaflet = []
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        elif leaflet == "upper" or leaflet == "lower":
            do_leaflet.append(leaflet)
        else:
            # unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'", leaflet, "\' from the LeafletGrids of frame ", self.myframe
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        out_name = out_path + "leaflet_grid_f" + str(self.myframe) + "_"
        for leaf in do_leaflet:
            out_name += leaf[0]
        out_name += ".xyz"
        # Open up the file to write to
        xyz_out = open(out_name, "w")
        npoints = (self.nbins_x * self.nbins_y) * len(do_leaflet)
        comment = "Leaflet Grid in xyz coordinate format"
        xyz_out.write(str(npoints))
        xyz_out.write("\n")
        xyz_out.write(comment)
        xyz_out.write("\n")

        for leaf in do_leaflet:
            cx = 0
            for x in self.leaf_grid[leaf].x_centers:
                cy = 0
                for y in self.leaf_grid[leaf].y_centers:
                    # get the z coordinate
                    z = self.leaf_grid[leaf].get_z_at(cx, cy)
                    ic = self.leaf_grid[leaf].get_index_at(cx, cy)
                    # optionally pull z value from lipid index dictionary
                    if zvalue_dict is not 'Default':
                        z = zvalue_dict[ic]
                        # get the lipid resname

                    oname = self.frame.lipidcom[ic].type
                    # write to file


                    line = str(oname) + " " + str(x) + " " + str(y) + " " + str(z)

                    xyz_out.write(line)
                    xyz_out.write("\n")
                    cy += 1
                cx += 1
        xyz_out.close()
        return

    def get_integer_type_arrays(self):
        grids_dict = {}
        # build the integer type id's
        group_to_int = {}
        groups = []
        i = 0
        for leaf in self.leaflets.keys():
            # groups in the leaflet
            lgroups = self.leaflets[leaf].get_group_names()
            ngroups = len(groups)
            # build dictionary for string group name to integer type
            for name in lgroups:
                if name not in groups:
                    groups.append(name)
                    group_to_int[name] = i
                    i += 1
        for leaf in self.leaf_grid.keys():
            nxbins = self.leaf_grid[leaf].x_nbins
            nybins = self.leaf_grid[leaf].y_nbins
            type_array = np.zeros((nxbins, nybins), dtype=np.int)
            cx = 0
            for x in self.leaf_grid[leaf].x_centers:
                cy = 0
                for y in self.leaf_grid[leaf].y_centers:
                    # get the z coordinate
                    ic = self.leaf_grid[leaf].get_index_at(cx, cy)
                    oname = self.frame.lipidcom[ic].type
                    itype = group_to_int[oname]
                    type_array[cx, cy] = itype
                    cy += 1
                cx += 1
            grids_dict[leaf] = np.copy(type_array)

        return grids_dict, group_to_int

    def get_one_array_per_leaflet(self):
        grids_dict = {}
        # build the integer type id's
        group_to_int = {}
        groups = []
        i = 0
        for leaf in self.leaf_grid.keys():
            # groups in the leaflet
            lgroups = self.leaflets[leaf].get_group_names()
            ngroups = len(groups)
            # build dictionary for string group name to integer type
            for name in lgroups:
                if name not in groups:
                    groups.append(name)
                    group_to_int[name] = i
                    i += 1
        for leaf in self.leaf_grid.keys():
            nxbins = self.leaf_grid[leaf].x_nbins
            nybins = self.leaf_grid[leaf].y_nbins
            # type_array = np.zeros((nxbins, nybins, 4))
            type_array = []
            cx = 0
            for x in self.leaf_grid[leaf].x_centers:
                cy = 0
                for y in self.leaf_grid[leaf].y_centers:
                    # get the z coordinate
                    ic = self.leaf_grid[leaf].get_index_at(cx, cy)
                    oname = self.frame.lipidcom[ic].type
                    itype = group_to_int[oname]
                    # type_array[cx,cy,0] = x
                    # type_array[cx,cy,1] = y
                    # type_array[cx,cy,2] = self.leaflets[leaf].lipid_grid[cx,cy]
                    # type_array[cx,cy,3] = itype
                    type_array.append((x, y, self.leaf_grid[leaf].lipid_grid[cx, cy], itype, oname))
                    cy += 1
                cx += 1
            type_array = np.array(type_array, dtype='f8, f8, i4, i4, |S10')
            grids_dict[leaf] = np.copy(type_array)
        return grids_dict
