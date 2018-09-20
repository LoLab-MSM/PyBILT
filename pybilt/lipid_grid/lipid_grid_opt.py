"""Build lipid grids derived from COMFrame objects.
Classes and functions to implement lipid COM gridding and analysis for
lipid bilayers. This module defines version that build grids off of
COMFrame objects and is meant primarily for internal use by the
BilayerAnalyzer class. The gridding and anlaysis procedures are based on
the descriptions given in Gapsys et al. J Comput Aided Mol Des (2013)
27:845-858, which is itself a modified version of the GridMAT-MD method
by Allen et al. Vol. 30, No. 12 Journal of Computational Chemistry.
However, I have currently left out bits of the extra functionality like
the handling of an embedded proteins.
"""

# TODO (blakeaw1102@gmail): Add embedded protein functionality to lipid grid.

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from builtins import object
import numpy as np
import sys
from scipy.ndimage.filters import gaussian_filter
# pybilt imports
from pybilt.common.running_stats import RunningStats
from pybilt.common.distance_cutoff_clustering import distance_euclidean_pbc
from six.moves import range




def grid_curvature(x_vals, y_vals, zgrid):
    """Compute the Mean and Gaussian curvature across a grid.
    Args:
        x_vals (np.array): The bin labels along the x-axis of the gridded data.
        y_vals (np.arrray): The bin labels along the y-axis of the gridded data.
        zgrid (np.array): The 2d grid of bin values.

    Returns:
        tuple: Returns a 2 item tuple with the 2d numpy arrays of the curvatures with
            format (mean curvature, Gaussian curvature).
    """
    nxb = len(x_vals)
    nyb = len(y_vals)
    x_incr = x_vals[1]-x_vals[0]
    y_incr = y_vals[1]-y_vals[0]
    # print("x_incr {} y_incr {}".format(x_incr, y_incr))
    [sy, sx] = np.gradient(zgrid, y_incr, x_incr)
    [syy, syx] = np.gradient(sy, y_incr, x_incr)
    [sxy, sxx] = np.gradient(sx, y_incr, x_incr)
    #now get curvatures
    curv_mean_u = np.zeros((nxb, nyb))
    curv_gauss_u = np.zeros((nxb, nyb))
    for ix in range(nxb):
        for iy in range(nyb):
            #upper
            sx_c = sx[ix, iy]
            sy_c = sy[ix, iy]
            ssx = sxx[ix, iy]
            ssy = syy[ix, iy]
            ssxy = sxy[ix, iy]
            sx_v = np.array([x_incr, 0.0, sx_c])
            sy_v = np.array([0.0, y_incr, sy_c])
            ssx_v = np.array([x_incr, 0.0, ssx])
            ssy_v = np.array([0.0, y_incr, ssy])
            ssxy_v = np.array([0.0, y_incr, ssxy])
            E = np.dot(sx_v, sx_v)
            F = np.dot(sx_v, sy_v)
            G = np.dot(sy_v, sy_v)
            n = np.cross(sx_v, sy_v)
            n /=np.linalg.norm(n)
            L = np.dot(ssx_v, n)
            M = np.dot(ssxy_v, n)
            N = np.dot(ssy_v, n)
            #mean curvature
            J = (E*N+G*L-2.0*F*M)/(2.0*(E*G-F)**2)
            #Gaussian curvature
            K = (L*N-M**2)/(E*G-F**2)
            curv_mean_u[ix, iy] = J
            curv_gauss_u[ix, iy] = K
            # print("ix: {} iy: {} J: {} K: {}".format(ix,iy,J,K))


    return (curv_mean_u, curv_gauss_u)

def grid_surface_area(x_vals, y_vals, zgrid):
    """Compute the surface area across a regular 2d grid.
    Args:
        x_vals (np.array): The bin labels along the x-axis of the gridded data.
        y_vals (np.arrray): The bin labels along the y-axis of the gridded data.
        zgrid (np.array): The 2d grid of bin values.

    Returns:
        float: Returns a the surface area estimate.
    """
    nxb = len(x_vals)
    nyb = len(y_vals)
    x_incr = x_vals[1]-x_vals[0]
    y_incr = y_vals[1]-y_vals[0]
    # print("x_incr {} y_incr {}".format(x_incr, y_incr))
    [sy, sx] = np.gradient(zgrid, y_incr, x_incr)
    #now get curvatures
    sa = 0.0
    for ix in range(nxb):
        for iy in range(nyb):
            #upper
            sx_c = sx[ix, iy]
            sy_c = sy[ix, iy]
            sx_v = np.array([1.0, 0.0, sx_c])
            sy_v = np.array([0.0, 1.0, sy_c])
            cross = np.cross(sx_v, sy_v)
            dA = np.sqrt(np.dot(cross, cross))*x_incr*y_incr
            sa += dA
    return sa

class LipidGrid2d(object):
    """A 2d lipid grid object.

    This object is used by the LipidGrids object to construct a 2d grid for
    a bilayer leaflet and assign lipids to it using the coordinates derived
    from a COMFrame representation object.

    Attributes:
        frame (COMFrame): Stores a local copy of the COMFrame object from
            which the the coordinate data and lipid type data is derived from.
        x_nbins (int): The number of grid bins in the 'x' dimension.
        y_nbins (int): The number of grid bins in the 'y' dimension.
        x_min (float): The lower boundary of the grid range in the 'x'
            dimension.
        x_max (float): The upper boundary of the grid range in the 'x'
            dimension.
        x_incr (float): The size of grid boxes (or spacing between grid
            points) in the 'x' dimension.
        x_centers (np.array): The center points of the grid boxes in the 'x'
            dimension.
        x_edges (np.array): The edges of the grid boxes in the 'x' dimension.
        x_centers (np.array): The centers of the grid boxes in the 'x'
            dimension.
        y_min (float): The lower boundary of the grid range in the 'y'
            dimension.
        y_max (float): The upper boundary of the grid range in the 'y'
            dimension.
        y_incr (float): The size of grid boxes (or spacing between grid
            points) in the 'y' dimension.
        y_centers (np.array): The center points of the grid boxes in the 'y'
            dimension.
        y_edges (np.array): The edges of the grid boxes in the 'y' dimension.
        y_centers (np.array): The centers of the grid boxes in the 'y'
            dimension.
        lipid_grid (np.array): A 2d array of size x_nbins*y_nbins that stores
            the index of lipids from the COMFrame that are assigned to each
            grid box.
        lipid_grid_z (np.array): A 2d array of size x_nbins*y_nbins that stores
            the z coordinate of the lipids assigned to each grid box.

    """
    def __init__(self, com_frame, com_frame_indices, plane, nxbins=50,
                 nybins=50):
        """Initialize the LipidGrid2d object.
        This version of the LipidGrid2d object uses a k-neareset neighbors
        (knn) based algorithm to assign lipids to grid points. The method runs
        knn over the lipid com positions first. Then the first grid point (0,0)
        is assigned by normal distance search. Then subsequent grid points are
        assigned by distance search over nearest neighbor list of the lipid
        assigned to the previous corner grid point
        (i.e. [max(0,i-1), max(0,j-1)] ). The overall speed improvement is
        roughly 2-3x with a system of 300 lipids per leaflet.

        Args:
            com_frame (COMFrame): The instance of COMFrame from which to
                pull the coordinates for lipids to use when building the grid.
            com_frame_indices (list): A list COMFrame lipid indices to
                include when building the grid.
            plane (list): The indices from the 3d coordinates for the
                coordinates that correspond to the bilayer lateral plane.
            nxbins (Optional[int]): The number of bins along the 'x'
                dimension, i.e. along the dimension corresponding to plane[0].
                Defaults to 50.
            nybins (Optional[int): The number of bins along the 'y'
                dimension, i.e. along the dimension corresponding to
                plane[1]. Defaults to 50.
        """
        # store the frame and leaflet
        self.frame = com_frame
        # self.leaflet = ms_leaflet
        # get the x and y indices
        ix = plane[0]
        iy = plane[1]
        iz = [i for i in [0, 1, 2] if i not in plane][0]
        # get the box dimemsions
        box = com_frame.box[plane]
        boxx = box[ix]
        boxy = box[iy]
        # save the numbers of bins
        self.x_nbins = nxbins
        self.y_nbins = nybins
        # initialize the edges of the and centers of the gridpoints
        # x
        self.x_min = 0.0
        self.x_max = boxx
        self.x_edges = np.linspace(self.x_min, self.x_max, (nxbins + 1),
                                   endpoint=True)
        self.x_incr = self.x_edges[1] - self.x_edges[0]
        x_incr_h = self.x_incr / 2.0
        self.x_centers = np.zeros(nxbins)
        self._x_nedges = len(self.x_edges)
        for i in range(1, self._x_nedges):
            j = i-1
            self.x_centers[j] = self.x_edges[j] + x_incr_h

        # y
        self.y_min = 0.0
        self.y_max = boxy
        self.y_edges = np.linspace(self.y_min, self.y_max, (nybins + 1),
                                   endpoint=True)
        self.y_incr = self.y_edges[1] - self.y_edges[0]
        y_incr_h = self.y_incr / 2.0
        self.y_centers = np.zeros(nybins)
        self.y_nedges = len(self.y_edges)
        for i in range(1, self._x_nedges):
            j = i - 1
            self.y_centers[j] = self.y_edges[j] + y_incr_h
        self.x_length = self.x_max - self.x_min
        self.y_length = self.y_max - self.y_min
        # get the lipid indices for this leaflet
        indices = com_frame_indices
        # now assign lipids to the gridpoints
        self.lipid_grid = np.zeros((nxbins, nybins), dtype=np.int)
        self.lipid_grid_z = np.zeros((nxbins, nybins))
        bxh = boxx / 2.0
        byh = boxy / 2.0
        knn = self._knn(com_frame, com_frame_indices, plane)
        # assign the first point
        cx = 0
        for x in self.x_centers:
            cy = 0
            for y in self.y_centers:
                r_min = 1.0e10
                i_min = 0
                z_min = 0.0
                #check lipid COMs
                if cx == 0 and cy == 0:
                    for i in indices:
                        xi = com_frame.lipidcom[i].com[ix]
                        yi = com_frame.lipidcom[i].com[iy]
                        zi = com_frame.lipidcom[i].com_unwrap[iz]
                        #print "iz ",iz," zi ",zi
                        dx = x - xi
                        dy = y - yi
                        #Minimum image -- coordinates must be pre-wrapped
                        if np.absolute(dx) > bxh:
                            dx = boxx - np.absolute(x - bxh) - np.absolute(xi -
                                                                           bxh)
                        if np.absolute(dy) > bxh:
                            dy = boxy - np.absolute(y - byh) - np.absolute(yi -
                                                                           byh)
                        rxy = np.sqrt(dx**2 + dy**2)
                        if rxy < r_min:
                            r_min = rxy
                            i_min = i
                            z_min = zi

                    self.lipid_grid[cx,cy] = i_min
                    self.lipid_grid_z[cx,cy] = z_min
                else:
                    nn_ix = max([0, cx-1])
                    nn_iy = max([0, cy-1])
                    nn_ind = self.lipid_grid[nn_ix, nn_iy]
                    r_min = 1.0e10
                    i_min = nn_ind
                    z_min = self.lipid_grid[nn_ix, nn_iy]
                    #check neighboring grids
                    neighbor_cells = []
                    #print(nn_ix, cx)
                    for xx in range(cx-1, cx):
                        xx = max([0, xx])
                        for yy in range(cy-1, cy):
                            yy = max([0, yy])
                            iii = self.lipid_grid[xx, yy]
                            #print(xx, yy, iii)
                            if iii not in neighbor_cells:
                                neighbor_cells.append(iii)
                    # print(neighbor_cells)
                    # if cx > 1:
                    #     quit()
                    for ind in neighbor_cells:

                        xi = com_frame.lipidcom[ind].com[ix]
                        yi = com_frame.lipidcom[ind].com[iy]
                        zi = com_frame.lipidcom[ind].com_unwrap[iz]
                        #print "iz ",iz," zi ",zi
                        dx = x - xi
                        dy = y - yi
                        #Minimum image -- coordinates must be pre-wrapped
                        if np.absolute(dx) > bxh:
                            dx = boxx - np.absolute(x - bxh) - np.absolute(xi -
                                                                           bxh)
                        if np.absolute(dy) > bxh:
                            dy = boxy - np.absolute(y - byh) - np.absolute(yi -
                                                                           byh)
                        rxy = np.sqrt(dx**2 + dy**2)
                        if rxy < r_min:
                            r_min = rxy
                            i_min = ind
                            z_min = zi
                    for neighbors in knn[i_min]:
                        n_ind = neighbors[0]
                        xi = com_frame.lipidcom[n_ind].com[ix]
                        yi = com_frame.lipidcom[n_ind].com[iy]
                        zi = com_frame.lipidcom[n_ind].com_unwrap[iz]
                        #print "iz ",iz," zi ",zi
                        dx = x - xi
                        dy = y - yi
                        #Minimum image -- coordinates must be pre-wrapped
                        if np.absolute(dx) > bxh:
                            dx = boxx - np.absolute(x - bxh) - np.absolute(xi -
                                                                           bxh)
                        if np.absolute(dy) > bxh:
                            dy = boxy - np.absolute(y - byh) - np.absolute(yi -
                                                                           byh)
                        rxy = np.sqrt(dx**2 + dy**2)
                        if rxy < r_min:
                            r_min = rxy
                            i_min = n_ind
                            z_min = zi
                    self.lipid_grid[cx, cy] = i_min
                    self.lipid_grid_z[cx,cy] = z_min
                cy += 1
            cx += 1

    def get_index_at(self, ix, iy):
        """Returns the COMFrame index of the lipid at the specified position
        in the lipid_grid.

        Args:
            ix (int): The 'x' index in the lipid_grid.
            iy (): The 'y' index in the lipid_grid.

        Returns:
            int: The index of the lipid.

        """
        return self.lipid_grid[ix, iy]

    def get_z_at(self, ix, iy):
        """Returns the z coordinate of the lipid at the specified position
        in the lipid_grid.

        Args:
            ix (int): The 'x' index in the lipid_grid.
            iy (): The 'y' index in the lipid_grid.

        Returns:
            float: The z coordinate of the lipid.

        """
        return self.lipid_grid_z[ix, iy]

    def z_perturb_grid(self):
        """Returns the array with z coordinates shifted by the mean.

        Returns:
            np.array: The mean shifted z coordinate array.

        """
        z_grid = self.lipid_grid_z
        z_avg = z_grid.mean()
        z_pert = z_grid - z_avg
        return z_pert


    # Outputs the grid as an xyz coordinate file
    def write_xyz(self, xyz_name):
        """Write out the lipid grid as an xyz coordinate file.

        Args:
            xyz_name (str): File path and name for the output file.


        """
        # Open up the file to write to
        xyz_out = open(xyz_name, "w")
        npoints = self.x_nbins*self.y_nbins
        comment = "Leaflet Grid " + self.leaflet.name
        xyz_out.write(str(npoints))
        xyz_out.write("\n")
        xyz_out.write(comment)
        xyz_out.write("\n")

        cx=0
        for x in self.x_centers:
            cy=0
            for y in self.y_centers:
                #get the z coordinate
                z = self.lipid_grid_z[cx,cy]
                #get the lipid resname
                ic = self.lipid_grid[cx,cy]
                oname = self.frame.lipidcom[ic].type
                #write to file


                line = str(oname)+" "+str(x)+" "+str(y)+" "+str(z)

                xyz_out.write(line)
                xyz_out.write("\n")
                cy+=1
            cx+=1
        xyz_out.close()
        return

    @staticmethod
    def _knn(com_frame, com_frame_indices, plane):
        # get the x and y indices
        ix = plane[0]
        iy = plane[1]
        iz = [i for i in [0, 1, 2] if i not in plane][0]
        # get the box dimemsions
        box = com_frame.box[plane]
        boxx = box[ix]
        boxy = box[iy]
        l_box = np.array([boxx, boxy])
        n_ind = len(com_frame_indices)
        # initialize knn dict
        knn = {key: [] for key in com_frame_indices}
        nn_k = 24
        for i in range(n_ind-1):
            i_i = com_frame_indices[i]
            pos_a = com_frame.lipidcom[i_i].com[[ix, iy]]
            for j in range(i+1, n_ind):
                j_j = com_frame_indices[j]
                pos_b = com_frame.lipidcom[j_j].com[[ix, iy]]
                dist = distance_euclidean_pbc(pos_a, pos_b, l_box,
                                              center='box_half')
                #ltype = ba_reps['com_frame'].lipidcom[j].type
                #print "ltype: ",ltype," dist ",dist
                knn[i_i].append([j_j, dist])
                knn[j_j].append([i_i, dist])
        for key in knn.keys():
            neighbors = knn[key]
            neighbors.sort(key=lambda x: x[1])
            nn_neighbors = neighbors[:nn_k]
            knn[key] = nn_neighbors

        return knn

class LipidGrids(object):
    def __init__(self, com_frame, leaflets, plane, nxbins=50, nybins=50):
        #store the frame and leaflet
        self.frame = com_frame
        self.leaflets = leaflets
        self.plane = plane
        self.norm = [i for i in [0,1,2] if i not in plane][0]
        self.nbins_x = nxbins
        self.nbins_y = nybins
        self.leaf_grid = {}
        self.myframe = com_frame.mdnumber
        #initialize the grids
        #upper
        upper_indices = leaflets['upper'].get_member_indices()
        self.leaf_grid['upper'] = LipidGrid2d(com_frame, upper_indices,
                                              plane, nxbins=nxbins,
                                              nybins=nybins)
        #lower
        lower_indices = leaflets['lower'].get_member_indices()
        self.leaf_grid['lower'] = LipidGrid2d(com_frame, lower_indices,
                                              plane, nxbins=nxbins,
                                              nybins=nybins)
        return

    def thickness_grid(self):
        tgrid = np.zeros((self.nbins_x, self.nbins_y))
        tgrid = self.leaf_grid['upper'].lipid_grid_z - self.leaf_grid['lower'].lipid_grid_z
        # for ix in range(self.nbins_x):
        #     for iy in range(self.nbins_y):
        #         zu = self.leaf_grid['upper'].get_z_at(ix, iy)
        #         zl = self.leaf_grid['lower'].get_z_at(ix, iy)
        #         dz = zu - zl
        #         tgrid[ix,iy] = dz
        #         if dz < 0.0:
        #             print "Warning!!--MD frame number ",self.myframe," --Value thickness less than zero (",dz,") at grid point ",ix," ",iy
        return tgrid

    def average_thickness(self, return_grid=False):
        trun = RunningStats()
        tgrid = self.thickness_grid()
        avg_out = (tgrid.mean(), tgrid.std())
        # for ix in range(self.nbins_x):
        #     for iy in range(self.nbins_y):
        #         tc = tgrid[ix, iy]
        #         trun.push(tc)
        # avg_out = (trun.mean(), trun.deviation())
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
            #unknown option--use default "both"
            print("!! Warning - request for unknown leaflet name \'",leaflet,"\' from the LeafletGrids of frame ",self.myframe)
            print("!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\"")
            do_leaflet.append('upper')
            do_leaflet.append('lower')

        out_dict = {}
        for leaf in do_leaflet:
            out_dict[leaf] = np.zeros((self.nbins_x,self.nbins_y))
            for ix in range(self.nbins_x):
                for iy in range(self.nbins_y):
                    com_ind=self.leaf_grid[leaf].get_index_at(ix,iy)
                    value = com_values_dict[com_ind]
                    out_dict[leaf][ix,iy]=value

        return out_dict

    def area_per_lipid(self):
        do_leaflet = []
        do_leaflet.append('upper')
        do_leaflet.append('lower')
        #get the unique type/resnames in the system
        resnames = []
        for leaf in do_leaflet:
            for group in self.leaflets[leaf].groups:
                gname = group.name()
                if gname not in resnames:
                    resnames.append(gname)
        #initialize counters for each residue/type
        area_run_per_res_type = {}

        for name in resnames:
            area_run_per_res_type[name]=RunningStats()

        area_per_lipid = {}

        area_run = RunningStats()
        for leaf in do_leaflet:
            area_per_bin = self.leaf_grid[leaf].x_incr*self.leaf_grid[leaf].y_incr
            lip_ind = self.leaflets[leaf].get_member_indices()
            for i in lip_ind:
                rname = self.frame.lipidcom[i].type
                locations = np.where(self.leaf_grid[leaf].lipid_grid == i)
                nlocs = len(locations[0])
                #print locations
                #print 'nlocs ',nlocs
                area = area_per_bin*nlocs
                area_per_lipid[i]=area
                area_run_per_res_type[rname].push(area)

                area_run.push(area)

        average_per_res = {}
        for name in resnames:
            average = area_run_per_res_type[name].mean()
            std = area_run_per_res_type[name].deviation()
            average_per_res[name] = (average,std)
        system_average = area_run.mean()
        #system_dev = area_run.deviation()

        output = (system_average, average_per_res, area_per_lipid)
        return output

    def curvature(self, use_gaussian_filter=True, filter_sigma=10.0, filter_mode='nearest'):
        x_vals = self.leaf_grid['upper'].x_centers
        y_vals = self.leaf_grid['upper'].y_centers
        z_grid = self.leaf_grid['upper'].lipid_grid_z
        if use_gaussian_filter:
            z_grid = gaussian_filter(z_grid, filter_sigma, mode=filter_mode)
        curv_upper = grid_curvature(x_vals, y_vals, z_grid)
        x_vals = self.leaf_grid['lower'].x_centers
        y_vals = self.leaf_grid['lower'].y_centers
        z_grid = self.leaf_grid['lower'].lipid_grid_z
        if use_gaussian_filter:
            z_grid = gaussian_filter(z_grid, filter_sigma, mode=filter_mode)
        curv_lower = grid_curvature(x_vals, y_vals, z_grid)

        return (curv_upper, curv_lower)

    def surface_area(self, use_gaussian_filter=True, filter_sigma=10.0, filter_mode='nearest'):
        x_vals = self.leaf_grid['upper'].x_centers
        y_vals = self.leaf_grid['upper'].y_centers
        z_grid = self.leaf_grid['upper'].lipid_grid_z
        if use_gaussian_filter:
            z_grid = gaussian_filter(z_grid, filter_sigma, mode=filter_mode)
        sa_upper = grid_surface_area(x_vals, y_vals, z_grid)
        x_vals = self.leaf_grid['lower'].x_centers
        y_vals = self.leaf_grid['lower'].y_centers
        z_grid = self.leaf_grid['lower'].lipid_grid_z
        if use_gaussian_filter:
            z_grid = gaussian_filter(z_grid, filter_sigma, mode=filter_mode)
        sa_lower = grid_surface_area(x_vals, y_vals, z_grid)
        return (sa_upper, sa_lower)

    def grid_to_dict(self,in_grid,leaflet='upper'):
        out_dict = {}
        for ix in range(self.nbins_x):
            for iy in range(self.nbins_y):
                l_i = self.leaf_grid[leaflet].get_index_at(ix,iy)
                grid_val = in_grid[ix,iy]
                out_dict[l_i]=grid_val
        return out_dict

    def get_xyzc(self,leaflet='both',zvalue_dict=None,color_dict=None,color_grid=None, color_type_dict=None):
        do_leaflet = []
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        elif leaflet == "upper" or leaflet == "lower":
            do_leaflet.append(leaflet)
        else:
            #unknown option--use default "both"
            print("!! Warning - request for unknown leaflet name \'",leaflet,"\' from the LeafletGrids of frame ",self.myframe)
            print("!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\"")
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        out_dict = {}
        npoints = (self.nbins_x*self.nbins_y)
        X = np.zeros(npoints)
        Y = np.zeros(npoints)
        Z = np.zeros(npoints)
        C = np.zeros(npoints)
        if color_dict is not None:
            if len(color_dict.shape)==2:
                C = np.zeros((npoints, color_dict.shape[1]))
        if color_type_dict is not None:
            #dict_type = type(color_type_dict[color_type_dict.keys()[0]])
            C = list()
        for leaf in do_leaflet:
            npt = 0
            cx=0
            for x in self.leaf_grid[leaf].x_centers:
                cy=0
                for y in self.leaf_grid[leaf].y_centers:
                    #get the z coordinate
                    z = self.leaf_grid[leaf].get_z_at(cx,cy)
                    ic = self.leaf_grid[leaf].get_index_at(cx,cy)
                    #optionally pull z value from lipid index dictionary
                    if zvalue_dict is not None:
                        z = zvalue_dict[ic]
                    X[npt]=x
                    Y[npt]=y
                    Z[npt]=z
                    if color_dict is not None:
                        C[npt]=color_dict[ic]
                    if color_grid is not None:
                        C[npt]=color_grid[cx,cy]
                    if color_type_dict is not None:
                        ltype = self.frame.lipidcom[ic].type
                        color_curr = color_type_dict[ltype]
                        C.append(color_curr)

                    npt+=1
                    cy+=1
                cx+=1

            if color_type_dict is not None:
                C = np.array(C)
            out_dict[leaf]=(X,Y,Z,C)
        return out_dict

    def write_xyz(self,leaflet='both',zvalue_dict='Default',out_path="./"):
        do_leaflet = []
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        elif leaflet == "upper" or leaflet == "lower":
            do_leaflet.append(leaflet)
        else:
            #unknown option--use default "both"
            print("!! Warning - request for unknown leaflet name \'",leaflet,"\' from the LeafletGrids of frame ",self.myframe)
            print("!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\"")
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        out_name = out_path+"leaflet_grid_f"+str(self.myframe)+"_"
        for leaf in do_leaflet:
            out_name+=leaf[0]
        out_name+=".xyz"
        # Open up the file to write to
        xyz_out = open(out_name, "w")
        npoints = (self.nbins_x*self.nbins_y)*len(do_leaflet)
        comment = "Leaflet Grid in xyz coordinate format"
        xyz_out.write(str(npoints))
        xyz_out.write("\n")
        xyz_out.write(comment)
        xyz_out.write("\n")

        for leaf in do_leaflet:
            cx=0
            for x in self.leaf_grid[leaf].x_centers:
                cy=0
                for y in self.leaf_grid[leaf].y_centers:
                    #get the z coordinate
                    z = self.leaf_grid[leaf].get_z_at(cx,cy)
                    ic = self.leaf_grid[leaf].get_index_at(cx,cy)
                    #optionally pull z value from lipid index dictionary
                    if zvalue_dict is not 'Default':
                        z = zvalue_dict[ic]
                    #get the lipid resname

                    oname = self.frame.lipidcom[ic].type
                    #write to file


                    line = str(oname)+" "+str(x)+" "+str(y)+" "+str(z)

                    xyz_out.write(line)
                    xyz_out.write("\n")
                    cy+=1
                cx+=1
        xyz_out.close()
        return

    def write_pdb(self, pdb_name, leaflet='both', z_grid_upper=None,
                  z_grid_lower=None, beta_grid_upper=None,
                  beta_grid_lower=None,
                  use_gaussian_filter=False, filter_sigma=10.0,
                  filter_mode='nearest'):
        """Write out the lipid grid as an PDB coordinate file.

        Args:
            pdb_name (str): File path and name for the output file.
            leaflet (Optional[str]): Specify which leaflets to write to the PDB
                file. The options are 'both', 'upper', or 'lower'. Defaults to
                'both'.
            z_grid_upper (Optional[np.array]): A 2d grid of values corresponding
                to the elements of the upper leaflet of the lipid grid that are
                to be written as the z-coordinate in the PDB file for upper
                leaflet members. Defaults to None.
            z_grid_lower (Optional[np.array]): A 2d grid of values corresponding
                to the elements of the lower leaflet of the lipid grid that are
                to be written as the z-coordinate in the PDB file for lower
                leaflet members. Defaults to None.
            beta_grid_upper (Optional[np.array]): A 2d grid of values corresponding
                to the elements of the upper leaflet of the lipid grid that are
                to be written in the Beta column of the PDB file for upper
                leaflet members. Defaults to None.
            beta_grid_lower (Optional[np.array]): A 2d grid of values corresponding
                to the elements of the lower leaflet of the lipid grid that are
                to be written in the Beta column of the PDB file for upper
                leaflet members. Defaults to None.
            use_gaussian_filter (Optional[bool]): Use SciPy's Gaussian filter
                to filter the z-coordinates before outputting the PDB file.
                Defaults to False. This option overrides inputs for
                z_grid_upper and z_grid_lower.
            filter_sigma (Optional[bool]): Set the sigma value for the Gaussian
                filter. Defaults to 10.0
            filter_mode (Optional['str']): Set the mode for Gaussian filter.
                Defaults to 'nearest'

        """

        do_leaflet = []
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        elif leaflet == "upper" or leaflet == "lower":
            do_leaflet.append(leaflet)
        else:
            #unknown option--use default "both"
            print("!! Warning - request for unknown leaflet name \'",leaflet,"\' from the LeafletGrids of frame ",self.myframe)
            print("!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\"")
            do_leaflet.append('upper')
            do_leaflet.append('lower')

        # Open up the file to write to
        pdb_out = open(pdb_name, "w")
        npoints = (self.nbins_x*self.nbins_y)*len(do_leaflet)
        # First line to write is the box dimensions
        box = self.frame.box
        box_data = "CRYST1  {:06.3f}  {:06.3f}  {:06.3f}".format(box[0],
                                                                 box[1],
                                                                 box[2])
        box_data += "  90.000  90.000 90.00 P 1           1"
        pdb_out.write(box_data)
        pdb_out.write("\n")
        if use_gaussian_filter:
            z_grid_upper = self.leaf_grid['upper'].lipid_grid_z
            z_max = z_grid_upper.max()
            z_grid_upper = gaussian_filter(z_grid_upper, filter_sigma, mode=filter_mode)
            z_max_f = z_grid_upper.max()
            z_diff = z_max_f - z_max
            # print("upper: {}".format(z_diff))
            z_grid_upper -= z_diff
            z_grid_lower = self.leaf_grid['lower'].lipid_grid_z
            z_min = z_grid_lower.min()
            z_grid_lower = gaussian_filter(z_grid_lower, filter_sigma, mode=filter_mode)
            z_min_f = z_grid_lower.min()
            z_diff = z_min_f - z_min
            # print("lower: {}".format(z_diff))
            z_grid_lower -= z_diff
        index = 1
        for leaf in do_leaflet:
            cx=0
            for x in self.leaf_grid[leaf].x_centers:
                cy=0
                for y in self.leaf_grid[leaf].y_centers:
                    ic = self.leaf_grid[leaf].get_index_at(cx, cy)
                    #get the z coordinate
                    z = self.leaf_grid[leaf].get_z_at(cx,cy)
                    if (leaf == 'upper') and (z_grid_upper is not None):
                        z = z_grid_upper[cx, cy]
                    elif (leaf == 'lower') and (z_grid_lower is not None):
                        z = z_grid_lower[cx, cy]
                    #get the lipid resname
                    oname = self.frame.lipidcom[ic].type
                    # Compose the elements of the line
                    # ATOM, columns 1-4, char
                    ATOM = "ATOM"
                    # Atom serial number, columns 7-11, right justified, int
                    serial = "{:5d}".format(index)
                    # Atom name, columns 13-16, left justified, char
                    name = "{:4.4}".format(oname)
                    # Alternate location indicator, column 17, char
                    alternate = "{:1.1}".format(" ")
                    # Residue name, columns 18-20, right justified, char
                    resname = "{:>3.3}".format(oname)
                    # Chain identifier, column 22, char
                    chain = "{:1.1}".format(leaf)
                    # Residue sequence number, columns 23-26, right justified, int
                    rsn = "{:>4d}".format(index)
                    if len(rsn) > 4: rsn = rsn[:(4-len(rsn))]
                    # Code for insertions of residues, column 27, char
                    cir = " "
                    # X orthogonal Angstrom coordinate, columns 31-38, right justified, real (8.3)
                    x_pos = "{:8.3f}".format(x)
                    # Y orthogonal Angstrom coordinate, columns 39-46, right justified, real (8.3)
                    y_pos = "{:8.3f}".format(y)
                    # Z orthogonal Angstrom coordinate, columns 31-38, right justified, real (8.3)
                    z_pos = "{:8.3f}".format(z)
                    # Occupancy, columns 55-60, right justified, real (6.2)
                    occupancy = "{:6.2f}".format(1.0)
                    # Temperature factor (Beta), columns 61-66, right justified, real (6.2)
                    beta = "{:6.2f}".format(0.0)
                    if (leaf == 'upper') and (beta_grid_upper is not None):
                        b_val = beta_grid_upper[cx, cy]
                        beta = "{:6.2f}".format(b_val)
                    elif (leaf == 'lower') and (beta_grid_lower is not None):
                        b_val = beta_grid_lower[cx, cy]
                        beta = "{:6.2f}".format(b_val)
                    # Segment identifier, columns 73-76, left justified, char
                    segid = "{:4.3}".format(leaf)
                    # Element symbol, columns 77-78, right justified, char
                    element = "{:>2.1}".format(oname)
                    # Charge, columns 79-80, char
                    charge = "{:02.1f}".format(0.0)
                    # Now put them all together to build the line
                    #      1-4     5-6    7-11     12   13-16    17
                    line = ATOM + "  " + serial + " " + name + alternate
                    #        18-20    19     22    23-26  27    28-30   31-38   39-46
                    line += resname + " " + chain + rsn + cir + "   " + x_pos + y_pos
                    #       47-54    55-60      61-66    67-72    73-76
                    line += z_pos + occupancy + beta + "      " + segid
                    #        77-78     79-80
                    line += element + charge
                    pdb_out.write(line)
                    pdb_out.write("\n")
                    index += 1
                    cy+=1
                cx+=1
        pdb_out.write('END')
        pdb_out.close()

        return

    def get_integer_type_arrays(self):
        grids_dict = {}
        # build the integer type id's
        group_to_int = {}
        groups = []
        i = 0
        for leaf in self.leaflets.keys():
            #groups in the leaflet
            lgroups = self.leaflets[leaf].get_group_names()
            # build dictionary for string group name to integer type
            for name in lgroups:
                if name not in groups:
                    groups.append(name)
                    group_to_int[name] = i
                    i+=1
        for leaf in self.leaf_grid.keys():
            nxbins = self.leaf_grid[leaf].x_nbins
            nybins = self.leaf_grid[leaf].y_nbins
            type_array = np.zeros((nxbins, nybins), dtype=np.int)
            cx=0
            for dummy_x in self.leaf_grid[leaf].x_centers:
                cy=0
                for dummy_y in self.leaf_grid[leaf].y_centers:
                    #get the z coordinate
                    ic = self.leaf_grid[leaf].get_index_at(cx,cy)
                    oname = self.frame.lipidcom[ic].type
                    itype = group_to_int[oname]
                    type_array[cx,cy] = itype
                    cy+=1
                cx+=1
            grids_dict[leaf] = np.copy(type_array)

        return grids_dict, group_to_int

    def get_one_array_per_leaflet(self):
        grids_dict = {}
        # build the integer type id's
        group_to_int = {}
        groups = []
        i = 0
        for leaf in self.leaf_grid.keys():
            #groups in the leaflet
            lgroups = self.leaflets[leaf].get_group_names()
            # build dictionary for string group name to integer type
            for name in lgroups:
                if name not in groups:
                    groups.append(name)
                    group_to_int[name] = i
                    i+=1
        for leaf in self.leaf_grid.keys():
            #type_array = np.zeros((nxbins, nybins, 4))
            type_array = []
            cx=0
            for x in self.leaf_grid[leaf].x_centers:
                cy=0
                for y in self.leaf_grid[leaf].y_centers:
                    #get the z coordinate
                    ic = self.leaf_grid[leaf].get_index_at(cx,cy)
                    oname = self.frame.lipidcom[ic].type
                    itype = group_to_int[oname]
                    #type_array[cx,cy,0] = x
                    #type_array[cx,cy,1] = y
                    #type_array[cx,cy,2] = self.leaflets[leaf].lipid_grid[cx,cy]
                    #type_array[cx,cy,3] = itype
                    type_array.append((x, y, self.leaf_grid[leaf].lipid_grid[cx,cy], itype, oname))
                    cy+=1
                cx+=1
            type_array = np.array(type_array, dtype='f8, f8, i4, i4, |S10')
            grids_dict[leaf] = np.copy(type_array)
        return grids_dict
