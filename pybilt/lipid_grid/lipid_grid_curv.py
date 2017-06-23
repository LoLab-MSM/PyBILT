'''
    Classes and functions to implement gridding and curvature correlation analysis for lipid bilayers.
    The gridding and anlaysis procedures are based on
    the decription given in section "Correlation between bilayer surface curvature and the
clustering of lipid molecules" of Koldso H, Shorthouse D, He lie J, Sansom MSP (2014) Lipid Clustering Correlates with Membrane Curvature as Revealed by Molecular Simulations of
Complex Lipid Bilayers. PLoS Comput Biol 10(10): e1003911. doi:10.1371/journal.pcbi.1003911
However, this implementation currently uses the z position (or normal position) of the lipids' centers of mass, while
their implementaion uses "the z coordinate of the interface between the head groups of the
lipids (excluding the current species being calculated and tails in
that box."

'''
import numpy as np
#import my running stats class
from pybilt.common.running_stats import *

class LipidGrid_2d:
    def __init__(self, com_frame, com_frame_indices,plane,nxbins=20,nybins=20):
        #store the frame and leaflet
        self.frame = com_frame
        #self.leaflet = ms_leaflet
        #get the x and y indices
        ix = plane[0]
        iy = plane[1]
        iz = [i for i in [0,1,2] if i not in plane][0]
        #get the box dimemsions
        box = com_frame.box
        boxx = box[ix]
        boxy = box[iy]
        box_com = com_frame.mem_com
        box_com_x = box_com[ix]
        box_com_y = box_com[iy]
        #save the numbers of bins
        self.x_nbins = nxbins
        self.y_nbins = nybins
        #initialize the edges of the and centers of the gridpoints
        # x
        #self.x_min = -box_com_x
        #self.x_max = boxx - box_com_x
        self.x_min = 0.0
        self.x_max = boxx
        self.x_edges = np.linspace(self.x_min,self.x_max,(nxbins+1),endpoint=True)
        self.x_incr = self.x_edges[1]-self.x_edges[0]
        x_incr_h = self.x_incr/2.0
        self.x_centers = np.zeros(nxbins)
        self.x_nedges = len(self.x_edges)
        for i in xrange(1,self.x_nedges):
            j=i-1
            self.x_centers[j]=self.x_edges[j]+x_incr_h

        # y
        #self.y_min = -box_com_y
        #self.y_max = boxy - box_com_y
        self.y_min = 0.0
        self.y_max = boxy
        self.y_edges = np.linspace(self.y_min,self.y_max,(nybins+1),endpoint=True)
        self.y_incr = self.y_edges[1]-self.y_edges[0]
        y_incr_h = self.y_incr/2.0
        self.y_centers = np.zeros(nybins)
        self.y_nedges = len(self.y_edges)
        for i in xrange(1,self.y_nedges):
            j=i-1
            self.y_centers[j]=self.y_edges[j]+y_incr_h
        self.x_length = self.x_max-self.x_min
        self.y_length = self.y_max-self.y_min
        # get the lipid indices for this leaflet
        indices = com_frame_indices
        void_ind = max(indices)+1
        #now assign lipids to the gridpoints

        self.lipid_grid = []
        #cx = 0
        #print self.x_edges
        mx_x = -1000.0
        mn_x = 1000.0
        for cx in range(len(self.x_edges)-1):
            self.lipid_grid.append([])
            x_lower = self.x_edges[cx]
            x_upper = self.x_edges[cx+1]
            #print "x_upper ",x_upper, " x_lower ",x_lower
            for cy in range(len(self.y_edges)-1):
                self.lipid_grid[cx].append([])
                y_lower = self.y_edges[cy]
                y_upper = self.y_edges[cy+1]
                #check lipid COMs
                for i in indices:
                    xi = com_frame.lipidcom[i].com[ix]
                    yi = com_frame.lipidcom[i].com[iy]
                    zi = com_frame.lipidcom[i].com_unwrap[iz]
                    x_box = xi > x_lower and xi < x_upper
                    y_box = yi > y_lower and yi < y_upper
                   # print "x_lower: ", x_lower, " xi: ", xi, " x_ upper: ", x_upper, " x_box: ", x_box
                   #quit()
                    if xi < mn_x:
                        mn_x = xi
                    if xi > mx_x:
                        mx_x = xi
                    if x_box and y_box:

                       # print "y_lower: ", y_lower, " yi: ", yi, " y_upper: ", y_upper, " y_box: ",y_box
                       # print
                        #add to this grid
                        self.lipid_grid[cx][cy].append((i, com_frame.lipidcom[i].type, zi))
                        #print "lipid index ",i," of type ",com_frame.lipidcom[i].type, " added to grid (",cx," ",cy,")"
        #print "minimum x coord: ", mn_x
        #print "maximum x coord: ", mx_x

    def get_index_at(self,ix,iy):
        return self.lipid_grid[ix][iy][:,0]

    def get_z_at(self,ix,iy):
        return self.lipid_grid[ix][iy][:,2]



class LipidGrids:
    def __init__(self, com_frame, leaflets,plane,nxbins=3,nybins=3):
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
        self.leaf_grid['upper'] = LipidGrid_2d(com_frame,upper_indices,plane,nxbins=nxbins,nybins=nybins)
        #lower
        lower_indices = leaflets['lower'].get_member_indices()
        self.leaf_grid['lower'] = LipidGrid_2d(com_frame,lower_indices,plane,nxbins=nxbins,nybins=nybins)
        return


    def norm_displacement_cross_correlation(self):
        output = dict()
        for leaf in self.leaflets.keys():
            output[leaf] = dict()
            ll_types = self.leaflets[leaf].get_group_names()
           # print "ll_types: ", ll_types
            for l_type in ll_types:
                #loop over grid boxes
                count = []
                z_vals = []
                n_box = 0.0
                #print (self.leaf_grid[leaf].lipid_grid)
               # print(len(self.leaf_grid[leaf].lipid_grid))
                for xb in self.leaf_grid[leaf].lipid_grid:
                    #print(len(xb))
                    #print "xb: "
                    #print xb
                    for yb in xb:
                        #print "yb :"
                        #print yb
                        #print(len(yb))
                        box_count = 0
                        box_z_vals = []
                        for lipid in yb:
                         #   print(lipid)
                            lipid_type = lipid[1]
                            lipid_z = lipid[2]
                            #print "lipid_z: ",lipid_z
                            if lipid_type == l_type:
                                box_count+=1
                            else:
                                box_z_vals.append(lipid_z)
                        #if len(yb) > 0:
                        #    n_box+=1.0
                        n_box+=1
                        if len(box_z_vals) > 0:
                            #n_box+=1.0
                            box_z_avg = box_z_vals[0]
                            if box_z_vals > 1:
                                box_z_avg = np.array(box_z_vals).mean()
                            #print " box_z_vals: "
                            #print box_z_vals
                           # print "box_z_avg: ",box_z_avg
                            #print "yb: "
                            #print yb
                            count.append(float(box_count))
                            z_vals.append(box_z_avg)
                cross_corr = 0.0
                if len(count) > 1 and len(z_vals) >1:
                   # print "l_type: ",l_type
                    count = np.array(count)
                   # print "count: ", count
                    z_vals = np.array(z_vals)
                   # print "z_vals: ", z_vals
                    count_mean = count.mean()
                    count_std = count.std()
                    z_mean = z_vals.mean()
                    z_std = z_vals.std()
                    cross_sum = np.dot(count-count_mean, z_vals-z_mean)
                    cross_corr = cross_sum/(count_std*z_std*n_box)
                   # print "n_box: ",n_box

                if np.isnan(cross_corr):
                   # print"cross_corr: ",cross_corr," cross_sum ", cross_sum," n_box ",n_box, " z_mean ",z_mean, " z_std ", z_std
                   # print "count:"
                   # print count
                   # print "z_vals:"
                   # print z_vals
                   # quit()
                    cross_corr = 0.0
                output[leaf][l_type] = cross_corr
        #quit()
        return output
