from pybilt.bilayer_analyzer import BilayerAnalyzer
import numpy as np
from pybilt.plot_generation import plot_generation_functions as pgf
from pybilt.plot_generation.plot_generation_functions import _color_list
from scipy.ndimage.filters import gaussian_filter
from six.moves import range
def test_lipid_grid_curvature():
    sel_string = "resname POPC DOPE TLCL2"
    name_dict = {'DOPE':['P'],'POPC':['P'],'TLCL2':['P1','P3']}
    analyzer = BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="resname POPC DOPE TLCL2")
    analyzer.remove_analysis('msd_1')
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    analyzer._analysis_protocol.use_objects['lipid_grid'] = True
    # adjust the number of bins for gridding
    nbins = 75
    analyzer.rep_settings['lipid_grid']['n_xbins'] = nbins
    analyzer.rep_settings['lipid_grid']['n_ybins'] = nbins
    analyzer.set_frame_range(0,1,1)
    i =0
    zgrid = None
    for _frame in analyzer:
        if i == 0:
            zgrid = analyzer.reps['lipid_grid'].leaf_grid['upper'].lipid_grid_z
        else:
            zgrid += analyzer.reps['lipid_grid'].leaf_grid['upper'].lipid_grid_z
        i+=1
    zgrid_i = zgrid/float(i)
    zgrid_i_f = gaussian_filter(zgrid_i, 5.0, mode="nearest")
    x_vals = analyzer.reps['lipid_grid'].leaf_grid['upper'].x_centers
    y_vals = analyzer.reps['lipid_grid'].leaf_grid['upper'].y_centers
    curvature_grids_i = grid_curvature(x_vals,
                                    y_vals,
                                    zgrid_i_f)
    xyzc_mean = analyzer.reps['lipid_grid'].get_xyzc(leaflet='upper',
                                                     color_grid=curvature_grids_i[0])['upper']
    #xyzc_u_mean_i = (x_vals, y_vals, x_vals, curvature_grids_i[0].flatten())
    pgf.plot_grid_as_scatter(xyzc_mean, save=False, show=False, colorbar=True)
    analyzer.reps['lipid_grid'].curvature()

    return


def grid_curvature(x_vals, y_vals, zgrid):
    nxb = len(x_vals)
    nyb = len(y_vals)
    x_incr = x_vals[1]-x_vals[0]
    y_incr = y_vals[1]-y_vals[0]
    # print("x_incr {} y_incr {}".format(x_incr, y_incr))
    [sy, sx] = np.gradient(zgrid, y_incr, x_incr)
    [syy, syx] = np.gradient(sy, y_incr, x_incr)
    [sxy, sxx] = np.gradient(sx, y_incr, x_incr)
    #now get curvatures
    curv_mean_u = np.zeros((nxb,nyb))
    curv_gauss_u = np.zeros((nxb,nyb))
    for ix in range(nxb):
        for iy in range(nyb):
            #upper
            sx_c = sx[ix,iy]
            sy_c = sy[ix,iy]
            ssx = sxx[ix,iy]
            ssy = syy[ix,iy]
            ssxy = sxy[ix,iy]
            sx_v = np.array([x_incr,0.0,sx_c])
            sy_v = np.array([0.0,y_incr,sy_c])
            ssx_v = np.array([x_incr,0.0,ssx])
            ssy_v = np.array([0.0,y_incr,ssy])
            ssxy_v = np.array([0.0,y_incr,ssxy])
            E = np.dot(sx_v,sx_v)
            F = np.dot(sx_v,sy_v)
            G = np.dot(sy_v,sy_v)
            n = np.cross(sx_v,sy_v)
            n /=np.linalg.norm(n)
            L = np.dot(ssx_v,n)
            M = np.dot(ssxy_v,n)
            N = np.dot(ssy_v,n)
            #mean curvature
            J = (E*N+G*L-2.0*F*M)/(2.0*(E*G-F)**2)
            #Gaussian curvature
            K = (L*N-M**2)/(E*G-F**2)
            curv_mean_u[ix,iy] = J
            curv_gauss_u[ix,iy] = K
            # print("ix: {} iy: {} J: {} K: {}".format(ix,iy,J,K))


    return (curv_mean_u,curv_gauss_u)


if __name__ == '__main__':
    test_lipid_grid_curvature()
