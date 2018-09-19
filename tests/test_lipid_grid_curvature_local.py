from __future__ import print_function
from pybilt.bilayer_analyzer import BilayerAnalyzer
import numpy as np
from pybilt.plot_generation import plot_generation_functions as pgf
from pybilt.plot_generation.plot_generation_functions import _color_list
import scipy.interpolate
from scipy.ndimage.filters import gaussian_filter
from six.moves import range
def test_lipid_grid_curvature():
    sel_string = "resname POPC DOPE TLCL2"
    name_dict = {'DOPE':['P'],'POPC':['P'],'TLCL2':['P1','P3']}
    #analyzer = BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
    #                              trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
    #                              selection="resname POPC DOPE TLCL2")
    analyzer = BilayerAnalyzer(structure='/home/blake/Downloads/popc_membrane_nowat.psf',
                                  trajectory='/home/blake/Downloads/popc_membrane_nowat_curve_e.pdb',
                                  selection="resname POPC")
    analyzer.remove_analysis('msd_1')
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    analyzer._analysis_protocol.use_objects['lipid_grid'] = True
    # adjust the number of bins for gridding
    analyzer.rep_settings['lipid_grid']['n_xbins'] =100
    analyzer.rep_settings['lipid_grid']['n_ybins'] = 100
    analyzer.adjust_rep_setting('leaflets', 'assign_method', 'orientation')
    analyzer.adjust_rep_setting('leaflets', 'orientation_atoms', {'DOPE': ['C218','P'],
                                'POPC': ['C218', 'P'], 'TLCL2': ['CA18', 'P1']})
    #analyzer.set_frame_range(0,1,1)
    i =0
    zgrid = None
    for _frame in analyzer:
        if i == 0:
            zgrid = analyzer.reps['lipid_grid'].leaf_grid['upper'].lipid_grid_z
        else:
            zgrid += analyzer.reps['lipid_grid'].leaf_grid['upper'].lipid_grid_z
        i+=1
        #thickgrid = analyzer.reps['lipid_grid'].thickness_grid()
        #xyzc = analyzer.reps['lipid_grid'].get_xyzc(leaflet='lower', color_grid=thickgrid)['lower']

        # with sns.color_palette("PuBuGn_d"):
        #pgf.plot_lipid_grid_thickness_map_2d(xyzc[0], xyzc[1], thickgrid, save=False, show=False,
        #                                     vmin=10.0, vmax=35.0, interpolation='gaussian')
        #curvature_grids = analyzer.reps['lipid_grid'].curvature()
        #break
        #xyzc_u_mean = analyzer.reps['lipid_grid'].get_xyzc(leaflet='upper',color_grid=curvature_grids[0][0])['upper']
        # print(xyzc_u_mean)
        #pgf.plot_grid_as_scatter(xyzc_u_mean, save=False, show=False, colorbar=True)
    analyzer.reps['lipid_grid'].write_xyz()
    analyzer.reps['lipid_grid'].write_pdb('./lg_local_t_e2.pdb', use_gaussian_filter=True, filter_sigma=10.0)
    xmax = max(analyzer.reps['lipid_grid'].leaf_grid['upper'].x_edges)
    xmin = min(analyzer.reps['lipid_grid'].leaf_grid['upper'].x_edges)
    ymax = max(analyzer.reps['lipid_grid'].leaf_grid['upper'].y_edges)
    grid_x, grid_y = np.mgrid[xmin:xmax:100j, 0:ymax:100j]
    xyzc = analyzer.reps['lipid_grid'].get_xyzc(leaflet='upper', color_grid=zgrid/float(i))['upper']
    #zgrid_i = scipy.interpolate.griddata((xyzc[0], xyzc[1]), xyzc[3],(grid_x, grid_y), method='nearest')
    zgrid_i = zgrid/float(i)
    print(zgrid_i)
    zgrid_i_f = gaussian_filter(zgrid_i, 5.0, mode="nearest")
    print(zgrid_i_f)
    #grid_x_b, grid_y_b = np.mgrid[0:xmax:400j, 0:ymax:400j]
    #zgrid_i = scipy.interpolate.griddata((grid_x.flatten(), grid_y.flatten()), zgrid_i.flatten(),(grid_x_b, grid_y_b), method='linear')
    #print(zgrid_i)
    #print(zgrid_i.shape)
    #print(len(grid_x[:,0]))
    #print(grid_y[0,:])
    curvature_grids_i = curvature_b(grid_x[:,0],
                                grid_y[0,:],
                                zgrid_i_f)
    #quit()
    #print(curvature_grids[0])
    curvature_grids = curvature_b(analyzer.reps['lipid_grid'].leaf_grid['upper'].x_centers,
                                analyzer.reps['lipid_grid'].leaf_grid['upper'].y_centers,
                                zgrid/float(i))
    xyzc_u_mean = analyzer.reps['lipid_grid'].get_xyzc(leaflet='upper',color_grid=curvature_grids[0])['upper']
    print((len(curvature_grids_i[0].flatten())))
    print((len(grid_x.flatten())))
    print((max(grid_x.flatten())))
    print((grid_x.flatten()[-1]))
    xyzc_u_mean_i = (grid_x.flatten(), grid_y.flatten(), grid_x.flatten(), curvature_grids_i[0].flatten())
    pgf.plot_grid_as_scatter(xyzc_u_mean_i, save=False, show=False, colorbar=True)
    pgf.plot_lipid_grid_thickness_map_2d(xyzc_u_mean_i[0], xyzc_u_mean_i[1], curvature_grids_i[0], save=False, show=False, vmin=-0.08, vmax=0.08)
    pgf.plot_lipid_grid_thickness_map_2d(xyzc_u_mean[0], xyzc_u_mean[1], curvature_grids[0], save=False, show=False, vmin=-0.08, vmax=0.08)
    xyzc = (grid_x.flatten(), grid_y.flatten(), grid_x.flatten(), zgrid_i.flatten())
    pgf.plot_grid_as_scatter(xyzc, save=False, show=False, colorbar=True)
    xyzc = (grid_x.flatten(), grid_y.flatten(), grid_x.flatten(), zgrid_i_f.flatten())
    pgf.plot_grid_as_scatter(xyzc, save=False, show=False, colorbar=True)
    return


def curvature_b(x_vals, y_vals, zgrid):
    nxb = len(x_vals)
    nyb = len(y_vals)
    x_incr = x_vals[1]-x_vals[0]
    y_incr = y_vals[1]-y_vals[0]
    print(("x_incr {} y_incr {}".format(x_incr, y_incr)))
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

def curvature(x_vals, y_vals, zgrid):
    nxb = len(x_vals)
    nyb = len(y_vals)
    x_incr = x_vals[1]-x_vals[0]
    y_incr = y_vals[1]-y_vals[0]
    #first order derivtives
    sx_u = np.zeros((nxb,nyb))
    sy_u = np.zeros((nxb,nyb))
    for ix in range(nxb):
        for iy in range(nyb):
            ixp = ix-1
            if ixp < 0:
                ixp+=nxb
            ixn = ix+1
            if ixn >= nxb:
                ixn-=nxb
            iyp = ix-1
            if iyp < 0:
                iyp+=nyb
            iyn = iy+1
            if iyn >= nyb:
                iyn-=nyb
            #upper
                ## using central difference for numerical first derivative
            sx = zgrid[ixn,iy]-zgrid[ixp,iy]
            sx/= (x_incr)**2
            sy = zgrid[ix,iyn]-zgrid[ix,iyp]
            sy/= (y_incr)**2
            sx_u[ix,iy]=sx
            sy_u[ix,iy]=sy
    #now do second order derivatives - central difference numerical derivative of the first derivative
    ssx_u = np.zeros((nxb,nyb))
    ssy_u = np.zeros((nxb,nyb))
    ssxy_u = np.zeros((nxb,nyb))
    for ix in range(nxb):
        for iy in range(nyb):
            ixp = ix-1
            if ixp < 0:
                ixp+=nxb
            ixn = ix+1
            if ixn >= nxb:
                ixn-=nxb
            iyp = ix-1
            if iyp < 0:
                iyp+=nyb
            iyn = iy+1
            if iyn >= nyb:
                iyn-=nyb
            #upper
                ## using central difference for numerical first derivative
            ssx = sx_u[ixn,iy]-sx_u[ixp,iy]
            ssx/= (x_incr)**2
            ssy = sy_u[ix,iyn]-sy_u[ix,iyp]
            ssy/= (y_incr)**2
            ssxy = sx_u[ix,iyn]-sx_u[ix,iyp]
            ssxy/=(y_incr)**2
            ssx_u[ix,iy]=ssx
            ssy_u[ix,iy]=ssy
            ssxy_u[ix,iy]=ssxy

    #now get curvatures
    curv_mean_u = np.zeros((nxb,nyb))
    curv_gauss_u = np.zeros((nxb,nyb))
    for ix in range(nxb):
        for iy in range(nyb):
            #upper
            sx = sx_u[ix,iy]
            sy = sy_u[ix,iy]
            ssx = ssx_u[ix,iy]
            ssy = ssy_u[ix,iy]
            ssxy = ssxy_u[ix,iy]
            sx_v = np.array([x_vals[ix]+x_incr,0.0,sx])
            sy_v = np.array([0.0,y_vals[iy]+y_incr,sy])
            ssx_v = np.array([x_vals[ix]+x_incr,0.0,ssx])
            ssy_v = np.array([0.0,y_vals[iy]+y_incr,ssy])
            ssxy_v = np.array([0.0,y_vals[iy]+y_incr,ssxy])
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
            print(("ix: {} iy: {} J: {} K: {}".format(ix,iy,J,K)))


    return (curv_mean_u,curv_gauss_u)

if __name__ == '__main__':
    test_lipid_grid_curvature()
