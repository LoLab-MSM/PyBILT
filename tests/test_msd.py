import vorbilt.bilayer_analyzer.bilayer_analyzer as ba
import vorbilt.com_trajectory.COMTraj as comtraj
import MDAnalysis as mda
import numpy as np
def test_msd():
    analyzer = ba.BilayerAnalyzer(psf_file='../vorbilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../vorbilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="not resname CLA and not resname TIP3 and not resname POT")


    analyzer.print_analysis_protocol()

    analyzer.run_analysis()
    msd_dat_a = analyzer.get_analysis_data('msd_1')
    print("MSD from bilayer analyzer:")
    print(msd_dat_a)

    u = mda.Universe('../vorbilt/sample_bilayer/sample_bilayer.psf', '../vorbilt/sample_bilayer/sample_bilayer_10frames.dcd')
    bilayer_sel = u.select_atoms("not resname CLA and not resname TIP3 and not resname POT")
    ct = comtraj.COMTraj(u.trajectory,bilayer_sel)
    msd_dat_b = ct.calc_msd()
    print("MSD from COMTraj:")
    print(msd_dat_b)

    print("matching: {}".format(np.isclose(msd_dat_a, msd_dat_b, rtol=0.001)))
if __name__ == '__main__':
    test_msd()



