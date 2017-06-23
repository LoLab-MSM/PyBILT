import pybilt.bilayer_analyzer.bilayer_analyzer as ba
import pybilt.com_trajectory.COMTraj as comtraj
import MDAnalysis as mda
import numpy as np
def test_msd():
    analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="not resname CLA and not resname TIP3 and not resname POT")


    analyzer.print_analysis_protocol()

    analyzer.run_analysis()
    msd_dat_a = analyzer.get_analysis_data('msd_1')
    print("MSD from BilayerAnalyzer:")
    print(msd_dat_a)

    u = mda.Universe('../pybilt/sample_bilayer/sample_bilayer.psf', '../pybilt/sample_bilayer/sample_bilayer_10frames.dcd')
    bilayer_sel = u.select_atoms("not resname CLA and not resname TIP3 and not resname POT")
    ct = comtraj.COMTraj(u.trajectory,bilayer_sel)
    msd_dat_b = ct.calc_msd()
    print("MSD from COMTraj:")
    print(msd_dat_b)

    print("BilayerAnalyzer and COMTraj match: {}".format(np.isclose(msd_dat_a, msd_dat_b, rtol=0.001)))

    #redo bilayer_anlayzer calc, but use the iterator looper
    analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="not resname CLA and not resname TIP3 and not resname POT")


    analyzer.print_analysis_protocol()
    for _frame in analyzer:
        pass
    msd_dat_c = analyzer.get_analysis_data('msd_1')
    print("MSD from BilayerAnalyzer obtained using iterator:")
    print(msd_dat_c)
    print("matches first BilayerAnalyzer comp: {}".format(np.isclose(msd_dat_a, msd_dat_b, rtol=0.001)))
    #print(analyzer.mda_data.bilayer_sel.residues[0].center_of_mass())


if __name__ == '__main__':
    test_msd()
