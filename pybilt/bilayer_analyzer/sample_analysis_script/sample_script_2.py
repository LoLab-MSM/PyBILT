#Same protocol as sample_script_1 without using the input setup file.

import pybilt.bilayer_analyzer.bilayer_analyzer as ba

analyzer = ba.BilayerAnalyzer(structure='../../sample_bilayer/sample_bilayer.psf',
    trajectory='../../sample_bilayer/sample_bilayer_10frames.dcd',
    selection='resname POPC DOPE TLCL')
analyzer.add_analysis('msd msd_1')
analyzer.add_analysis('msd msd_2 leaflet upper resname POPC')
analyzer.add_plot('msd msd_p msd_1 DOPE-U msd_2 POPC-U')
analyzer.run_analysis()

analyzer.show_plot('msd_p')
