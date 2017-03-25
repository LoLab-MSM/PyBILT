#Setup using an input file

import vorbilt.bilayer_analyzer.bilayer_analyzer as ba


analyzer = ba.BilayerAnalyzer(input_file='../sample_input_script/sample_1.in')

analyzer.run_analysis()

analyzer.show_plot('msd_p')

