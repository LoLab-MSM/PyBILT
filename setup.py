from distutils.core import setup

setup(name='pybilt',
      version='0.2.0',
      description='Lipid bilayer analysis toolkit.',
      author='Blake A. Wilson',
      author_email='blake.a.wilson@vanderbilt.edu',
      url='http://pybilt.readthedocs.io/en/latest/index.html',
      packages=['pybilt', 'pybilt.bilayer_analyzer', 'pybilt.common',
                'pybilt.com_trajectory', 'pybilt.diffusion',
                'pybilt.lipid_grid', 'pybilt.mda_tools',
                'pybilt.plot_generation'],
      license='MIT',
      keywords=['lipid bilayer', 'molecular dynamics', 'analysis']            
     )
