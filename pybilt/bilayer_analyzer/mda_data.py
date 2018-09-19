"""Defines the MDAData class that is used to store the MDAnalysis objects built from the structure and trajectory file.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from builtins import object
import MDAnalysis as mda

class MDAData(object):
    """Store the MDAnalysis objects used in the analyses.
     Attributes:
         mda_universe (
    """
    def __init__(self,psf, trajectory, bilayer_sel_string):

        self.mda_universe = mda.Universe(psf, trajectory)
        self.mda_trajectory = self.mda_universe.trajectory
        self.psf_file = psf
        self.trajectory_file = trajectory
        self.bilayer_sel_string = bilayer_sel_string
        self.bilayer_sel = self.mda_universe.select_atoms(bilayer_sel_string)
        self.nframes = len(self.mda_trajectory)
        self.natoms = len(self.bilayer_sel)
        self.indices = self.bilayer_sel.indices
        self.residues = self.bilayer_sel.residues
        self.n_residues = len(self.residues)
        print((list(self.__dict__.keys())))

    def update_trajectory(self, trajectory_file):
        self.mda_universe.load_new(trajectory_file)
        self.trajectory_file = trajectory_file
        self.mda_trajectory = self.mda_universe.trajectory
        self.nframes = len(self.mda_trajectory)
        return

    def __getstate__(self):
        odict = dict()
        for key in self.__dict__.keys():
            odict[key] = self.__dict__[key]
        odict['mda_universe'] = None
        odict['mda_trajectory'] = None
        odict['bilayer_sel'] = None
        odict['residues'] = None
        #print('pickling....')
        return odict

    def __setstate__(self, in_dict):
        #print("unpickling...")
        #mda_universe = mda.Universe(dict['psf_file'], dict['trajectory_file'])
        #mda_universe = None
        #print(mda_universe)
        #mda_trajectory = mda_universe.trajectory
        # print(mda_trajectory)
        #bilayer_sel = mda_universe.select_atoms(dict['bilayer_sel_string'])
        #print(bilayer_sel)
        #residues = bilayer_sel.residues

        #dict['mda_universe'] = mda_universe
        #dict['mda_trajectory'] = mda_trajectory
        #dict['bilayer_sel'] = bilayer_sel
        #dict['residues'] = residues
        self.__dict__ = in_dict

    def reduced(self):
        return MDADataReduced(self)

    #def __getinitargs__(self):
    #    return (self.psf_file, self.trajectory_file, self.bilayer_sel_string)

class MDADataReduced(object):
    def __init__(self, mda_data):
        self.psf_file = mda_data.psf_file
        self.trajectory_file = mda_data.trajectory_file
        self.bilayer_sel_string = mda_data.bilayer_sel_string
        self.nframes = mda_data.nframes
        self.natoms = mda_data.natoms
        self.indices = mda_data.indices
        self.n_residues = mda_data.n_residues
        return
