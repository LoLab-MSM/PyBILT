"""Defines the MDAData class that is used to store the MDAnalysis objects built from the structure and trajectory file.
"""

import MDAnalysis as mda

class MDAData:
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

    def update_trajectory(self, trajectory_file):
        self.mda_universe.load_new(trajectory_file)
        self.trajectory_file = trajectory_file
        self.mda_trajectory = self.mda_universe.trajectory
        self.nframes = len(self.mda_trajectory)
        return
