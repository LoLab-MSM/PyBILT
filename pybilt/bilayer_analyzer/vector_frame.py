#imports
import MDAnalysis as mda
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import pybilt.plot_generation.plot_generation_functions as pgf
import matplotlib.pyplot as plt
import sys

#range/xrange fix
if sys.version_info < (3,0):
    def range(*args, **kwargs):
        return xrange(*args, **kwargs)

class LipidVec(object):
    """ A lipid vector object.

    This object stores the coordinates corresponding to a vector representation
    of a lipid (e.g. the PN vector); the object uses the wrapped coordinates.
    This object also stores information about the type of lipid as well as the
    total mass of the lipid.
    """
    def __init__(self):
        """LipidVec initialization

        Attributes:
            resname (str): The lipid resname.
            start (np.array): The length three vector holding the wrapped xyz
                coordinates of the vector starting point.
            end (np.array): The length three vector holding the wrapped xyz
                coordinates of the vector end point.
            vector (np.array): The length three vector (i.e end-start).
            mass (float): The total mass of the atoms used to define LipidVec.
            leaflet (str): The leaflet that the lipid is in can be stored here.
            resid (int): The integer resid of the lipid.
        """
        # lipid type/resname or other name
        self.resname="UNK"
        # wrapped coordinates
        self.start=np.zeros(3)
        self.end=np.zeros(3)
        self.vector = np.zeros(3)
        self.start_mass=1.0
        self.end_mass = 1.0
        self.mass = 1.0
        self.leaflet = "UNK"
        self.resid = 0
        return
    # The name of this function could be changed to be more desriptive, e.g.
    # extract_com_mda_residue
    def extract(self, mda_residue, ref_atoms):
        """ Get the center of mass coordinates from an MDAnalysis residue

        This function calls the MDAnalysis member function center_of_mass() of the residue
        to compute the center of mass of the atoms constituting the residue.

        Args:
            mda_residue (MDAnalysis.residue): An MDAnalysis residue object from
                which to extract (build) a vector representation.
            ref_atoms (dict): A dictionary storing the references atoms to be
                used in building the vector representation. The outer keys of
                ref_atoms should be the resnames of the lipids.
        """
        self.resname = mda_residue.resname
        self.resid = mda_residue.resid
        if ref_atoms is not None:
            self.ref_atoms = ref_atoms[self.resname]
            #start
            start_atoms = self.ref_atoms['start']
            if isinstance(start_atoms, (list, tuple)):

                atom_group = mda.core.AtomGroup.AtomGroup([eval("mda_residue.atoms."+start_atoms[0])])
                n_names = len(start_atoms)
                for i in range(1, n_names):
                    atom_group+=eval("mda_residue.atoms."+start_atoms[i])
                self.start = atom_group.center_of_mass()
                self.start_mass = atom_group.total_mass()
            else:
                atom_group = mda.core.AtomGroup.AtomGroup([eval("mda_residue.atoms."+start_atoms)])
                self.start = atom_group.center_of_mass()
                self.start_mass = atom_group.total_mass()

            end_atoms = self.ref_atoms['end']
            if isinstance(end_atoms, (list, tuple)):

                atom_group = mda.core.AtomGroup.AtomGroup([eval("mda_residue.atoms."+end_atoms[0])])
                n_names = len(end_atoms)
                for i in range(1, n_names):
                    atom_group+=eval("mda_residue.atoms."+end_atoms[i])
                self.end = atom_group.center_of_mass()
                self.end_mass = atom_group.total_mass()
            else:
                atom_group = mda.core.AtomGroup.AtomGroup([eval("mda_residue.atoms."+end_atoms)])
                self.end = atom_group.center_of_mass()
                self.end_mass = atom_group.total_mass()
        else:
            atom_group = mda_residue.atoms
            min_atom_index = np.argmin(atom_group.positions[:, 2])
            max_atom_index = np.argmax(atom_group.positions[:, 2])
            self.ref_atoms = {'start':atom_group[min_atom_index].name,
                              'end': atom_group[max_atom_index].name}
            self.end = atom_group[max_atom_index].position
            self.end_mass = atom_group[max_atom_index].mass
            self.start = atom_group[min_atom_index].position
            self.start_mass = atom_group[max_atom_index].mass
            # print(self.ref_atoms)
        self.vector = self.end - self.start
        self.mass = mda_residue.atoms.total_mass()
        return

# a Vector frame object
class VectorFrame(object):
    """ A molecular dynamics style Frame object for LipidVec objects.

    Atrributes:
        lipidvec (list of obj:LipidVec): A list of the LipidVec objects assigned to the COMFrame.
        box (np.array): A 3 element vector containing the (rectangular) xyz box edge lengths.
        time (float): The simulation time that this Frame represents.
        number (int): The frame number of this Frame.
        mdnumber (int): The corresponding frame number in the original MD trajectory


    """

    # does not check that nlipids is an int
    def __init__(self, mda_frame, mda_bilayer_selection, ref_atoms):
        """ Frame initialization.

        Args:
            mda_frame (MDAnalysis.Timestep): The MDAnalysis frame for the coordinates to use in computing the
                lipid centers of mass.
            mda_bilayer_selection (MDAnalysis.AtomSelection): The MDAnalsysis atom selection object for the
                containing the atoms in the bilayer.
            unwrap_coords (numpy.Array): A numpy array containing the the unwrapped coordinates of the bilayer selection
                atoms.
        """
        # list to store the nlipids LipidVec objects
        self.lipidvec = []
        self.ref_atoms = ref_atoms
        # box dimensions -- assumes the box originates a 0,0,0
        # It might be worth adding functionality to specifiy the box origin (or center)
        # This also assumes a rectangular box
        self.box = mda_frame.dimensions[0:3]
        # simulation time
        self.time = mda_frame.time
        # frame number in the MD trajectory
        self.mdnumber = mda_frame.frame
        self.number = self.mdnumber
        nlipids = len(mda_bilayer_selection.residues)
        # initialize all the LipidVec objects
        for dummy_i in range(nlipids):
            self.lipidvec.append(LipidVec())
        #atom indices in mda selection/frame
        index = mda_bilayer_selection.indices
        # loop over the residues (lipids) and get the centers of mass
        ## do the wrapped coordinates

        r=0
        for res in mda_bilayer_selection.residues:
            #print(res," ",res.center_of_mass())
            self.lipidvec[r].extract(res, ref_atoms)
            #self.lipidvec[r].mass = res.total_mass()
            r+=1

        return

    def __repr__(self):
        return 'VecFrame for frame %s with %s lipids' % (self.number, len(self.lipidvec))

    def set_box(self, box_lengths):
        """ Set the rectangular xyz box edge lengths.

        Args:
            box_lengths (numpy.array): A 1d, 3 element numpy.array containing the x,y,z box sizes (or edge lengths)

        """
        self.box = box_lengths
        return

    def set_time(self, time):
        """ Set the simulation time.

        Args:
            time (float): The simulation time to assign to this Frame.

        """
        self.time = time
        return

    def __len__(self):
        """ Returns the number of LipidVec objects assigned to this Frame

        Returns:
            int: Number of LipidVec objects currently assigned to this Frame
        """
        return len(self.lipidvec)

    def vectors(self):
        coords = []
        for item in self.lipidvec:
            coords.append(item.vector)
        return np.array(coords)

    def masses(self):
        output = []
        for lipid in self.lipidvec:
            output.append(lipid.mass)
        return np.array(output)

    def resids(self):
        output = []
        for lipid in self.lipidvec:
            output.append(lipid.resid)
        return np.array(output)

    def resnames(self):

        output = []
        for lipid in self.lipidvec:
            output.append(lipid.resname)
        return output

    def leaflets(self):
        output = []
        for lipid in self.lipidvec:
            output.append(lipid.leaflet)
        return output

    def unique_resnames(self):

        output = []
        for lipid in self.lipidvec:
            if lipid.resname not in output:
                output.append(lipid.resname)
        return output

    # def write_xyz(self, xyz_name, wrapped=True):
    #     # Open up the file to write to
    #     xyz_out = open(xyz_name, "w")
    #
    #     comment = "COMFrame "+str(self.number)+" MD Frame "+str(self.mdnumber)
    #     xyz_out.write(str(len(self.lipidvec)))
    #     xyz_out.write("\n")
    #     xyz_out.write(comment)
    #     xyz_out.write("\n")
    #
    #
    #     i=0
    #     for dummy_lip in self.lipidvec:
    #         #get the coordinates
    #         x = self.lipidvec[i].com[0]
    #         y = self.lipidvec[i].com[1]
    #         z = self.lipidvec[i].com_unwrap[2]
    #         if not wrapped:
    #             x = self.lipidvec[i].com_unwrap[0]
    #             y = self.lipidvec[i].com_unwrap[1]
    #             #z = self.lipidvec[i].com_unwrap[2]
    #
    #         #get the lipid resname
    #         oname = self.lipidvec[i].type
    #
    #         #write to file
    #         line = str(oname)+" "+str(x)+" "+str(y)+" "+str(z)
    #         xyz_out.write(line)
    #         xyz_out.write("\n")
    #         i+=1
    #
    #     xyz_out.close()
    #     return

    def to_quiver(self):
        pass
