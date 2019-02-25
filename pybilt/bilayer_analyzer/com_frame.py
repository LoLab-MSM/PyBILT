#imports
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from builtins import object
import MDAnalysis as mda
import numpy as np

import sys
from six.moves import range


class LipidCOM(object):
    """ A lipid center of mass (COM) object.

    This object stores the COM coordinates of a lipid (or other molecule or group
    of atoms) computed from both the wrapped and unwrapped atomic coordinates. This
    object also stores information about the type of lipid as well as the total mass
    of the lipid.
    """
    def __init__(self):
        """LipidCOM initialization

        Attributes:
            type (str): The lipid type (e.g. the lipid could be typed b resname).
            com (np.array): The length three vector holding the wrapped xyz coordinates.
            com_unwrap (np.array): The length three vector holding the unwrapped xyz coordinates.
            mass (float): The total mass of the atoms used to define LipidCOM.
        """
        # lipid type/resname or other name
        self.type="UNK"
        # wrapped coordinates
        self.com=np.zeros(3)
        # unwrapped coordinates
        self.com_unwrap=np.zeros(3)
        # total mass
        self.mass=1.0
        self.leaflet = "UNK"
        self.resid = 0
        return
    # The name of this function could be changed to be more desriptive, e.g.
    # extract_com_mda_residue
    def extract(self, mda_residue, unwrap=False, box=None, name_dict=None, make_whole=True):
        """ Get the center of mass coordinates from an MDAnalysis residue

        This function calls the MDAnalysis member function center_of_mass() of the residue
        to compute the center of mass of the atoms constituting the residue.

        Args:
            mda_residue (MDAnalysis.residue): An MDAnalysis residue object from which to extract
                a center of masss coordinates.
            unwrap (bool, optional): Define which com container to store coordiates in.
                False (default) - The COM coordinates are stored in the
                container designated for the wrapped coordinate representation.
                True - The COM coordinates are stored in the container designated
                for the unwrapped coordinate representation.
        """
        self.type = mda_residue.resname
        self.resid = mda_residue.resid
        self.atom_names = mda_residue.atoms.names
        atom_group = mda_residue.atoms

        if isinstance(name_dict, dict):

            names = name_dict[self.type]

            self.atom_names = names
            n_names = len(names)
            atom_group = mda.core.AtomGroup.AtomGroup([eval("mda_residue.atoms."+names[0])])
            #if (not unwrap) and (n_names > 1):
            #    mda.lib.mdamath.make_whole(mda_residue.atoms, reference_atom=atom_group)

            for i in range(1, n_names):
                atom_group+=eval("mda_residue.atoms."+names[i])

        #else:
        if (not unwrap) and make_whole:
            mda.lib.mdamath.make_whole(mda_residue.atoms)
        if unwrap:
            self.com_unwrap = atom_group.center_of_mass()

        else:

            if box is not None:
                self.com = atom_group.center_of_mass()
                self.com_unwrap = self.com[:]
            else:
                self.com = atom_group.center_of_mass()
                self.com_unwrap = self.com[:]


        self.mass = atom_group.total_mass()
        return

# a Center of Mass frame object
class COMFrame(object):
    """ A molecular dynamics style Frame object for LipidCOM objects.

    Atrributes:
        lipidcom (list of obj:LipidCOM): A list of the LipidCOM objects assigned to the COMFrame.
        box (np.array): A 3 element vector containing the (rectangular) xyz box edge lengths.
        time (float): The simulation time that this Frame represents.
        number (int): The frame number of this Frame.
        mdnumber (int): The corresponding frame number in the original MD trajectory


    """

    # does not check that nlipids is an int
    def __init__(self, mda_frame, mda_bilayer_selection, unwrap_coords,
                 name_dict=None, multi_bead=False, rewrap=False,
                 make_whole=False):
        """ Frame initialization.

        Args:
            mda_frame (MDAnalysis.Timestep): The MDAnalysis frame for the coordinates to use in computing the
                lipid centers of mass.
            mda_bilayer_selection (MDAnalysis.AtomSelection): The MDAnalsysis atom selection object for the
                containing the atoms in the bilayer.
            unwrap_coords (numpy.Array): A numpy array containing the the unwrapped coordinates of the bilayer selection
                atoms.
        """
        # list to store the nlipids LipidCOM objects
        self.lipidcom = []
        # box dimensions -- assumes the box originates a 0,0,0
        # It might be worth adding functionality to specifiy the box origin (or center)
        # This also assumes a rectangular box
        self.box = mda_frame.dimensions[0:3]
        # simulation time
        self.time = mda_frame.time
        # frame number in the MD trajectory
        self.mdnumber = mda_frame.frame
        self.number = self.mdnumber
        self._multi_bead = multi_bead
        self._name_dict = name_dict
        self._rewrapped = rewrap
        self._make_whole = make_whole
        if (name_dict is not None) and multi_bead:
            self._build_multi_bead(mda_frame, mda_bilayer_selection,
                                   unwrap_coords, name_dict)
        else:
            self._build_single_bead(mda_frame, mda_bilayer_selection,
                                    unwrap_coords, name_dict=name_dict)
        if rewrap:
            self._rewrap()

        return

    def _build_single_bead(self, mda_frame, mda_bilayer_selection,
                           unwrap_coords, name_dict=None):
        nlipids = len(mda_bilayer_selection.residues)
        # initialize all the LipidCOM objects
        for dummy_i in range(nlipids):
            self.lipidcom.append(LipidCOM())
        #atom indices in mda selection/frame
        index = mda_bilayer_selection.indices
        # loop over the residues (lipids) and get the centers of mass
        ## do the wrapped coordinates

        r=0
        for res in mda_bilayer_selection.residues:
            #print(res," ",res.center_of_mass())
            self.lipidcom[r].extract(res, name_dict=name_dict, make_whole=self._make_whole)
            #self.lipidcom[r].mass = res.total_mass()
            r+=1
        #now unwrapped coordinates


        mda_frame._pos[index] = unwrap_coords[:]
        #now we need to adjust for the center of mass motion of the membrane -- for simplicity set all frames to (0,0,0)
        # to remove center of mass motion of the membrane
        total_mass = mda_bilayer_selection.masses.sum()
        mem_com_x = (mda_frame.positions[index][:,0]*mda_bilayer_selection.masses).sum()/total_mass
        mem_com_y = (mda_frame.positions[index][:,1]*mda_bilayer_selection.masses).sum()/total_mass
        mem_com_z = (mda_frame.positions[index][:,2]*mda_bilayer_selection.masses).sum()/total_mass
        mem_com = np.array([mem_com_x, mem_com_y, mem_com_z])
        mda_frame._pos[index] -= mem_com
        self.mem_com = mem_com
        r=0
        for res in mda_bilayer_selection.residues:
            self.lipidcom[r].extract(res, unwrap=True, name_dict=name_dict, make_whole=self._make_whole)
            r+=1
        return

    def _build_multi_bead(self, mda_frame, mda_bilayer_selection, unwrap_coords, name_dict):
        nlipids = len(mda_bilayer_selection.residues)

        #atom indices in mda selection/frame
        index = mda_bilayer_selection.indices
        # loop over the residues (lipids) and get the centers of mass
        ## do the wrapped coordinates

        r=0
        for res in mda_bilayer_selection.residues:
            #print(res," ",res.center_of_mass())
            resname = res.resname
            # initialize all the LipidCOM objects
            for atom in name_dict[resname]:
                name_dict_single = {resname:[atom]}
                self.lipidcom.append(LipidCOM())
                self.lipidcom[r].extract(res, name_dict=name_dict_single, make_whole=self._make_whole)
            #self.lipidcom[r].mass = res.total_mass()
                r+=1
        #now unwrapped coordinates


        mda_frame._pos[index] = unwrap_coords[:]
        #now we need to adjust for the center of mass motion of the membrane -- for simplicity set all frames to (0,0,0)
        # to remove center of mass motion of the membrane
        total_mass = mda_bilayer_selection.masses.sum()
        mem_com_x = (mda_frame.positions[index][:,0]*mda_bilayer_selection.masses).sum()/total_mass
        mem_com_y = (mda_frame.positions[index][:,1]*mda_bilayer_selection.masses).sum()/total_mass
        mem_com_z = (mda_frame.positions[index][:,2]*mda_bilayer_selection.masses).sum()/total_mass
        mem_com = np.array([mem_com_x, mem_com_y, mem_com_z])
        mda_frame._pos[index] -= mem_com
        self.mem_com = mem_com
        r=0
        for res in mda_bilayer_selection.residues:

            #print(res," ",res.center_of_mass())
            resname = res.resname
            # initialize all the LipidCOM objects
            for atom in name_dict[resname]:
                name_dict_single = {resname:[atom]}
                self.lipidcom[r].extract(res, unwrap=True, name_dict=name_dict_single, make_whole=self._make_whole)
            #self.lipidcom[r].mass = res.total_mass()
                r+=1
        return

    def _rewrap(self):
        ncom = len(self.lipidcom)
        box_low = np.array([0.0, 0.0, 0.0])
        box_high = self.box
        # print(box_low)
        # print(box_high)
        # print(self.box)
        # print(self.mem_com)
        n_re = 0
        for r in range(ncom):
            pos_c = self.lipidcom[r].com_unwrap + self.mem_com
            pos_n = np.zeros(3)
            diff = pos_c - self.lipidcom[r].com
            dist = np.sqrt(np.dot(diff, diff))
            # is_wrap = False
            if dist > 1.0:
                n_re += 1
                for i in range(3):
                    p_i = pos_c[i]
                    b_l = box_low[i]
                    b_h = box_high[i]
                    b = self.box[i]

                    while p_i < b_l:
                        p_i += b
                        # is_wrap = True
                    while p_i > b_h:
                        p_i -= b
                        # is_wrap = True
                    pos_n[i] = p_i
                # print(r, pos_c, pos_n)
                # print('com', self.lipidcom[r].com)
                # print('com_u', self.lipidcom[r].com_unwrap)
                # print('bl bh',box_low, box_high)
                # print('mc b', self.mem_com, self.box)
                self.lipidcom[r].com = pos_n
        print(('n_re', n_re))
        return

    def __repr__(self):
        return 'COMFrame for frame %s with %s lipids' % (self.number, len(self.lipidcom))

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
        """ Returns the number of LipidCOM objects assigned to this Frame

        Returns:
            int: Number of LipidCOM objects currently assigned to this Frame
        """
        return len(self.lipidcom)

#    def COG(self,unwrapped=False):
#        cog_out = np.zeros(3)
#        for lipid in self.lipidcom:
#            if not unwrapped:
#                cog_out+=lipid.com
#            else:
#                cog_out+=lipid.com_unwrap
#        cog_out/=len(self)
#        return com_out

    def com(self, wrapped=True):
        """ Computes the center of mass (COM) for the Frame

        This member function is used to compute the overall center of mass (COM) of the
        COMFrame using the LipidCOM object coordinates and masses.

        Args:
            wrapped (bool, optional): Define which set of coordinates to use in the computation.
                True (default) - The wrapped LipidCOM coordinates are used to compute
                the COM of the frame.
                False - The unwrapped LipidCOM coordinates are used to compute
                the COM of the frame.

        Returns:
            np.array: A 3 element vector containing the xyz coordinates of the Frame's COM
        """
        com_out = np.zeros(3)
        total_mass = 0.0
        for lipid in self.lipidcom:
            if wrapped:
                com_out+=lipid.com*lipid.mass
                total_mass+=lipid.mass
            else:
                com_out+=lipid.com_unwrap*lipid.mass
                total_mass+=lipid.mass
        com_out/=total_mass
        return com_out

    def cog(self, wrapped=True):
        """ Computes the center of geometry (COG) for the Frame

        This member function is used to compute the overall center of geometry (COG) of the
        COMFrame using the LipidCOM object coordinates.

        Args:
            wrapped (bool, optional): Define which set of coordinates to use in the computation.
                True (default) - The wrapped LipidCOM coordinates are used to compute
                the COG of the frame.
                False - The unwrapped LipidCOM coordinates are used to compute
                the COG of the frame.

        Returns:
            np.array: A 3 element vector containing the xyz coordinates of the Frame's COG
        """
        com_out = np.zeros(3)
        total_mass = 0.0
        for lipid in self.lipidcom:
            if wrapped:
                com_out+=lipid.com
                total_mass+=1.0
            else:
                com_out+=lipid.com_unwrap
                total_mass+=1.0
        com_out/=total_mass
        return com_out

    def coordinates(self, wrapped=True, leaflet=None):
        coords = []
        if leaflet is None:

            for item in self.lipidcom:
                if wrapped:
                    coords.append(item.com)
                else:
                    coords.append(item.com_unwrap)
        else:
            for item in self.lipidcom:
                if item.leaflet == leaflet:
                    if wrapped:
                        coords.append(item.com)
                    else:
                        coords.append(item.com_unwrap)            
        return np.array(coords)

    def masses(self):
        output = []
        for lipid in self.lipidcom:
            output.append(lipid.mass)
        return np.array(output)

    def resids(self):
        output = []
        for lipid in self.lipidcom:
            output.append(lipid.resid)
        return np.array(output)

    def resnames(self):

        output = []
        for lipid in self.lipidcom:
            output.append(lipid.type)
        return output

    def leaflets(self):
        output = []
        for lipid in self.lipidcom:
            output.append(lipid.leaflet)
        return output

    def unique_resnames(self):

        output = []
        for lipid in self.lipidcom:
            if lipid.type not in output:
                output.append(lipid.type)
        return output

    def write_xyz(self, xyz_name, wrapped=True, name_by_leaflet=False):
        # Open up the file to write to
        xyz_out = open(xyz_name, "w")

        comment = "COMFrame "+str(self.number)+" MD Frame "+str(self.mdnumber)
        xyz_out.write(str(len(self.lipidcom)))
        xyz_out.write("\n")
        xyz_out.write(comment)
        xyz_out.write("\n")


        i=0
        for dummy_lip in self.lipidcom:
            #get the coordinates
            x = self.lipidcom[i].com[0]
            y = self.lipidcom[i].com[1]
            z = self.lipidcom[i].com_unwrap[2]
            if not wrapped:
                x = self.lipidcom[i].com_unwrap[0]
                y = self.lipidcom[i].com_unwrap[1]
                #z = self.lipidcom[i].com_unwrap[2]

            #get the lipid resname
            oname = self.lipidcom[i].type
            if name_by_leaflet:
                oname = self.lipidcom[i].leaflet
            #write to file
            line = str(oname)+" "+str(x)+" "+str(y)+" "+str(z)
            xyz_out.write(line)
            xyz_out.write("\n")
            i+=1

        xyz_out.close()
        return
