"""Create a LAMMPS input data file from another

This script return a LAMMPS input data file for running crystalline simulations
Allows the user to replicate a unit cell similar to the replicate  command in LAMMPS. 
The script also generates the map file required by fix-phonon
The user specifies an LAMMPS data file which serves as the unit cell.
The unit cell is repeated in the X,Y,Z directions.

.. todo:: add option for translating (centering), changing to another basis
.. todo:: add option to generate other types of data files (charge, dihedrals, ...)
.. todo:: print the mass section to the output file
"""

import numpy as np
import re
from ..Csystem import Csystem

def read_lmp(fin:str, verbose=False) -> Csystem : 
    """ Reads a LAMMPS data file and returns a csystem

    .. todo:: Fix atoms section to auto detect charge information

    :param fin: filename
    :type fin: str
    :param verbose: Print some info about what's going on
    :type verbose: bool
    :return: A :py:class:Csystem
    """

    def find_line_number(lst:list,key:str) -> int : 
        """ Find line number starting with a certain substring
        
        Scan a list of strings for a matching beginning.
        
        :param lst: List containing lines (from readlines())
        :type lst: list
        :type key: str
        :param key: Substring to look for
        :return: the index of the first line that has the key
        :rtype: int
        """
        for i in range(len(lst)):
            elem_str= lst[i]
            if elem_str.startswith(key): return i
        return -1

    ucelltext = open(fin).readlines()

    for i in range(len(ucelltext)): 
        ucelltext[i] = ucelltext[i].strip() # remove whitespaces
        ucelltext[i] = re.sub(' +', ' ', ucelltext[i]) # remove consecutive whitespaces

    # General section
    # Line numbers are fixed
    natoms = int(ucelltext[2].split()[0])
    ntypes = int(ucelltext[3].split()[0])
    cell = np.zeros((3,3))
    origin = np.zeros((3,1))
    cell[0,0] = float(ucelltext[5].split()[1])
    cell[1,1] = float(ucelltext[6].split()[1])
    cell[2,2] = float(ucelltext[7].split()[1])
    origin[0,0] = float(ucelltext[5].split()[0])
    origin[1,0] = float(ucelltext[6].split()[0])
    origin[2,0] = float(ucelltext[7].split()[0])

    # Check if triclinic information exists
    if ucelltext[8] != "":
        cell[0,1] = float(ucelltext[8].split()[0]) #xy
        cell[0,2] = float(ucelltext[8].split()[1]) #xz
        cell[1,2] = float(ucelltext[8].split()[2]) #yz
    # Print some useful information on the box
    if verbose:
        print("... {:d} basis atoms in unit cell".format(natoms)) 
        print("... {:d} atoms in unit cell".format(ntypes))
        print("... unit cell = \n"+
                "     ┌        ┐\n" +
                "     |  |  |  | = ({0[0]:f} {0[1]:f} {0[2]:f})\n".format(cell[0,:])+
                "     |a1|a2|a3| = ({0[0]:f} {0[1]:f} {0[2]:f})\n".format(cell[1,:])+
                "     |  |  |  | = ({0[0]:f} {0[1]:f} {0[2]:f})\n".format(cell[2,:])+
                "     └        ┘ "
        )
        print("... origin = "+ "{0:f}, {1:f}, {2:f}".format(origin[0,0],origin[1,0],origin[2,0]))

    
    # Masses esction
    masses = {}
    ui = find_line_number(ucelltext,"Masses")
    if ui >0: # No information on masses in file
        for ln in ucelltext[ui+2:ui+ntypes+2]:
            ln = ln.split()
            masses[int(ln[0])] = float(ln[1])

    
    # Atoms section
    ui = find_line_number(ucelltext,"Atoms")
    atoms = np.zeros((natoms,5))
    if ui <0:
        print("Error: no atom section -- stopping")
        quit()

    iatom = 0
    for ln in ucelltext[ui+2:ui+natoms+2]:
            ln = ln.split()
            atoms[iatom,0] = int(ln[0]) #id
            atoms[iatom,1] = int(ln[1]) #type
            atoms[iatom,2] = float(ln[2]) #x
            atoms[iatom,3] = float(ln[3]) #y
            atoms[iatom,4] = float(ln[4]) #z
            iatom += 1

    return Csystem(cell,origin,atoms,masses)


def lmpbox(oldbox:np.array ):
    """ Makes a non-LAMMPS box into LAMMPS compatible

    Boxes or unit vectors in LAMMPS have to obey certain rules. See section on `triclinic boxes <https://lammps.sandia.gov/doc/Howto_triclinic.html>`_
    The function rotates a 3x3 numpy array so that it is LAMMPS compliant following the method given on the website.

    :type oldbox : numpy.array
    :param oldbox: 3x3 array with basis vectors as the columns 
    :return: the new box rotated to comply with lammps
    :rtype: numpy.array
    """

    newbox = np.zeros((3,3))
    Ahat = oldbox[:,0]/np.linalg.norm(oldbox[:,0])
    AxB = np.cross(oldbox[:,0],oldbox[:,1])
    AxBhat = AxB/np.linalg.norm(AxB)

    newbox[0,0] = np.linalg.norm(oldbox[:,0])
    
    newbox[0,1] = np.dot(oldbox[:,1],Ahat)
    newbox[1,1] = np.linalg.norm(np.cross(Ahat,oldbox[:,1]))

    newbox[0,2] = np.dot(oldbox[:,2],Ahat)
    newbox[1,2] = np.dot(oldbox[:,2], np.cross(AxBhat,Ahat))
    newbox[2,2] = np.dot(oldbox[:,2],AxBhat)

    return newbox

def replicate(xsys: Csystem, nrep: list) -> tuple[np.array, np.array, np.array, np.array] :
    """ Replicates a crystalline system to create an array of atoms and a box
    
    :param xsys: crystalline system containing basis vectors and atoms
    :type xsys: :py:class:Csystem
    :param nrep: 3 integers representing repetition in the x,y,z directions
    :type nrep: list
    :return: the atoms, the map information, the cell or box, and the origin
    :rtype: tuple
    """
    nx = nrep[0]
    ny = nrep[1]
    nz = nrep[2]
    # Replicate
    natoms = nrep[0]*nrep[1]*nrep[2]*xsys.natoms
    atoms = np.zeros((natoms,5))
    amap = np.zeros((natoms,5),dtype=int)

    ia = 0
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for ik in range(xsys.natoms):
                    shft = np.array([[float(ix)],[float(iy)],[float(iz)]])
                    amap[ia,:] = np.array([ix,iy,iz,ik,ia+1])
                    atoms[ia,0:2] = np.array([ia+1, xsys.atoms[ik,1]])
                    atoms[ia,2:5] = (xsys.atoms[ik,2:5] +
                    np.transpose(np.dot(xsys.cell , shft)))
                    ia +=1

    cell = np.zeros((3,3))
    cell[:,0] = float(nx)*xsys.cell[:,0]
    cell[:,1] = float(ny)*xsys.cell[:,1]
    cell[:,2] = float(nz)*xsys.cell[:,2]
    origin = xsys.origin

    return atoms,amap,cell,origin

def write_from_file(fin, nx, ny, nz, foutn) -> None:
    """ Builds a crystal from input, replicating as necessary, and outputs the result to a file
    :param fin: input file name
    :type fin: str
    :param nx: number of unit cells in x-direction
    :type nx: int
    :param ny: number of unit cells in y-direction
    :type ny: int
    :param nz: number of unit cells in z-direction
    :type nz: int
    :param fout: LAMMPS output file name
    :type fout: str
    
    .. note:: function always outputs a triclinic lammps box
    .. note:: function does not check if the input lammps file has a valied box
    """
    xsys = read_lmp(fin)
    #print(xsys)
    natoms = nx*ny*nz*xsys.natoms

    atoms,amap,cell,origin = replicate(xsys, [nx,ny,nz])

    fout = open(foutn,'w')
    relpath = foutn.rfind('/')
    if relpath==-1:
        fmap = open("map."+foutn, 'w')
    else:
        fmap = open(foutn[:relpath+1]+"map."+foutn[relpath+1:], 'w')

    fmap.write("{:d} {:d} {:d} {:d}".format(nx, ny, nz, xsys.natoms))
    fmap.write("\n#l1 l2 l3 k atom_id")

    fout.write("# FCC primitive Ar")
    fout.write("\n")
    fout.write("\n\t\t{:d} atoms\n\t\t{:d} atom types".format(natoms, xsys.ntypes))
    fout.write("\n")
    fout.write("\n\t{0[0]:f} {1[0]:f} xlo xhi".format(origin[0],cell[0,:]))
    fout.write("\n\t{0[0]:f} {1[1]:f} ylo yhi".format(origin[1],cell[1,:]))
    fout.write("\n\t{0[0]:f} {1[2]:f} zlo zhi".format(origin[2],cell[2,:]))
    fout.write("\n\t{0[1]:f} {0[2]:f} {1[2]:f} xy xz yz".format(cell[0],cell[1]))
    fout.write("\n")
    fout.write("\nAtoms # atomic")
    fout.write("\n")
    for i in range(natoms):
        fout.write("\n\t\t{0:d}  {1:d}       {2[2]:f} {2[3]:f} {2[4]:f}".format(
            int(atoms[i,0]),int(atoms[i,1]),atoms[i]))
        fmap.write("\n{0[0]:d} {0[1]:d} {0[2]:d} {0[3]:d} {0[4]:d}".format(amap[i]))

    fmap.close()
    fout.close()
