# -*- coding: UTF-8 -*-.
import re, numpy as np
class csystem:
    def __init__(self,acell, aorigin, abasis, amasses=None):
        self.cell = acell
        self.atoms = abasis
        self.origin = aorigin
        self.natoms = abasis.shape[0]
        self.ntypes = int(max(abasis[:,1]))
        self.masses = dict(enumerate(range(self.ntypes)))
        if amasses is not None:
            for i in amasses.keys():
                self.masses[i] = amasses[i]

    def __str__(self):
        return ("{:d} atoms, {:d} atom types\n".format(self.natoms, self.ntypes)+
        "\n cell = \n" + str(self.cell) + 
        "\n origin = \n" + str(self.origin) +
        "\n atoms = \n" +str(self.atoms))


    
"""
Reads a LAMMPS data file
"""
def read_lmp(fin, verbose=False):
    def find_line_number(lst,key):
        for i in range(len(lst)):
            elem_str= lst[i]
            if elem_str.startswith(key): return i
        return -1

    ucelltext = open(fin).readlines()

    for i in range(len(ucelltext)): 
        ucelltext[i] = ucelltext[i].strip()
        ucelltext[i] = re.sub(' +', ' ', ucelltext[i])

    """
    Read information on the unit cell
    """
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

    if ucelltext[8] != "":
        cell[0,1] = float(ucelltext[8].split()[0]) #xy
        cell[0,2] = float(ucelltext[8].split()[1]) #xz
        cell[1,2] = float(ucelltext[8].split()[2]) #yz

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

        """
        Read information on the basis atoms (masses)
        """
    masses = {}
    ui = find_line_number(ucelltext,"Masses")
    if ui >0:
        for ln in ucelltext[ui+2:ui+ntypes+2]:
            ln = ln.split()
            masses[int(ln[0])] = float(ln[1])

    """
    Read information on the basis atoms (coordinates)
    """
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

    return csystem(cell,origin,atoms,masses)


"""
Change a box to LAMMPS specifications
"""
def lmpbox(oldbox):

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
