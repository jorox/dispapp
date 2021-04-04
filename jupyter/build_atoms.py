#!/usr/bin/python

import argparse
import numpy as np
from utils import read_lmp,csystem


parser = argparse.ArgumentParser()
parser.add_argument('fin', metavar='ucell.lmp', type=str, help="Unit cell in LAMMPS format")
parser.add_argument('--options', metavar="OPT", type=str, nargs="+")
parser.add_argument('nx', metavar="N", type=int, help="Unit cells in x")
parser.add_argument('ny', metavar="N", type=int, help="Unit cells in y")
parser.add_argument('nz', metavar="N", type=int, help="Unit cells in z")
parser.add_argument('--fout', metavar='fname', 
type=str, help="Output file name", default="sucell.lmp")

# Read in the data file
args = parser.parse_args()

# Get the system
xsys = read_lmp(args.fin)
#print(xsys)

# Replicate
natoms = args.nx*args.nx*args.nx*xsys.natoms
atoms = np.zeros((natoms,5))
amap = np.zeros((natoms,5),dtype=int)

ia = 0
for ix in range(args.nx):
    for iy in range(args.ny):
        for iz in range(args.nz):
            for ik in range(xsys.natoms):
                shft = np.array([[float(ix)],[float(iy)],[float(iz)]])
                amap[ia,:] = np.array([ix,iy,iz,ik,ia+1])
                atoms[ia,0:2] = np.array([ia+1, xsys.atoms[ik,1]])
                atoms[ia,2:5] = (xsys.atoms[ik,2:5] +
                np.transpose(np.dot(xsys.cell , shft)))
                ia +=1

cell = np.zeros((3,3))
cell[:,0] = float(args.nx)*xsys.cell[:,0]
cell[:,1] = float(args.ny)*xsys.cell[:,1]
cell[:,2] = float(args.nz)*xsys.cell[:,2]
origin = xsys.origin
fout = open(args.fout,'w')
fmap = open("map."+args.fout, 'w')

fmap.write("{:d} {:d} {:d} {:d}".format(args.nx, args.ny, args.nz, xsys.natoms))
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