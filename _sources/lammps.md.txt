# The LAMMPS stuff

The LAMMPS input script for generating the dynamical matrix at different temperatures and pressures.

## Index variables

- `NUC` number of lattice cells in each direction. Assumed to be the same for all directions. Required for monitering the lattice constant change and calculating the stress-strain relationship
- `TEMP` the temperature in Kelvins
- `ALAT0` lattice constant at zero pressure. Required to calculate strain
- `DLAT` lattice constant scaling factor depending on the hkl direction of the unit cell
- `PRESS` the pressure along the x-axis of the box

## Building the argon crystal

### The map file

[fix-phonon](https://lammps.sandia.gov/doc/fix_phonon.html) requires a **map** file in order to work. The map file is the unit cell and index for each atom. It has a specific format

   nx ny nz nk
   #comment line
   ix iy iz ik id

   20 20 20 1
   #l1 l2 l3 k atom_id
   0 0 0 0 1

- nx, ny, nz: are the number of unit cells in the x y z directions.
- nk : number of atoms in the unit cell
- ix, iy, iz, ik: indices for an atom (which unit cell it belongs to and its basis index)
- id: the id of the atom

### The primitive cell

There are many ways to generate a LAMMPS data file. We are goin to use Python and our own code.

The file was obtained by rotating the standard FCC primitive unit cell

```python

from md.utils import lmp_box
import numpy as np
lmp_box( 5.12*np.array([[0.5,0.0,0.5], [0.5, 0.5, 0], [0.0, 0.5, 0.5]]) )

```

The primitive FCC argon unit cell is then

    # fcc_prim.lmp - Fcc Ar oriented X=[110] Y=[011] Z=[101].
 
           1  atoms
           1  atom types
 
      0.00000000       3.62745779  xlo xhi
      0.00000000       3.14147060  ylo yhi
      0.00000000       2.96180688  zlo zhi
      1.81372889       1.81372889  1.04715687 xy xz yz

    Masses
 
           1   39.94800000    # Ar
            
    Atoms # atomic
 
         1    1        0.00000000       0.00000000       0.00000000

To generate a 20x20x20 super cell from this unit cell we use

```bash

python md.pre.build_atoms fcc_prim.lmp 20 20 20 ar_prim.lmp

```

## Simulation steps

1. The simulation reads the data file **argon_prim.lmp** which has the coordinates of the atoms.
2. AN LJ interatomic potential is set up using values obtained from literature
3. The atoms are relaxed at the required temperature and pressure for 200ps
4. The positions of the atoms are sampled every 10 steps (0.05ps) and the correlation calculated over 5000steps (25ps)

## The fix-phonon parameters

fix-phonon produces better results (smoother data) when run for very long times. After taking a look at the examples in the LAMMPS source files we set the following 

    fix 1 all npt temp ${TEMP} ${TEMP} 0.5 iso ${PRESS} ${PRESS} 1.0 nreset 10000 pchain 8 drag 1.0
    ...
    variable  nequil equal ${TEQ}/dt
    variable  ncorr  equal ${TCOR}/dt
    variable  nsteps equal ${nequil}+5*(10*${ncorr})
    ...
    fix 3 all phonon   10 ${ncorr} ${nequil} map.argon_prim.lmp ArPrim.${PRESS} nasr 100
    ...
    run ${nsteps}

- `TEMP` and `PRESS` are obviously the temperature and pressure for the simulation. The temperature was set at 10K, and the pressure was varied between -500bar and 500bar.

- `TEQ` and `TCOR` are the durations for the equilibration and correlation periods. They were chosen to be 500ps and 100ps respectively.

- The simulation is ran for one `TEQ` plus five (10*`TEQ`) picoseconds. The factor of 10 is because the positions are samples every 10 steps, so the sampling window is actually 1 nanosecond.

- This will give you 5 measurments of the dynamical matrix, but they are not averaged in post. I don't believe this is done in phana.
