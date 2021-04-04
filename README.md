# Introduction
We are going to use Python and LAMMPS to look at how the normal modes of vibration change as a function of the structural properties of Argon.

# LAMMPS
## Setup
Make sure that you install the USER-Phonon package needed to calculate the phonon properties

## Interatomic potential
We are going to use the LJ potential for Argon. The parameters are from @whiteLennardJonesModelArgon1999. The values of $\sigma$ and $\epsilon$ are

| Parameter | Value       |
| --------- | ----------- |
| $\sigma$   | 3.419 Å |
| $\epsilon$ | 0.01015122 eV |

## Chemocal properties of argon
- Melting point 
- Unit cell [WebElements](https://www.webelements.com/argon/crystal_structure.html)
  - Space group: Fm3-m
  - Space group number: 225
  - Structure: ccp (cubic close-packed) (FCC)
  - Cell parameters:
    - a: 525.6 pm
    - b: 525.6 pm
    - c: 525.6 pm
    - α: 90.000°
    - β: 90.000°
    - γ: 90.000°

We are going to use the conventional cell for FCC unit cell with a supercell of 10x10x10

## Molecular Dynamics Simulation
The equilibrium unit cell for Argon at