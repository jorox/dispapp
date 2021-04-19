# Introduction

The following project is going to look at the effect of stress on the normal modes of argon. We are using LAMMPS to get the data and Python and the Bokeh library for pre and post processing.

## Setup

The  [fix-phonon](https://lammps.sandia.gov/doc/fix_phonon.html) is used to get the dynamical matrix at a set temperature and pressure. The eigenvalues of that matrix represent the frequencies of the normal modes. There should be 3 acoustic and 3N-3 optiocal modes where N is the number of atoms in our unit call. fix-phonon allows us to evaluate the normal modes at different temperatures and pressures.

The following libraries are needed:

- fftw
- open MPI

To compile LAMMPS with fix-phonon run the following under **lammps/src**

    cp MAKE/OPTIONS/lmp_g++_openmpi
    make yes-kspace
    make yes-user-phonon
    make -j4 g++_openmpi

Once you have the LAMMPS executable we are ready to write the LAMMPS script.

## Post-processing

Once the data is available we use our own python module to post-processing. This invloves 

1. Parsing the log files of fix-phonon

2. Getting the eigenvalues ei for the qi in the log file

3. Building a dispersion curve along the symmetru points for FCC

4. Interpolating the data using splines

5. Creating the interactive plot

## What about phana ?

Phana is a program that reads the binary file produced by fix-phonon. It produces dispersion/dos data and interpolates them. It is written in C++ by the author of fix-phonon.

Compiling phana was a bit of a problem on a standard linux machine. It required several libraries, modification of a Makefile, and a lot of headaches. I gave up on that. You can also clone the git repository and use the pre-compiled executable there.

    git clone https://github.com/lingtikong/phana.git

This will probably save you 3 steps out of 5. Personally, I prefer the use of Python and Jupyter notebooks since it allows more flexibility and insight to the raw data. 

## Inter-atomic potential

We are going to use the LJ potential for Argon. The parameters are from @whiteLennardJonesModelArgon1999. The values of $\sigma$ and $\epsilon$ are

| Parameter | Value       |
| --------- | ----------- |
| $\sigma$   | 3.419 Å |
| $\epsilon$ | 0.01015122 eV |

## Chemical properties of argon

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