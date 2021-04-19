# Interactive Phonon Dispersion Curve

Create an interactive dispersion curve for an FCC crystal using LAMMPS, Python, and Bokeh. Useful for studying the effect of pressure (or another parameter) on the normal modes.

## Requirements

- A valid LAMMPS installation with fix-phonon enabled
- Python
- Bokeh

## Output

The project creates a standalone web page that lets the user explore the change in the dispersion curve as a function of pressure

## Installation

    git clone https://github.com/jorox/dispapp.git
    cd disp_app
    mkdir workspace
    ./main.sh

## main.sh

The project source files are under ./src in two folders

- lammps
- python

The first holds the input file needed by LAMMPS to run the simulation. It only requires a data file and a map file to work. These can be generated using the python module developed for this task

## pylmp

A python wrapper to use the modules for pre and post processing

## mdlib.fixphonon

The fixphonon submodule has been tested using the CuPhonon data provided by LAMMPS

## Documentation

To build the documentation for the project

    cd doc
    make clean
    sphinx-apidoc -f -o source/ ../../src/python
    make html
