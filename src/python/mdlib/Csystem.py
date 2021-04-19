import numpy as np
class Csystem:
    """  Represent a cyrstalline system: basis vectors + basis atoms
        
    :param cell: 3x3 array column vectors are the three basis vectors
    :type cell: np.array
    :param atoms: 3xn array, column vectors are positions of basis atoms
    :type atoms: np.array
    :param origin: a 3x1 array for the origin of the unit cell
    :type origin: np.array  
    :param natoms: Number of basis atoms
    :type natoms: int
    :param ntypes: Number of different atom types 1 - natoms
    :type ntypes: int
    :attr masses: A dictionay for the mass of each atom type
    :type masses: dict

    .. todo:: add a method to transform the basis vectors for the crystalline system

    """

    def __init__(self, acell:np.array, aorigin:np.array, abasis:np.array, amasses=None):
        """Constructor function for the class

        Initialize the main parameters

        :param acell: the basis vectors
        :type acell: numpy.array
        :param aorigin: the origin of the system
        :type aorigin: numpy.array
        :param abasis: the basis atoms 3xn array
        :type abasis: numpy.array
        :param amasses: the masses for the different types, defaults to None
        :type amasses: dict, optional
        """
        self.cell = acell
        self.atoms = abasis
        self.origin = aorigin
        self.natoms = abasis.shape[0]
        self.ntypes = int(max(abasis[:,1]))
        self.masses = dict(enumerate(range(self.ntypes)))
        if amasses is not None:
            for i in amasses.keys():
                self.masses[i] = amasses[i]

    def __str__(self) -> str:
        """Return string representation of the object

        Returns the cell vectors, atoms, and the origin

        :return: a string for quick info on the system
        :rtype: str
        """
        return ("{:d} atoms, {:d} atom types\n".format(self.natoms, self.ntypes)+
        "\n cell = \n" + str(self.cell) + 
        "\n origin = \n" + str(self.origin) +
        "\n atoms = \n" +str(self.atoms))
