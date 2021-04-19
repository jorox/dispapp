# -*- coding: UTF-8 -*-.
"""Process data from log files of fix-phonon. 
    
    Use it to parse log files, find intermediate eigenvector-eigenvalue pairs, sort eigenvalues, and more.
    Cannot be used as a standalone module
    
    .. todo:: Add possibility to process log file to disp file directly
"""
import numpy as np
from itertools import permutations
import re
import numpy as np

def get_data(fin :str) -> tuple[np.array, np.array]:
    """Parse fix-phonon log file for eigenvectors and eigenvalues

    This function only looks for the last dynamical matrix in the log file. It sorts the eigenvalues for each q-vector by increasing order

    :param fin: log file name
    :type fin: str
    :return: numpy arrays of eigenvectors and eigenvalues unsorted
    :rtype: tuple[np.array, np.array]
    """
    def find_line_number(lst:list,key:int) -> int:
        """Find the *last* line starting with a substring

        Change this function to return the first or last dynamical matrix outputed by fix-phonon

        :param lst: a list of strings from file.readlines()
        :type lst: list
        :param key: the substring ket
        :type key: str
        :return: the index of the last line we are looking for
        :rtype: int
        """
        for i in range(len(lst)-1,0,-1):
            elem_str= lst[i]
            if elem_str.startswith(key): return i
        return -1

    lines = open(fin).readlines()
    # Fix the text for processing (whitespaces)
    for i in range(len(lines)):
        lines[i] = lines[i].strip()            # remove trailing white spaces
        lines[i] = re.sub('\t', ' ', lines[i]) # replace tab with white space
        lines[i] = re.sub(' +', ' ', lines[i]) # remove consec white spaces

    # Get the Dynamical matrix values
    ipStart = find_line_number(lines,"# qx qy qz") # find last Phi section
    ipEnd = len(lines)  
    phitext = lines[ipStart+1:ipEnd]
    nlines = len(phitext)
    nvals = (len(phitext[0].split())-3)//2

    # Create matrices to hold Phi, q, and eigenvalues
    phi = np.zeros((nlines,nvals), dtype=np.complex)
    q = np.zeros((nlines,3))
    eig = np.zeros((nlines,3))
    for i in range(nlines):
        tmp = np.fromstring(phitext[i], sep=" ")
        q[i,:] = tmp[:3]
        phi[i,:] = np.array( [complex(x,y) for 
            x,y in zip(tmp[3::2], tmp[4::2])] )
        eig[i,:],eigvec = np.linalg.eig(np.reshape(phi[i,:], (3,3)))
        eig[i,:] = np.sort(eig[i,:])
    return eig,q

def collinear(stpt: np.array, endpt: np.array, pt: np.array) -> bool:
    """Determine whether three points are collinear and in in order. Useful for finding intermediate eigenvectors

    Check if the three points are parallel and if the midpoint is between the start and end exclusively.

    :param stpt: The starting point
    :type stpt: np.array
    :param endpt: The end point
    :type endpt: np.array
    :param pt: The point to check
    :type pt: np.array
    :return: true or false
    :rtype: bool
    """
    a = endpt-stpt
    b = pt-stpt
    anorm = np.linalg.norm(a)
    if anorm < 0.0001: print("ERROR: anrom cannot be zero")
    bnorm = np.linalg.norm(b)
    ahat = a/anorm
    if bnorm>0.0000001: 
        bhat = b/bnorm 
    else: 
        bhat = b
    return np.abs(np.dot(ahat[0],bhat[0]) - 1) < 0.00001 and bnorm < anorm #same strtpoint or collinear


def sort_eigvals(qsord:np.array, eigenqs:np.array) -> np.array:
    """Sort eigenvalues by distance from previous point 

    The function tries to sort the eigenvalues for each q_i vector by 
    finding the closest eigenvalue to the previous e_j to e_i-1. 
    The algorithm does not always produce correct results, but is fast
    
    .. note:: 
        Assumes that the q-vectors have to be in increasing order

    .. note:: don't use this. Just sort eigenvalues by increasing
    
    :param qsord: nx3 array of eigenvectors in order
    :type qsord: np.array
    :param eigenqs: nx3 array of eigenvalues
    :type eigenqs: np.array
    :return: return an eigenvalue array now sorted 
    :rtype: np.array
    """
    neig = eigenqs.shape[1]
    eigenqsord = eigenqs
    for iq in range(1,qsord.shape[0]):
        x = np.linalg.norm(qsord[iq,:]-qsord[iq-1,:])
        for i in range(neig):
            y = eigenqs[iq,i]-eigenqs[iq-1,i]
            w = np.linalg.norm([x,y])
            for j in range(neig):
                y2 = eigenqs[iq,j]-eigenqs[iq-1,i]
                w2 = np.linalg.norm([x,y2])
                if w2<w:
                    #print("!!! SWAPPING !!!")
                    eigenqsord[iq,i], eigenqsord[iq,j] = eigenqsord[iq,j], eigenqsord[iq,i]
    return eigenqsord


def sort_eigvals_lsq(qsord:np.array, eigenqs:np.array) -> np.array:
    """ Sort eigenvalues by sum of least-squares method

    Calculate the distance between (q_i,e_i,j) and (q_i+1,e_i+1,k) for all
    permutations of k, and the find the permutation that produces the least
    sum of squares of the distances j-k.
    
    :param qsord: nx3 array of eigenvectors in order
    :type qsord: np.array
    :param eigenqs: nx3 array of eigenvalues
    :type eigenqs: np.array
    :return: return an eigenvalue array now sorted 
    :rtype: np.array

    .. note:: Assumes that the q vectors are in inceasing order
    """
    neigs = eigenqs.shape[1]
    perm = list(permutations(range(neigs)))
    eigenqsord = eigenqs
    for iq in range(1,qsord.shape[0]):
        dstmin = np.linalg.norm( 
            (qsord[iq,:]-qsord[iq-1,:], eigenqs[iq,perm[0]]-eigenqs[iq-1,perm[0]])
        )
        ipmin = 0
        for ip in range(1,len(perm)):
            dst = np.linalg.norm( 
            (qsord[iq,:]-qsord[iq-1,:], eigenqs[iq,perm[ip]]-eigenqs[iq-1,perm[0]])
            )
            if dst < dstmin: 
                ipmin = ip
                dstmin = dst 
        #print(iq,perm[ipmin])
        eigenqsord[iq,:] = eigenqs[iq,perm[ipmin]]
    
    return eigenqsord


def find_eigvals(qs:np.array, q:np.array, D:np.array) -> np.array or np.NaN :
    """Find the eigenvalues for a specific eigenvector

    Return the eigenvalues for a specific vector. If the vector is not in the list return np.NaN

    :param qs: mx3 array of eigenvector(s) we are looking for
    :type qs: np.array
    :param q: nx3 array of eigenvectors
    :type q: np.array
    :param D: nx3 array of eigenvectors
    :type D: np.array
    :return: return the eigenvalues or NaN if not found
    :rtype: np.array or np.NaN
    """
    eigsqs = np.zeros((qs.shape[0], D.shape[1]))
    for iqs in range(qs.shape[0]):
        found = False
        for iq in range(q.shape[0]):
            if np.linalg.norm(qs[iqs]-q[iq]) < 0.000001:
                eigsqs[iqs] = D[iq]
                found = True
        if not found: eigsqs[iqs] = np.NaN
    return eigsqs


def find_intermediate_q(qstart:np.array, qend:np.array, qarray:np.array) -> list:
    """Return intermediate q's in increasing order between qstart and qend with their distance from qstart

    From a list of eigenvectors find the list which are between two points (exclusively)

    :param qstart: starting q-vector
    :type qstart: np.array
    :param qend: ending q-vectors
    :type qend: np.array
    :param qarray: nx3 array of all q-vectors
    :type qarray: np.array
    :return: list of q-vectors which are in between
    :rtype: list
    """
    qint = []
    for q in qarray:
        if collinear(qstart, qend, q):
            # add q and its distance to the list
            qint.append( q )

    res = np.array(  sorted( qint, key=lambda qp : np.linalg.norm(qp-qstart) ) ) #sort according to distance from qstart
    return res


def build_curve(SPs:list, eig:np.array, q:np.array) -> tuple[np.array, np.array, list]:
    """Build a specific dispersion curve using q vectors and eigen values
    .. todo:: make return a single np.array and adjust phana.ipynb

    Main function that builds a dispersion curve which is the eigenvalues and eigenvectors along a trajectory of special points
    .. note:: the function does not check if the special points form a valid line

    :param SPs: List of special points to find the dispersion curve along
    :type SPs: list
    :param eig: nx3 array of eigenvalues
    :type eig: np.array
    :return: an n array representing the q-points of the trajecotry, the eigenvalues/vectors, and the q-points of the special points
    :rtype: tuple[np.array, np.array, list]
    """
    # Initialize data structures
    neigvals   = len(eig[0])
    
    disp_data  = []  # dispersion data - <nlines

    iSP = [0]
    for i in range(0,len(SPs)-1):
        qstart = np.reshape(SPs[i], (1,3))     # Start point
        qend   = np.reshape(SPs[i+1], (1,3))
        qint   = find_intermediate_q(qstart,qend,q) # q's are sorted
        if (qint.shape[0] > 0): qstint = np.append(qstart, qint, axis=0)
        else: qstint = qstart

        eigqstint = find_eigvals(qstint, q, eig)

        for iq in range(qstint.shape[0]):
            if np.any(np.isnan(eigqstint[iq,0])): continue
            disp_data.append(qstint[iq].tolist())
            for ie in range(neigvals):
                disp_data[-1].append( eigqstint[iq,ie] )
        iSP.append(len(disp_data))


    tmp = np.reshape(SPs[-1], (1,3))
    tmp2 = find_eigvals(tmp,q,eig).tolist()
    tmp3 = [qx for qx in SPs[-1]]
    for ie in tmp2: tmp3.extend(ie)
    if np.any(np.isnan(tmp2)): disp_data.append( tmp3 )


    disp_data = np.array(disp_data)
    #disp_data[:,3:] = sort_eigvals_lsq(disp_data[:,:3], disp_data[:,3:])

    qr = np.zeros((disp_data.shape[0]))
    for iq in range(1,qr.shape[0]):
        qr[iq] = np.linalg.norm(disp_data[iq,:3] - disp_data[iq-1,:3]) +qr[iq-1]
    
    SPqr=[qr[i] for i in iSP[:-1]] 

    return qr, disp_data, SPqr