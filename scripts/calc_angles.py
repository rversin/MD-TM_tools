"""

This script reads a structure file (pdb or gro file) and a gromacs trajectory file
(.xtc, is optional), and calculate either the crossing-angle between two helices,
or the tilt angle.
The results will be shown in the shell.

Get the help :
python3 calc_angles.py -h

Usage for crossing angle:
1) for a whole trajectory
python3 calc_angles.py structure.gro -xtc trajectory.xtc -bh1 1 -eh1 20 -bh2 21 -eh2 42 > CA.dat
2) for only one structure
python3 calc_angles.py structure.gro -bh1 1 -eh1 20 -bh2 21 -eh2 42 > CA.dat

Usage for tilt angle:
1) for a whole trajectory
python3 calc_angles.py structure.gro -xtc trajectory.xtc -bh1 1 -eh1 20 -mb > TILT.dat
python3 calc_angles.py structure.gro -xtc trajectory.xtc -bh1 21 -eh1 42 -mb -z -20 > TILT.dat
2) for only one structure
python3 calc_angles.py structure.gro -bh1 1 -eh1 20 -mb > TILT.dat
python3 calc_angles.py structure.gro -bh1 21 -eh1 42 -mb -z -20 > TILT.dat
WARNING : the tilt angle calculation assumes that the bilayer normal is aligned along Z.

"""
#!/usr/bin/env python

import MDAnalysis as mda
import numpy as np
import argparse
import math
import sys


RES_BEGIN = 0
RES_END = 50

def get_arg():
    """Get arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("gro", help=".gro file")
    parser.add_argument("-xtc", help=".xtc file, if \"none\" only the gro file is "+
                        "taken into account", default="none")
    parser.add_argument("-mb", help="If in command line, calculates the angle of the"+
                        " helix with the membrane. If not, calculates the crossing angle"+
                        " between the two helix.", action="store_true")
    parser.add_argument("-bh1", type=int,
                        help="Number of the first residue for the first helix")
    parser.add_argument("-eh1", type=int,
                        help="Number of the last residue for the first helix")
    parser.add_argument("-bh2", type=int, default=0,
                        help="Number of the first residue for the second helix")
    parser.add_argument("-eh2", type=int, default=0,
                        help="Number of the last residue for the second helix")
    parser.add_argument("-z", type=int, default=20,
                        help="Define the Z component of the axis normal to the membrane")
    args = parser.parse_args()
    return(args)


def calc_centroid(list_coord):
    nb_atom = len(list_coord)
    sum_x = np.sum([ele[0] for ele in list_coord])
    sum_y = np.sum([ele[1] for ele in list_coord])
    sum_z = np.sum([ele[2] for ele in list_coord])
    return [sum_x/nb_atom, sum_y/nb_atom, sum_z/nb_atom]


def get_vector_pept(pept, second_pept=False):
    """For one helical peptide get the main vector"""
    turn_helix_1 = [pept[RES_BEGIN].position, pept[RES_BEGIN+1].position,
                    pept[RES_BEGIN+2].position, pept[RES_BEGIN+3].position]
    centroid_1 = calc_centroid(turn_helix_1)
    turn_helix_2 = [pept[RES_END-3].position, pept[RES_END-2].position,
                    pept[RES_END-1].position, pept[RES_END].position]
    centroid_2 = calc_centroid(turn_helix_2)
    return([centroid_1[0]-centroid_2[0],
           centroid_1[1]-centroid_2[1],
           centroid_1[2]-centroid_2[2]])


def get_vector_angle(v1, v2):
    """Simple (wrong?) calculation of angles between two vector"""
    c = np.dot(v1,v2)/np.linalg.norm(v1)/np.linalg.norm(v2)
    angle = np.arccos(np.clip(c, -1, 1))
    return math.degrees(angle)


def scalar(A,B):
    """Returns the scalar (or inner) product between vectors A & B
    """
    return (A[0]*B[0]) + (A[1]*B[1]) + (A[2]*B[2])


def get_crossing_angle(b_pept1, t_pept1, t_pept2, b_pept2):
    """Calculates the crossing angle between two helices (dihedral angle)"""
    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4 = b_pept1+t_pept1+t_pept2+b_pept2
    ab_x = (x2 - x1)
    ab_y = (y2 - y1)
    ab_z = (z2 - z1)
    bc_x = (x3 - x2)
    bc_y = (y3 - y2)
    bc_z = (z3 - z2)
    cd_x = (x4 - x3)
    cd_y = (y4 - y3)
    cd_z = (z4 - z3)	
    d012  = ab_x * bc_x + ab_y * bc_y + ab_z * bc_z
    d123  = cd_x * bc_x + cd_y * bc_y + cd_z * bc_z
    d0123 = ab_x * cd_x + ab_y * cd_y + ab_z * cd_z
    d01   = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z
    d12   = bc_x * bc_x + bc_y * bc_y + bc_z * bc_z
    d23   = cd_x * cd_x + cd_y * cd_y + cd_z * cd_z
    num = d012 * d123 - d12 * d0123
    den = (d01*d12 - d012*d012)*(d12*d23 - d123*d123)
    arccos = num / math.sqrt(den)
    if arccos > 1.:
        arccos = 1.
    if arccos < -1.:
        arccos = -1.
    RS = math.acos(arccos)
    RS1 = cd_x * (ab_y * bc_z - ab_z * bc_y) + \
          cd_y * (bc_x * ab_z - ab_x * bc_z) + \
          cd_z * (ab_x * bc_y - ab_y * bc_x)
    if RS1 > 0:
        return math.degrees(RS)
    else:
        return -math.degrees(RS)


def get_inertia_axis(COOR,eigen_index=0):
    """This function computes the 3 inertia axis of a set of atoms (in COOR format) and
    returns the chosen axis (by default, the first one). The object returned is of type array.
    """
    # convert COOR to a 3rows*Ncol rows matrix (array object)
    X = np.transpose(np.array(COOR))
    # scale X around the origin
    X[0] = X[0] - np.mean(X[0])
    X[1] = X[1] - np.mean(X[1])
    X[2] = X[2] - np.mean(X[2])
    # get the inertia matrix (matrix product of X time its transpose)
    Y = np.dot(X, np.transpose(X))
    # calculate the eigen values and vectors
    eigenvalues = np.linalg.eig(Y)[0]
    eigenvectors = np.linalg.eig(Y)[1]
    # sort out eigenvalues
    index_eigenvalues_sorted = list(np.argsort(eigenvalues))
    # reverse order cause argsort sorts out in crescent order
    index_eigenvalues_sorted.reverse()
    # return the wanted eigenvector (generally the one of greatest eigen value)
    return eigenvectors[:,index_eigenvalues_sorted[eigen_index]]


def get_helix_axis(pept, scaling=20):
    """This function takes a list of atoms in COOR format, 2 integer indices (to indicate
    which is the first residue and the last residue of a given helix) and returns the helix
    axis as a vector (in array format). Basically, the vector is oriented from Nterm to Cterm.
    The scaling is by default 20, which gives an arrow slightly longer than a helix of 23 residues.
    """
    # get COOR from helix
    HCOOR=[]
    maxim = RES_END+1
    if RES_END+1 > len(pept):
        maxim = len(pept)
    for i in range(RES_BEGIN, maxim):
        HCOOR.append(pept[i].position)

    # compute the helix axis
    Haxis = get_inertia_axis(HCOOR)
    # get center of mass
    HCOM = pept[RES_BEGIN:RES_END+1].center_of_mass()
    # get the vector COM-lastCA
    vCOM_Cterm = [HCOOR[-1][0]-HCOM[0],
                  HCOOR[-1][1]-HCOM[1],
                  HCOOR[-1][2]-HCOM[2]] #vect_AB(HCOM,HCOOR[-1])
    # orient the axis from N to C term
    # see http://en.wikipedia.org/wiki/Scalar_resolute for scalar resolute (projection of a vector onto another)
    if get_vector_angle(Haxis, vCOM_Cterm) > get_vector_angle(-Haxis, vCOM_Cterm):
    #angle_vector(Haxis,vCOM_Cterm) > angle_vector(-Haxis,vCOM_Cterm):
        Haxis = -Haxis
    return Haxis


def get_angle_pept(univ, pept1, pept2):
    """Get all the angles between two helices, for a trajectory"""
    # Print columns names
    print("time angle")

    # For each 
    for ts in univ.trajectory:

        # Calculates axis of the two helices
        # axis_pept1 = get_vector_pept(pept1)
        # axis_pept2 = get_vector_pept(pept2)
        axis_pept1 = get_helix_axis(pept1)
        axis_pept2 = get_helix_axis(pept2)

        # Calculates the COM of the two helices
        cm_pept1 = pept1[RES_BEGIN:RES_END+1].center_of_mass()
        cm_pept2 = pept2[RES_BEGIN:RES_END+1].center_of_mass()

        maxim = RES_END
        if RES_END+1 > len(pept1) or RES_END+1 > len(pept2):
            maxim = -1

        # Calculates beginning and ending of the helices
        scal_b_pept1 = scalar(axis_pept1, pept1[RES_BEGIN].position-cm_pept1)
        scal_b_pept2 = scalar(axis_pept2, pept2[RES_BEGIN].position-cm_pept2)
        scal_e_pept1 = scalar(axis_pept1, pept1[maxim].position-cm_pept1)
        scal_e_pept2 = scalar(axis_pept2, pept2[maxim].position-cm_pept2)
        b_pept1 = [cm_pept1[i] + scal_b_pept1*axis_pept1[i] for i in range(3)]
        b_pept2 = [cm_pept2[i] + scal_b_pept2*axis_pept2[i] for i in range(3)]
        e_pept1 = [cm_pept1[i] + scal_e_pept1*axis_pept1[i] for i in range(3)]
        e_pept2 = [cm_pept2[i] + scal_e_pept2*axis_pept2[i] for i in range(3)]

        # Find the smallest distance between both helices
        w1 = scalar([b_pept1[i] - b_pept2[i] for i in range(3)],
                    [e_pept1[i] - b_pept1[i] for i in range(3)])
        w2 = scalar([b_pept1[i] - b_pept2[i] for i in range(3)],
                    [e_pept2[i] - b_pept2[i] for i in range(3)])
        u11 = scalar([e_pept1[i] - b_pept1[i] for i in range(3)],
                     [e_pept1[i] - b_pept1[i] for i in range(3)])
        u12 = scalar([e_pept1[i] - b_pept1[i] for i in range(3)],
                     [e_pept2[i] - b_pept2[i] for i in range(3)])
        u22 = scalar([e_pept2[i] - b_pept2[i] for i in range(3)],
                     [e_pept2[i] - b_pept2[i] for i in range(3)])

        det = u11*u22 - u12**2
        if det:
            s_pept1 = (-w1*u22 + w2*u12) / det
            s_pept2 = (w2*u11 - w1*u12) / det
        else:
            s_pept1, s_pept2 = 0.5, 0.5
        if s_pept1 > 1:
            s_pept1 = 1
        elif s_pept1 < 0:
            s_pept1 = 0.01 # <--- ??? should be zero but not correct to calculate the dihedral !!!
        if s_pept2 > 1:
            s_pept2 = 1
        elif s_pept2 < 0:
            s_pept2 = 0.01   # <--- ??? should be zero but not correct to calculate the dihedral !!!
        # calculate tk and tj
        t_pept1 = [b_pept1[k] + s_pept1*(e_pept1[k] - b_pept1[k]) for k in range(3)]
        t_pept2 = [b_pept2[k] + s_pept2*(e_pept2[k] - b_pept2[k]) for k in range(3)]

        # print("{:.0f}  {:.3f}".format(univ.trajectory.time, get_vector_angle(axis_pept1,
        #                                                                      axis_pept2)))
        # print("b_pept1 :", b_pept1, "b_pept2 :", b_pept2, "t_pept1", t_pept1, "t_pept2 :", t_pept2)
        print("{:.0f}  {:.3f}".format(univ.trajectory.time,
                                      get_crossing_angle(b_pept1, t_pept1, t_pept2, b_pept2)))


def get_angle_mb(univ, pept, Z):
    """Get all the angles between the normal of the membrane and the two helix,
    for a trajectory"""
    print("time angle")
    norm_mb = [0, 0, Z]
    for ts in univ.trajectory:
        axis_pept = get_helix_axis(pept)
        print("{:.0f}  {:.3f}".format(univ.trajectory.time,
                                         get_vector_angle(axis_pept, norm_mb)))


def main():
    # Get arguments
    args = get_arg()

    # Check if an xtc file was given
    # If not given, only use the structure file
    if args.xtc=="none":
        univ = mda.Universe(args.gro)
    # If the trajectory was given, use it
    else:
        univ = mda.Universe(args.gro, args.xtc)

    # Check if the tilt is asked
    if args.mb:
        # If an other helix is given
        if args.bh2 or args.eh2:
            sys.exit('Only one helix is expected for the tilt')

        pept1_command = "resid {0}-{1} and name CA".format(args.bh1, args.eh1)
        pept1 = univ.select_atoms(pept1_command)

        get_angle_mb(univ, pept1, args.z)

    # If the tilt is not asked, calculate the crossing-angle
    else:
        # If an helix is missing
        if not args.bh2 and not args.eh2:
            sys.exit('Two helix are expected for the crossing-angle')

        # Define the first peptide
        pept1_command = "resid {0}-{1} and name CA".format(args.bh1, args.eh1)
        pept1 = univ.select_atoms(pept1_command)

        # Define the first peptide
        pept2_command = "resid {0}-{1} and name CA".format(args.bh2, args.eh2)
        pept2 = univ.select_atoms(pept2_command)

        get_angle_pept(univ, pept1, pept2)


if __name__=="__main__":
    main()

