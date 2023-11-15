""" R@ph!! """

import MDAnalysis as mda
import numpy as np
import argparse
import os
import math
from numba import njit

PHOSPHATE_ATOM_CG="PO4"
LIPID_TAILS_CG=["C1A", "D2A", "C3A", "C4A", "C1B", "C2B", "C3B", "C4B"]
PHOSPHATE_ATOM_AA="P"
LIPID_TAILS_AA=['C23', 'H3R', 'H3S', 'C24', 'H4R', 'H4S', 'C25', 'H5R',
                'H5S', 'C26', 'H6R', 'H6S', 'C27', 'H7R', 'H7S', 'C28',
                'H8R', 'H8S', 'C29', 'H91', 'C210', 'H101', 'C211', 'H11R',
                'H11S', 'C212', 'H12R', 'H12S', 'C213', 'H13R', 'H13S',
                'C214', 'H14R', 'H14S', 'C215', 'H15R', 'H15S', 'C216',
                'H16R', 'H16S', 'C217', 'H17R', 'H17S', 'C218', 'H18R',
                'H18S', 'H18T', 'C33', 'H3X', 'H3Y', 'C34', 'H4X', 'H4Y',
                'C35', 'H5X', 'H5Y', 'C36', 'H6X', 'H6Y', 'C37', 'H7X',
                'H7Y', 'C38', 'H8X', 'H8Y', 'C39', 'H9X', 'H9Y', 'C310',
                'H10X', 'H10Y', 'C311', 'H11X', 'H11Y', 'C312', 'H12X',
                'H12Y', 'C313', 'H13X', 'H13Y', 'C314', 'H14X', 'H14Y',
                'C315', 'H15X', 'H15Y', 'C316', 'H16X', 'H16Y', 'H16Z']


def get_arg():
    """Get arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-xtc", help="trajectory")
    parser.add_argument("-gro", help="structure")
    parser.add_argument("-time_min", help="minimal time cutoff", type=int,
                        default=0)
    parser.add_argument("-time_max", help="maximum time cutoff", type=int,
                        default=10000000)
    parser.add_argument("-c", help="protrusion cutoff, in nm", type=float,
                        default=0.1)
    parser.add_argument("-d", help="calculates protrusion as a function of distances"
                        " to the protein",
                        action="store_true")
    parser.add_argument("-r", help="representation: aa (charmm36) or cg (martini2/3)",
                        choices=["aa", "cg"], default="aa")
    args = parser.parse_args()

    return(args)


def get_middle_of_mb(select_phosphate):
    """Finds the middle of the membrane on the Z axis"""
    return(np.average(np.array(select_phosphate)))


@njit
def dot_product(vect1, vect2):
    """Dot productor of 3D vectors"""
    return(vect1[0]*vect2[0]+vect1[1]*vect2[1]+vect1[2]*vect2[2])


def get_pbc_dist(center, x, box_x, box_y):
    """Calculates distance from x according to minimum image convention"""
    dx = x[0]-center[0]
    if dx > box_x/2:
        dx -= box_x
    elif dx <= -box_x/2:
        dx += box_x

    dy = x[1]-center[1]
    if dy > box_y/2:
        dy -= box_y
    elif dy <= -box_y/2:
        dy += box_y

    return(math.sqrt(dx*dx + dy*dy))


def compute_protrusion_simple(middle_mb, lipid_atoms, cutoff, tot_lipid, size_lipid, time):
    """computes for a selection the protrusion"""
    nb_of_protrusion = 0

    for i in range(tot_lipid):
        phosphate = lipid_atoms[0][i]
        
        # Check if the lipid is in the lower layer
        if phosphate[2] < middle_mb:
            normal_mb = [0, 0, -1]
        # is in the upper layer
        else:
            normal_mb = [0, 0, 1]

        is_protrusion = False
        for tail in lipid_atoms[1][i*size_lipid:(i+1)*size_lipid]:
            vect_lipid_phos = tail - phosphate
            # calculates the projection of vect_lipid_phos on normal_mb
            # this should be dot(normal_mb, vect_lipid_phos)/norm(normal_mb)
            # norm(normal_mb) is equal to 1
            projection = np.dot(normal_mb, vect_lipid_phos)

            if projection >= cutoff:
                is_protrusion = True

        if is_protrusion:
            nb_of_protrusion += 1

    print("{:8.0f}   {:.5f}".
          format(time, nb_of_protrusion/tot_lipid))
    # return(nb_of_protrusion)


def compute_protrusion_distance(middle_mb, lipid_atoms, cutoff, tot_lipid, size_lipid, time):
    """computes for a selection the protrusion"""
    nb_of_protrusion = 0
    protrusion_for_dist = [[i, 0, 0] for i in range(0, int((lipid_atoms[3]/2)+10),10)]

    for i in range(tot_lipid):
        phosphate = lipid_atoms[0][i]

        # Check if the lipid is in the lower layer
        if phosphate[2] < middle_mb:
            normal_mb = [0, 0, -1]
        # is in the upper layer
        else:
            normal_mb = [0, 0, 1]

        dist_phos = get_pbc_dist(lipid_atoms[2], phosphate, lipid_atoms[3], lipid_atoms[4])
        index_lipid = int(dist_phos//10)
        protrusion_for_dist[index_lipid][2] += 1

        is_protrusion = False
        for tail in lipid_atoms[1][i*size_lipid:(i+1)*size_lipid]:
            vect_lipid_phos = tail - phosphate
            # calculates the projection of vect_lipid_phos on normal_mb
            # this should be dot(normal_mb, vect_lipid_phos)/norm(normal_mb)
            # norm(normal_mb) is equal to 1
            projection = np.dot(normal_mb, vect_lipid_phos)

            if projection >= cutoff:
                is_protrusion = True
                break

        if is_protrusion:
            protrusion_for_dist[index_lipid][1] += 1

    if protrusion_for_dist[-1][2] == 0:
        protrusion_for_dist = protrusion_for_dist[:-1]
    protrusion_for_dist_final = [[protru_slice[0], protru_slice[1]/protru_slice[2]]
                                 for protru_slice in protrusion_for_dist]

    to_be_printed = '\n'.join("{:8.0f}   {}   {:.5f}".
                              format(time, protru_slice[0], protru_slice[1])
                                     for protru_slice in protrusion_for_dist_final)
    print(to_be_printed)


def protrusion(univ, middle_mb, lipid_atoms, args, size_lipid):
    """Calculates the protrusion for the whole bilayer as a function of time."""
    time_min = args.time_min
    time_max = args.time_max
    cutoff = args.c*10
    tot_lipid = len(lipid_atoms[0])

    if not args.d:
        print("time prop_protrusion")
        compute_protrusion = compute_protrusion_simple
    elif args.d:
        print("time distance prop_protrusion")
        compute_protrusion = compute_protrusion_distance
    
    for ts in univ.trajectory:
        if univ.trajectory.time >= time_min and univ.trajectory.time <= time_max:
            lipid_atoms_positions = [lipid_atoms[0].positions,
                                     lipid_atoms[1].positions,
                                     lipid_atoms[2].center_of_mass(),
                                     univ.dimensions[0],
                                     univ.dimensions[1]]
            nb_of_protrusion = compute_protrusion(middle_mb,
                                                  lipid_atoms_positions,
                                                  cutoff,
                                                  tot_lipid,
                                                  size_lipid,
                                                  univ.trajectory.time)


def main():
    args = get_arg()
    univ = mda.Universe(args.gro, args.xtc)

    if args.r == 'aa':
        select_phosphate = "name P"#"resid 0-256 and name P"
        select_lipid_tail = " or ".join("name {}".format(lipid_atom)
                                        for lipid_atom in LIPID_TAILS_AA)
        # select_lipid_tail = "resid 0-256 and (" + select_lipid_tail + ")"
        select_protein = "name CA N C O"
        lipid_atoms = [univ.select_atoms(select_phosphate),
                       univ.select_atoms(select_lipid_tail),
                       univ.select_atoms(select_protein)]
        size_lipid = len(LIPID_TAILS_AA)
    elif args.r == 'cg':
        select_phosphate = "name PO4"#"resid 30-229 and name PO4"
        select_lipid_tail = " or ".join("name {}".format(lipid_atom)
                                        for lipid_atom in LIPID_TAILS_CG)
        # select_lipid_tail = "resid 30-229 and (" + select_lipid_tail + ")"
        select_protein = "name BB"
        lipid_atoms = [univ.select_atoms(select_phosphate),
                       univ.select_atoms(select_lipid_tail),
                       univ.select_atoms(select_protein)]
        size_lipid = len(LIPID_TAILS_CG)

    middle_mb = get_middle_of_mb(lipid_atoms[0].positions[:,2])
    # print(len(lipid_atoms[0]), len(lipid_atoms[1]), len(lipid_atoms[2]))
    
    # if not args.d:
    # print("time prop_protrusion")
    protrusion(univ, middle_mb, lipid_atoms, args, size_lipid)
    # elif args.d:
    #     print("time distance prop_protrusion")
    #     protrusion_distance(univ, args, select_phosphate)


if __name__=="__main__":
    main()
