"""

This script reads a structure file (pdb or gro file) and a gromacs trajectory file
(.xtc, is optional), and can either calculate the distances for a definded cutoff
between all the defined atoms of the helix or calculates the contact statistics
used then to build contact maps.

To build a contact maps :
python3 calc_resi-dist.py structure.gro -xtc trajectory.xtc -a -d 7 -bh1 5 -eh1 25 -bh2 36 -eh2 56 -a

"""
#!/usr/bin/env python

import MDAnalysis as mda
import MDAnalysis.analysis.distances
import numpy as np
import argparse
from math import sqrt


def get_arg():
    """Get arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("gro", help=".gro file")
    parser.add_argument("-xtc", help=".xtc file, if \"none\" only the gro file is "+
                        "taken into account", default="none")
    parser.add_argument("-d", help="distance cutoff", type=int, default=1000000)
    parser.add_argument("-time_min", help="minimal time cutoff", type=int,
                        default=0)
    parser.add_argument("-time_max", help="maximum time cutoff, the default value is at"+
                        "10 microsec", type=int, default=10000000)
    parser.add_argument("-p", help="prints all the distances", action="store_true")
    parser.add_argument("-a", help="make contact analysis", action="store_true")
    parser.add_argument("-ele", help="elements checked for distances, bb id for backbone atoms"+
                        "all is for all the atoms", choices=["bb", "all"], default="bb")
    parser.add_argument("-bh1", type=int,
                        help="Number of the first residue for the first helix")
    parser.add_argument("-eh1", type=int,
                        help="Number of the last residue for the first helix")
    parser.add_argument("-bh2", type=int,
                        help="Number of the first residue for the second helix")
    parser.add_argument("-eh2", type=int,
                        help="Number of the last residue for the second helix")
    args = parser.parse_args()

    return(args)


def sort_key(element):
    """Define the condition to sort the list"""
    return element[1]


def contact_analysis(list_couples, total_time, args):
    """Count the number of uccurence of a couple"""
    # Make a list with all the possible combinations
    all_couples = [[i,j] for i in range(args.bh1, args.eh1+1)
                   for j in range(args.bh2, args.eh2+1)]
    # print(all_couples)

    # Count all the occurences of each combinaition [[couple1 couple2], occurences]
    contacts_statistics = []
    for couple in all_couples:
        if couple in list_couples:
            count = list_couples.count(couple)
            contacts_statistics.append([couple, count])

    # Trier la liste en fonction du nombre d'occurences
    contacts_statistics.sort(key=sort_key, reverse=True)

    print("atom1 atom2  count")
    nb_contact = len(contacts_statistics)
    to_be_printed = '\n'.join("{:5d} {:5d}  {:5f}".
                              format(contacts_statistics[i][0][0],
                                     contacts_statistics[i][0][1],
                                     float(contacts_statistics[i][1])/total_time)
                              for i in range(nb_contact))
    print(to_be_printed)


def get_dist_from_traj(univ, pept1, pept2, args):
    """Prints all the distances between all atoms of the backbone.
    The calculation starts at 1 micro sec"""

    # Parameters
    dist_cutoff = args.d
    pr = args.p
    ana = args.a
    time_min = int(args.time_min*len(univ.trajectory)//univ.trajectory[-1].time)
    time_max = int(args.time_max*len(univ.trajectory)//univ.trajectory[-1].time)
    if args.time_max > univ.trajectory[-1].time:
        time_max = len(univ.trajectory)
    print(time_min, time_max, len(univ.trajectory[time_min:time_max]))
    couples = []
    total_time = 0

    if pr:
        print("time   dist   atom1   atom2")

    # For each trajectory, calculates all the distances
    for ts in univ.trajectory[time_min:time_max]:
        matrix_dist = MDAnalysis.analysis.distances.capped_distance(pept1.positions,
                                                                    pept2.positions,
                                                                    dist_cutoff)
        matrix_length = len(matrix_dist[1])

        # Print for each distances one line with the time, distance, atom1, atom2
        if pr:
            to_be_printed = '\n'.join("{:8.0f}   {:.3f}   {}   {}".
                                      format(univ.trajectory.time,
                                             matrix_dist[1][i],
                                             matrix_dist[0][i][0]+1,
                                             matrix_dist[0][i][1]+1)
                                      for i in range(matrix_length))
            print(to_be_printed)

        if ana:
            to_be_added = [[pept1[list(ele)[0]].resid, pept2[list(ele)[1]].resid]
                           for ele in matrix_dist[0]]
            to_be_added_n = []
            [to_be_added_n.append(x) for x in to_be_added if x not in to_be_added_n]
            couples += to_be_added_n
            total_time += 1

    if ana:
        # Call to analysis function
        contact_analysis(couples, total_time, args)


def main():
    args = get_arg()

    if args.xtc == "none" :
        univ = mda.Universe(args.gro)
    else:
        univ = mda.Universe(args.gro, args.xtc)

    if args.ele == "bb":
        # The first peptide has 29 residues
        pept1_command = "resid {0}-{1} and name N CA C O".format(args.bh1, args.eh1)
        pept1 = univ.select_atoms(pept1_command)
        # The second peptide has 29 residues
        pept2_command = "resid {0}-{1} and name N CA C O".format(args.bh2, args.eh2)
        pept2 = univ.select_atoms(pept2_command)

    if args.ele == "all":
        # The first peptide has 29 residues
        pept1_command = "resid {0}-{1}".format(args.bh1, args.eh1)
        pept1 = univ.select_atoms(pept1_command)
        # The second peptide has 29 residues
        pept2_command = "resid {0}-{1}".format(args.bh2, args.eh2)
        pept2 = univ.select_atoms(pept2_command)

    # Prints the distances lower than the cutoff (if -p)
    # And/Or prints the contacts statistics (if -a)
    get_dist_from_traj(univ, pept1, pept2, args)


if __name__=="__main__":
    main()
