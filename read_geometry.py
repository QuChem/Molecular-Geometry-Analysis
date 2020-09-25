# Title:	read_geometry.py
# Author:	Reza Hemmati
# Created	09/25/2020
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#

# GEOMETRY ANALYSIS
# Usage: python3  read_geometry.py  xyz_file

import numpy as np
import math, sys
import argparse
from sympy import *

# X, X1 are dummy atoms.
atomic_masses = {'H' : 1.00794, 'He': 4.00260, 'Li': 6.94100, 'Be': 9.01218, 'B' : 10.8110,
		'C' : 12.0107, 'N' : 14.0067, 'O' : 15.9994, 'F' : 18.9984, 'Ne': 20.1797,
		'Na': 22.9898, 'Mg': 24.3050, 'Al': 26.9815, 'Si': 28.0855, 'P' : 30.9738,
		'S' : 32.0650, 'Cl': 35.4530, 'Ar': 39.9480, 'K' : 39.0983, 'Ca': 40.0780,
      		'Sc': 44.9559, 'Ti': 47.8670, 'V' : 50.9415, 'Cr': 51.9961, 'Mn': 54.9380,
	       	'Fe': 55.8450, 'Co': 58.9332, 'Ni': 58.6934, 'Cu': 63.5460, 'Zn': 65.4090,
		'X' : 0.0, 'X1' : 0.0}


# Atomic number of elements.
atomic_znumber = {'H' : 1, 'He': 2, 'Li': 3, 'Be': 4, 'B' : 5,
        	'C' : 6, 'N' : 7, 'O' : 8, 'F' : 9, 'Ne': 10,
          	'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P' : 15,
          	'S' : 16, 'Cl': 17, 'Ar': 18, 'K' : 19, 'Ca': 20,
          	'Sc': 21, 'Ti': 22, 'V' : 23, 'Cr': 24, 'Mn': 25,
          	'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
          	'X' : 0.0, 'X1': 0.0}



def read_xyz(filename):
	"""
	Read a file containing the geometry of molecules.
	"""
	
	xyz_file = open(filename, 'r')
	
	if not xyz_file.closed:
		
		# Read the first line
		natom = int(xyz_file.readline())
		atom_symboles = []
		
		# Make Nx3 matrix of coordinates
		xyz_arr = np.zeros([natom, 3])
		
		i = 0
		
		for line in xyz_file:
			words = line.split()
			if len(words) > 3:
				atom_symboles.append(words[0])
				xyz_arr[i][0] = float(words[1])
				xyz_arr[i][1] = float(words[2])
				xyz_arr[i][2] = float(words[3])
				i += 1
	return atom_symboles, xyz_arr


def print_geom(atom_num, xyz_coords):
	
	print (F'Number of atoms: ', atom_num)
	print (F'Input Cartesian coordinates:')
	
	for i in range(atom_num):
		for j in range(3):
			print ('%15.12f' % xyz_coords[i][j], end = '   ')
		print ('\n', end = '')



def bond_distance(coords1, coords2):
    """
    Calculate distance between two cartesian coordinates
    """
    r = 0.0
    for i in range(3):
        r += (coords2[i] - coords1[i]) ** 2
    d = round(math.sqrt(r), 8)
    return d


def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]
    return c



def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))


def length(v):
  return math.sqrt(dotproduct(v, v))


def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2))) * (180.0 / (math.pi))
	

def dihedral(v1, v2, v3):
	n12 = cross(v1, v2)
	l12 = length(n12)
	e12 = [x/l12 for x in n12]
	
	n23 = cross(v2, v3)
	l23 = length(n23)
	e23 = [x/l23 for x in n23]
	
	len_v2 = length(v2)
	e_v2 = [x/len_v2 for x in v2]
	m1 = cross(e_v2, e12)
	
	x = dotproduct(e12, e23)
	y = dotproduct(m1, e23)
	
	return math.atan2(y, x) * (180.0 / math.pi)


def out_of_plane_angle(v1, v2, v3):
	e_v1 = [x/length(v1) for x in v1]
	
	n23 = cross(v2, v3)
	e23 = [x/length(n23) for x in n23]
	
	cos_angle = math.acos(dotproduct(e23, e_v1) / (length(e23) * length(e_v1))) * (180.0 / (math.pi))
	
	return 90.0 - cos_angle


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='This script analyzes a user given xyz file.')
	parser.add_argument('xyz_file', help='The filepath for the xyz file to analyze.')
	parser.add_argument('-minimum_length', help = 'The minimum distance to consider atoms bonded.', type = float, default = 0.0)
	parser.add_argument('-maximum_length', help = 'The maximum distance to consider atoms bonded.', type = float, default = 4.0)
	args = parser.parse_args()

	xyz_filename = args.xyz_file
	symbols, coord = read_xyz(xyz_filename)
	
	print_geom(len(symbols), coord)
