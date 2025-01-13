#!/usr/bin/env python 
# Compares two slabs ans checks their cmnmensurability 
################################################################################
# Copyright Jose J. Plata, J. Amaya Suarez, E. R: Remesal, Antonio M. Marquez  #
# and Javier Fdez Sanz		   (2019)  								           #
#                                                                              #
# This is free software: you can                                               #
# redistribute it and/or modify it under the terms of the GNU General Public   #
# License as published by the Free Software Foundation, either version 3 of    #
# the License, or (at your option) any later version.                          #
# This program is distributed in the hope that it will be useful, but WITHOUT  #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for    #
# more details.                                                                #
# You should have received a copy of the GNU General Public License along with #
# this program.  If not, see <http://www.gnu.org/licenses/>.                   #
#                                                                              #
################################################################################

import os
import sys
import re
import numpy as np
import xml
import xml.etree.ElementTree as ET
import math
import argparse
from operator import itemgetter
from scipy.integrate import simpson
from scipy.integrate import trapezoid
from scipy.integrate import cumulative_trapezoid

listR = []

def matGen(N):

	genMatkk = []
	for ii in range(1,N+1):
		if N % ii == 0:
			jj = N/ii
			for j in range(int(jj)):
				kk = np.zeros((2,2))
				kk[0][0]=ii
				kk[0][1]=j
				kk[1][0]=0
				kk[1][1]=jj
				genMatkk.append(kk)	
		
	return genMatkk

def dotproduct(v1, v2):
	return sum((a*b) for a, b in zip(v1, v2))

def length(v):
	return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
	return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))


def getSuperSlabs(listMat,surfPrim):

	listOfentries = []
	for i in range(len(listMat)):
		entry = [0.0]*3
		#Generate new SuperSlab
		newLat =  np.dot(listMat[i],surfPrim)
		#np.matmul(listMat[i], surfPrim)
		#Compute modules and angles
		u =  np.linalg.norm(newLat[0])
		v =  np.linalg.norm(newLat[1])
		alfa = angle(newLat[0],newLat[1]) 		
		if u < v:
			entry[0] = u	
			entry[1] = v
		else:
			entry[0] = v	
			entry[1] = u
		entry[2] = alfa
		listOfentries.append(entry)

	return listOfentries	


def compareSuperLattice(nSub,nTF,surfSub,surfTF,errAMax,arrAlphaMax):
	
	genMatSub = [] 
	genMatTF=[]
	#Create superslabs generator matrixes
	genMatSub = matGen(nSub)
	genMatTF =  matGen(nTF)

	#Create superslab for Sub
	superSlabSub = getSuperSlabs(genMatSub,surfSub)
	#Create superslab for TF
	superSlabTF = getSuperSlabs(genMatTF,surfTF)

	areaPrimSub = np.linalg.norm(np.cross(surfSub[0],surfSub[1]))
	tempMCIA = nSub*areaPrimSub
	for i in range(len(superSlabSub)):
		for j in range(len(superSlabTF)):
			entry = []
			errA = 100*(superSlabSub[i][0]-superSlabTF[j][0])/superSlabTF[j][0]
			errB = 100*(superSlabSub[i][1]-superSlabTF[j][1])/superSlabTF[j][1]
			errAlpha = 100*(superSlabSub[i][2]-superSlabTF[j][2])/superSlabTF[j][2]
			if errA < errAMax and errB < errAMax and errAlpha < arrAlphaMax:
				entry.append(tempMCIA)
				entry.append(nSub)
				entry.append(nTF)
				entry.append(errA)
				entry.append(errB)
				entry.append(errAlpha)
				entry.append(genMatSub[i])
				entry.append(genMatTF[j])
				listR.append(entry)		

def readSurf(infile):

	if not os.path.isfile(infile):
		print("File path {} does not exist. Exiting...".format(infile))
		sys.exit()
	f = open(infile, "r")
	lines = f.readlines()
	f.close()
	#Read factor
	chain=lines[1].split()
	factor=float(chain[0])
	#Read Lattice
	latt = np.zeros((3,3))
	for i in range(2,5):
		chain=lines[i].split()
		x = np.array(chain)
		latt[i-2] = x.astype(float)
	#Facro*Lattice
	latt = factor * latt
    #SurfaceLatt
	latt2D = np.zeros((2,2))
	for i in range(2):
		for j in range(2):
			latt2D[i][j] = latt[i][j]
	return latt2D
	
		
def calcAreasList(areaPrimTF,areaPrimSub,maxMCIA,maxAreaErr):

	listAreas = []
	MCIA = 0.0
	it = 1
	while MCIA < maxMCIA:
		MCIA =  areaPrimSub*float(it)
		it2 = 1
		dummyArea = 0.0
		while dummyArea < maxMCIA:
			dummyArea = areaPrimTF*float(it2)
			difArea = 100.0*(dummyArea-MCIA)/dummyArea
			if abs(difArea) < maxAreaErr:
				dummyVec = [0]*2
				dummyVec[0] = it
				dummyVec[1]=it2
				listAreas.append(dummyVec)
			it2 = it2 +1
		it = it +1	
	return listAreas;


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("Substrate", help="File containing substarte POSCAR.")
	parser.add_argument("ThinFilm",  help="File containing thin-film POSCAR.")
	parser.add_argument("-m", "--maxMCIA", type=float, default=100.0, help="Maximum Minimum coincident interface area to be considered in Angstrom.")
	parser.add_argument("-a", "--maxAreaErr",    type=float, default=12.0, help="Maximum percentual error between surpercell surface areas.")
	parser.add_argument("-u", "--maxError_u",    type=float, default=5.0, help="Maximum percentual error between lattices.")
	parser.add_argument("-w", "--maxError_alpha",type=float, default=5.0, help="Maximum percentual error between angles.")
	args = parser.parse_args()

	#Read unit surface cell for thin film POSCAR 
	surfTF=readSurf(args.ThinFilm)
	#Calculate  unit surface cell  area for thin film POSCAR 
	areaPrimTF = np.linalg.norm(np.cross(surfTF[0],surfTF[1]))
	#Read unit surface cell for substrate film POSCAR 
	surfSub=readSurf(args.Substrate)
	#Calculate  unit surface cell  area for substrate POSCAR 
	areaPrimSub = np.linalg.norm(np.cross(surfSub[0],surfSub[1]))

	#Check tentative MCIA and buld N1/N2 lookup table
	listAreas = calcAreasList(areaPrimTF,areaPrimSub,args.maxMCIA,args.maxAreaErr);
	for i in range(len(listAreas)):
		compareSuperLattice(listAreas[i][0],listAreas[i][1],surfSub,surfTF,args.maxError_u,args.maxError_alpha)

	#Sort by MCIA
	sorted(listR, key=itemgetter(0))

	#Print results
	filename = "results.dat"
	f= open(filename,"w+")
	f.write("# Conmensurability between {:>18} as substrate and {:>18} as thin film using epitaxyNML \n".format(args.Substrate,args.ThinFilm))
	f.write("#####################################################################################################################################\n")
	f.write("#{:>10}{:>15}{:>15}{:>14}{:>14}{:>15}{:>15}{:>25}        #\n".format("MCIA","N_TF","N_Sub","Err_a","Err_b","Err_alpha","MAT_TF","MAT_Sub"))
	f.write("#####################################################################################################################################\n")

	for entry in listR:
		f.write('{:>14.6f}{:>14.6f}{:>14.6f}{:>14.6f}{:>14.6f}{:>14.6f}  [{:4d},{:4d},{:4d},{:4d}]  [{:4d},{:4d},{:4d},{:4d}]\n'.format(entry[0], entry[2], entry[1], entry[3], entry[4], entry[5], int(entry[7][0][0]), int(entry[7][0][1]), int(entry[7][1][0]), int(entry[7][1][1]), int(entry[6][0][0]), int(entry[6][0][1]), int(entry[6][1][0]), int(entry[6][1][1])))
	f.write("#####################################################################################################################################\n")
	f.close()



if __name__ == '__main__':
	main()

