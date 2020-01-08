# Copyright (C) 2019  Andrea Vázquez Varela

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


#Authors:
# Narciso López López
# Andrea Vázquez Varela
#Creation date: 09/12/2019
#Last update: 11/12/2019

import sys
import glob
import numpy as np 
import bundleTools as bT
import bundleMetrics as bM
import shutil
import os
from collections import OrderedDict 
import operator
from scipy.optimize import linear_sum_assignment

thr = 30

def isSquare (m): return all (len (row) == len (m) for row in m)


def file_per_bundle(path,hemi):
	fmap = {}
	snames = os.listdir(path)
	for sname in snames:
		bundle_files = glob.glob(path+"/"+sname+"/"+sname+"_"+hemi+"/*.bundles")
		for f in bundle_files:
			bname = os.path.basename(f).split("_")[1]
			if bname not in fmap:
				fmap[bname] = [f]
			else:
				fmap[bname].append(f)
	return fmap

def read_bundles(files):
	smap = OrderedDict()
	for f in files:
		sname = f.split("/")[-2].split("_")[0]
		bundle = bT.read_bundle(f)
		if len(bundle)>0:
			if sname not in smap:
				smap[sname] = [bundle]
			else:
				smap[sname].append(bundle)
	smap = OrderedDict(sorted(smap.items(), key = lambda x : len(x[1]), reverse = True))
	return smap

def write_bundles(labelmap,conection,hemi,outpath):
	if hemi == "left-hemi":
		h = "lh_"
	elif hemi == "rigth-hemi":
		h = "rh_"
	else:
		h = "bh_"
	if not os.path.exists(outpath+"/"+hemi):
		os.mkdir(outpath+"/"+hemi)
	for key,bundles in labelmap.items():
		splitted = key.split("_")
		sname = splitted[0]
		bindex = splitted[1]
		if not os.path.exists(outpath+"/"+hemi+"/"+sname):
			os.mkdir(outpath+"/"+hemi+"/"+sname)

		bT.write_bundle(outpath+"/"+hemi+"/"+sname+"/"+h+conection+"_"+bindex+".bundles",bundles)

def matching(smap):
	labeling_map = {}
	index = 0
	while smap:
		s1,bundles1 = smap.popitem(0)
		reference_labels = []
		for i in range(len(bundles1)):
			labeling_map[s1+"_"+str(index)] = bundles1[i]
			reference_labels.append(index)
			index+=1

		for s2,bundles2 in smap.items():
			centroids1 = [bM.calc_centroid(b) for b in bundles1]
			centroids2 = [bM.calc_centroid(b) for b in bundles2]
			matrix_dist = bM.matrix_dist_group(centroids1, centroids2)

			row_ind, col_ind = linear_sum_assignment(matrix_dist)

			delete_indices = []
			for i in range(len(row_ind)):
				if matrix_dist[row_ind[i]][col_ind[i]] < thr:
					labeling_map[s2+"_"+str(reference_labels[row_ind[i]])] = bundles2[col_ind[i]]
					delete_indices.append(col_ind[i])
			smap[s2] = [b for i,b in enumerate(bundles2) if i not in delete_indices]

	return labeling_map
			


			
def main():
	inpath = sys.argv[1]
	outpath = sys.argv[2]

	if os.path.exists(outpath):
		shutil.rmtree(outpath)
	os.mkdir(outpath)

	lfiles_map = file_per_bundle(inpath,"left-hemi")
	rfiles_map = file_per_bundle(inpath,"rigth-hemi")
	bfiles_map = file_per_bundle(inpath,"both-hemi")


	files_list = [lfiles_map,rfiles_map,bfiles_map]
	hemis = ["left-hemi","rigth-hemi","both-hemi"]
	for i, files_map in enumerate(files_list):
		hemi = hemis[i]
		for conection,files in lfiles_map.items():
			smap = read_bundles(files)
			labeling_map = matching(smap)
			write_bundles(labeling_map,conection,hemi,outpath)




if __name__=="__main__":
	main()