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
# Andrea Vázquez Varela
# Narciso López López
#Creation date: 15/10/2019
#Last update: 30/10/2019

from dipy.segment.metric import AveragePointwiseEuclideanMetric
from dipy.segment.clustering import QuickBundles #PARA LA MDF
import dipy
import bundleTools as bT
import numpy as np
import argparse
import time
import os
import shutil
import dipy.segment.metric as dipymetric
import itertools
import numpy as np
from dipy.segment.metric import ResampleFeature
from dipy.tracking.streamline import set_number_of_points
import shutil
from dipy.segment.metric import Metric

class MaxDirectFlipMetric(Metric):
    def __init__(self, feature):
        super(MaxDirectFlipMetric, self).__init__(
             feature=feature)

    @property
    def is_order_invariant(self):
        return True  # Ordering is handled in the distance computation

    def are_compatible(self, shape1, shape2):
            return shape1[0] == shape2[0]

    def dist(self, v1, v2):
        def max_euclidean(x, y):
            return np.max(norm(x-y, axis=1))
        dist_direct = max_euclidean(v1, v2)
        dist_flipped = max_euclidean(v1, v2[::-1])
        return min(dist_direct, dist_flipped)


def norm(x, ord=None, axis=None):
    if axis is not None:
        return np.apply_along_axis(np.linalg.norm, axis,
                x.astype(np.float64), ord)


def save_clusters_fibers(clusters, filename):
    text=""
    file = open(filename,"w+")
    for i,cluster in enumerate(clusters):
        text = text + str(i)+ "-"
        for fiber_index in cluster.indices:
            text += str(fiber_index)+" "
        file.write(text+"\n")
        text =""
    file.close()

# Parse input arguments
parser = argparse.ArgumentParser(description='QuickBundles')
parser.add_argument('--infile', help='Input streamlines file')
parser.add_argument('--outfile', type = str, help='output directory')
parser.add_argument('--thr', type = float, help='Distance threshold')

args = parser.parse_args()


streamlines = np.asarray(bT.read_bundle(args.infile))


feature = ResampleFeature(nb_points=12)
max_metric = MaxDirectFlipMetric(feature=feature)

#max_metric = MaxDistnace(feature)
qb = QuickBundles(threshold=args.thr) #PARA MDF
#qb = QuickBundles(threshold=args.thr, metric = max_metric)
start_time = time.time()
# Clusters is a ClusterMap object which contains attributes that provide information about the clustering result
clusters = qb.cluster(streamlines)
# Execution time
print("Time spent: --- %s seconds ---" % (time.time() - start_time))



if os.path.exists(args.outfile):
    shutil.rmtree(args.outfile)
os.mkdir(args.outfile)
save_clusters_fibers(clusters, args.outfile+str(args.thr)+".txt")
centroids = []
bundles_dir = args.outfile+"/finalBundles"
os.mkdir(bundles_dir)
for i,cluster in enumerate(clusters):
	if len(cluster) > 1:
		new_cluster = set_number_of_points(cluster,21)
		centroid = set_number_of_points(cluster.centroid,21)
		final_bundles_file = bundles_dir+"/"+str(i)+".bundles"
		#c = np.asarray(cluster) #[:]
		centroids.append(np.asarray(centroid))
		bT.write_bundle(final_bundles_file,np.asarray(new_cluster))
bT.write_bundle(args.outfile+"/centroids.bundles",np.asarray(centroids))

