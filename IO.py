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
#Creation date: 27/08/2019
#Last update: 11/12/2019

import os
import numpy as np
import shutil
import bundleTools as bT


def create_path(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


def read_intersection_file(infile):
    with open(infile) as f:
        content = f.readlines()
        total_triangles = int(content[0])
        InTri = [int(tri) for tri in (content[1].split(" ")[:-1])]
        FnTri = [int(tri) for tri in (content[2].split(" ")[:-1])]

        InPoints = []
        init_points = [float(tri) for tri in (content[3].split(" ")[:-1])]
        for i in range(0,len(init_points), 3):
            coords = np.zeros(3)
            coords[0] = init_points[i]
            coords[1] = init_points[i + 1]
            coords[2] = init_points[i + 2]
            InPoints.append(coords)

        FnPoints = []
        fin_points = [float(tri) for tri in (content[4].split(" ")[:-1])]
        for i in range(0,len(fin_points), 3):
            coords = np.zeros(3)
            coords[0] = fin_points[i]
            coords[1] = fin_points[i + 1]
            coords[2] = fin_points[i + 2]
            FnPoints.append(coords)

        fib_index = [int(tri) for tri in (content[5].split(" ")[:-1])]

    f.close();
    return np.asarray(InTri), np.asarray(FnTri), np.asarray(InPoints), np.asarray(FnPoints), np.asarray(fib_index);


def read_intersections(path):
    files = os.listdir(path)
    files = sorted(files)
    intersections = []
    for name in files:
        inter_path = path+"/"+name
        intersections.append(read_intersection(inter_path,name))
    return intersections

def read_intersection(infile,name):
    with open(infile) as f:
        content = f.readlines()
        total_triangles = int(content[0])
        InHemi = [hem for hem in (content[1].split(" ")[:-1])]
        FnHemi = [hem for hem in (content[2].split(" ")[:-1])]
        InTri = [int(tri) for tri in (content[3].split(" ")[:-1])]
        FnTri = [int(tri) for tri in (content[4].split(" ")[:-1])]

        InPoints = []
        init_points = [float(tri) for tri in (content[5].split(" ")[:-1])]
        for i in range(0, len(init_points), 3):
            coords = np.zeros(3)
            coords[0] = init_points[i]
            coords[1] = init_points[i + 1]
            coords[2] = init_points[i + 2]
            InPoints.append(coords)

        FnPoints = []
        fin_points = [float(tri) for tri in (content[6].split(" ")[:-1])]
        for i in range(0, len(fin_points), 3):
            coords = np.zeros(3)
            coords[0] = fin_points[i]
            coords[1] = fin_points[i + 1]
            coords[2] = fin_points[i + 2]
            FnPoints.append(coords)

        fib_index = [int(tri) for tri in (content[7].split(" ")[:-1])]

    f.close();
    return [name,np.asarray(InHemi), np.asarray(InTri), np.asarray(InPoints), np.asarray(FnHemi), np.asarray(FnTri), np.asarray(
        FnPoints), np.asarray(fib_index)];


def read_mesh_obj(mesh_path):
    vertex = []
    triangles = []
    with open(mesh_path) as file:
        lines = file.readlines()
        for line in lines:
            splitted = line.split(" ")
            id = str(splitted[0])
            if (id == 'v'):
                coords = np.asarray([float(splitted[1]), float(splitted[2]), float(splitted[3])])
                vertex.append(coords)
            if (id == 'f'):
                vs = np.asarray([int(splitted[1].split("//")[0]) - 1, int(splitted[2].split("//")[0]) - 1,
                                int(splitted[3].split("//")[0]) - 1])
                triangles.append(vs)
    file.close()
    return np.asarray(vertex), np.asarray(triangles);


def read_vertex_labels(labels_path):
    labels = []
    with open(labels_path, 'r') as file:
        for line in file:
            label = int(line)
            if label == -1:
                label = 0
            labels.append(label)
    return(np.asarray(labels))


def read_parcels(names_path):
    names = []
    bundles = []
    with open(names_path, 'r') as file:
        for line in file:
            spplited = line.split(" ")
            names.append(str(spplited[0].rstrip("\n\r")))
            bundles.append(str(spplited[1].rstrip("\n\r")))
    return(np.asarray(names), np.asarray(bundles))


def rename_clusters(parcel_map,clusters_dir,output,sname):
    for new_name,intersection in parcel_map.items():
        if new_name.split("_")[0] == "lh":
            path = output+"/"+sname+"/"+sname+"_left-hemi"
        elif new_name.split("_")[0] == "rh":
            path = output+"/"+sname+"/"+sname+"_right-hemi"
        else:
            path = output+"/"+sname+"/"+sname+"_both-hemi"
        for new_name,name in parcel_map.items():
            old_name = name.split(".")[0]
            bundles_old = clusters_dir+"/"+old_name+".bundles"
            shutil.copyfile(bundles_old,path+"/"+new_name+".bundles")
            shutil.copyfile(bundles_old+"data",path+"/"+new_name + ".bundlesdata")

def write_clusters(parcel_map,clusters_dir,output,sname):
    for new_name,intersection in parcel_map.items():
        if new_name.split("_")[0] == "lh":
            path = output+"/"+sname+"/"+sname+"_left-hemi"
        elif new_name.split("_")[0] == "rh":
            path = output+"/"+sname+"/"+sname+"_right-hemi"
        else:
            path = output+"/"+sname+"/"+sname+"_both-hemi"
        name = intersection[0]
        old_name = name.split(".")[0]
        bundles_old = clusters_dir+"/"+old_name+".bundles"
        old_fibers = bT.read_bundle(bundles_old)
        new_fibers = []
        for findex in intersection[7]:
            new_fibers.append(old_fibers[findex])
        bT.write_bundle(path+"/"+new_name+".bundles",new_fibers)
