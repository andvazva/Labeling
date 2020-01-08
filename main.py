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

import argparse
import IO
import numpy as np
from collections import Counter
import time
import os
import shutil

def label_triangles(triangle_vertex, labels):
    triangle_label = []
    for vertices in triangle_vertex:
        label0 = labels[vertices[0]]
        label1 = labels[vertices[1]]
        label2 = labels[vertices[2]]
        if (label0 == label1 == label2):
             triangle_label.append(label2)
        elif (label0 == label1 and label0!=label2):
            triangle_label.append(label0)
        elif (label1 == label2 and (label1!= label0)):
            triangle_label.append(label1)
        elif (label0 == label2) and (label0!=label1):
            triangle_label.append(label0)
        else:
            triangle_label.append(label0)
    return triangle_label


def get_mce (labels):
    if len(labels)!=0:
        return (Counter(labels).most_common(1)[0][0])
    else:
        return -1

def get_sce(init_labels,end_labels,first_mce):
    mce_init_perc = init_labels.count(first_mce)/len(init_labels)
    mce_end_perc = end_labels.count(first_mce)/len(end_labels)
    if mce_init_perc < mce_end_perc:
        selected_list = init_labels
    else:
        selected_list = end_labels
    sce = get_mce([l for l in selected_list if l != first_mce])
    sce_percent = selected_list.count(sce)/len(selected_list)
    if sce_percent > 0.4:
        return sce
    else:
        return first_mce


def get_bundle_name(intersection,Ltriangle_label,Rtriangle_label,parcel_names):
    init_hemi = intersection[1]
    init_tri = intersection[2]
    end_hemi = intersection[4]
    end_tri = intersection[5]

    init_labels, end_labels = [], []
    hemi1 = list(set(init_hemi))
    hemi2 = list(set(end_hemi))
    if (len(hemi1) == len(hemi2) == 1):
        if hemi1[0] == hemi2[0] == "R":
            type_con = "rh"
        elif hemi1[0] == hemi2[0] == "L":
            type_con = "lh"
        else:
            type_con = "bh"
    else:
        type_con = "bh"
    for i in range(len(init_tri)):

        if init_tri[i] > 81919:
            init_tri[i] -= 81920
        if end_tri[i] > 81919:
            end_tri[i] -= 81920

        if init_hemi[i] == "L":
            init_labels.append(Ltriangle_label[init_tri[i]])
        else:
            init_labels.append(Rtriangle_label[init_tri[i]])
        if end_hemi[i] == "L":
            end_labels.append(Ltriangle_label[end_tri[i]])
        else:
            end_labels.append(Rtriangle_label[end_tri[i]])
    init_label = get_mce(init_labels)
    end_label = get_mce(end_labels)
    if init_label == end_label:
        end_label = get_sce(init_labels,end_labels,init_label)
    if init_label < end_label:
        name = type_con+"_"+parcel_names[init_label] + "-" + parcel_names[end_label]
    else:
        name = type_con+"_"+parcel_names[end_label] + "-" + parcel_names[init_label]
    return name


def delete_points(points,indices):
    new_points = []
    for i,p in enumerate(points):
        if not i in indices:
            new_points.append(p)
    return (new_points)


def align_intersection(intersection, name,Ltriangle_label,Rtriangle_label,parcel_names,filter):
    name1 = name.split("_")[1].split("-")[0]
    name2 = name.split("_")[1].split("-")[1]
    init_hemi = intersection[1]
    init_tri = intersection[2]
    end_hemi = intersection[4]
    end_tri = intersection[5]
    bad_indices = []
    for i in range(len(init_tri)):
        if init_hemi[i] == "L":
            init_name = parcel_names[Ltriangle_label[init_tri[i]]]
        else:
            init_name = parcel_names[Rtriangle_label[init_tri[i]]]
        if end_hemi[i] == "L":
            end_name = parcel_names[Ltriangle_label[end_tri[i]]]
        else:
            end_name = parcel_names[Rtriangle_label[end_tri[i]]]

        if (init_name != name1) or (end_name != name2):
            if init_name == name2 and end_name == name1:
                intersection[1][i], intersection[2][i], intersection[3][i], intersection[4][i], intersection[5][i], intersection[6][i] = \
                intersection[4][i], intersection[5][i], intersection[6][i], intersection[1][i], intersection[2][i], intersection[3][i]
            else:
                bad_indices.append(i)

    if filter == 'y':
        intersection[1] = np.delete(intersection[1], bad_indices)
        intersection[2] = np.delete(intersection[2], bad_indices)
        intersection[3] = delete_points(intersection[3], bad_indices)
        intersection[4] = np.delete(intersection[4], bad_indices)
        intersection[5] = np.delete(intersection[5], bad_indices)
        intersection[6] = delete_points(intersection[6], bad_indices)
        intersection[7] = np.delete(intersection[7], bad_indices)
        return []

    return bad_indices


def calc_center(points,bad_indices):
    npoints = len(points)-len(bad_indices)
    if npoints != 0:
        x = sum([p[0] for i,p in enumerate(points) if i not in bad_indices]) / npoints
        y = sum([p[1] for i,p in enumerate(points) if i not in bad_indices]) / npoints
        z = sum([p[2] for i,p in enumerate(points) if i not in bad_indices]) / npoints
        return (x,y,z)
    else:
        return (0,0,0)


def labeling(intersections,Ltriangle_label,Rtriangle_label,parcel_names,filter, coord, desc_sort):
    names_map = {} 
    centroids = {}
    nindices_PrCPoC,nbadindices_PrCPoC = 0,0
    for i,intersection in enumerate(intersections):
        if len(intersection[2]) != 0:
            name = get_bundle_name(intersection,Ltriangle_label,Rtriangle_label,parcel_names)
            if name not in names_map:
                names_map[name] = [i]
            else:
                names_map[name].append(i)
            bad_indices = align_intersection(intersection, name,Ltriangle_label,Rtriangle_label,parcel_names,filter)
            if name.split("-")[0]+"-"+name.split("-")[1]  == "lh_PoC-PrC":
            	nindices_PrCPoC+=len(intersection[2])
            	nbadindices_PrCPoC += len(bad_indices)
            init_points = intersection[3]
            centroids[i] = (calc_center(init_points,bad_indices))
    parcel_map = {} 
    for name,clusters in names_map.items():
        name_init = name.split("-")[0]

        index = 0 
        y_values = [centroids[i][coord]  for i in clusters]
        if desc_sort == 'y':
            sorted_clusters = [x for _,x in sorted(zip(y_values,clusters),reverse = True)]
        else:
            sorted_clusters = [x for _, x in sorted(zip(y_values, clusters))]
        for cluster_index in sorted_clusters:
            cluster_name = intersections[cluster_index][0]
            if filter == 'y':
                parcel_map[name + "_" + str(index)] = intersections[cluster_index]
            else:
                parcel_map[name+"_"+str(index)] = cluster_name
            index+=1
    print("Número de fibras total PrC_PoC: "+str(nindices_PrCPoC))
    print("Número de fibras que se salen PrC_PoC: "+str(nbadindices_PrCPoC))
    print("Porcentaje de fibras que se salen PrC_PoC: "+str(round(((nbadindices_PrCPoC/nindices_PrCPoC)*100),2))+"%")
    print("Número de clusters etiquetados en la parcela PoC: "+str(len(parcel_map)))
    return parcel_map

def main():
    parser = argparse.ArgumentParser(description='Perform clustering on a dataset of streamlines')
    parser.add_argument('--intersection-dir',type= str,help='Input intersection file of both hemispheres')
    parser.add_argument('--clusters-dir', type=str, help="Input directory of clusters .bundles")
    parser.add_argument('--output-path', type = str, help="Output path")
    parser.add_argument('--subject-name', type = str, help="Name of input subject")
    parser.add_argument('--coord', type=str, default = 'y', help="Type coordinate to sort bundles: 'x' 'y' or 'z'")
    parser.add_argument('--desc-sort', type=str, default = 'n', help="Type 'y' to descendant sort")
    parser.add_argument('--filter', type=str, default = 'y',help="'y' to filter fibers 'n' to not filter fibers")

    args = parser.parse_args()
    start = time.time()
    if args.coord == 'x':
        selected_coord = 0
    elif args.coord == 'y':
        selected_coord = 1
    elif args.coord == 'z':
        selected_coord = 2

    LMesh = "input_data/lh.r.aims.pial.obj"
    RMesh = "input_data/rh.r.aims.pial.obj"
    Llavels = "input_data/lh_labels.txt"
    Rlavels = "input_data/rh_labels.txt"
    desikan_atlas = "input_data/desikan_atlas.txt"

    if not os.path.exists(args.output_path):
    	os.mkdir(args.output_path)
    out_path = args.output_path+"/"+args.subject_name
    if os.path.exists(out_path):
    	shutil.rmtree(out_path)
    IO.create_path(out_path)
    Loutput_filename = out_path+"/"+args.subject_name+"_left-hemi"
    Routput_filename = out_path+"/"+args.subject_name+"_right-hemi"
    Boutput_filename = out_path+"/"+args.subject_name+"_both-hemi"
    IO.create_path(Loutput_filename)
    IO.create_path(Routput_filename)
    IO.create_path(Boutput_filename)

    intersections = IO.read_intersections(args.intersection_dir)

    Lvertex_coord, Ltriangles_vertex = IO.read_mesh_obj(LMesh) 
    Lvertex_labels = IO.read_vertex_labels(Llavels)  
    Rvertex_coord, Rtriangles_vertex = IO.read_mesh_obj(RMesh)  
    Rvertex_labels = IO.read_vertex_labels(Rlavels)

    parcel_names, parcel_acr = IO.read_parcels(desikan_atlas)   
    Ltriangle_label = label_triangles(Ltriangles_vertex,Lvertex_labels) 
    Rtriangle_label = label_triangles(Rtriangles_vertex,Rvertex_labels)  

    parcel_map = labeling(intersections,Ltriangle_label,Rtriangle_label,parcel_acr,args.filter,selected_coord,args.desc_sort)
    if args.filter == 'y':
        IO.write_clusters(parcel_map, args.clusters_dir, args.output_path, args.subject_name)
    else:
        IO.rename_clusters(parcel_map,args.clusters_dir,args.output_path,args.subject_name)
    end = time.time()
    print("Tiempo de ejecución: "+str(end-start))
if __name__ == '__main__':
    main()