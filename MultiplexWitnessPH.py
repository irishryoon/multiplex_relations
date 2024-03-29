#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 23:16:41 2023

@author: stolz
"""


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial
from scipy.spatial import Delaunay
#from scipy.spatial import Voronoi, voronoi_plot_2d

import gudhi
import gudhi.wasserstein
import gudhi.representations

from sklearn.cluster import KMeans
#import ot


# Functions

def getDelaunayTriangulation(zero_simplex_coords):

    #define delauney triangulation on vertices

    tri = Delaunay(zero_simplex_coords)

    #plt.triplot(zero_simplex_coords[:,0], zero_simplex_coords[:,1], tri.simplices)
    #plt.show()


    #we now extract all simplices from the triangulation

    two_simplices_coords = zero_simplex_coords[tri.simplices]

    one_simplices_indices = np.empty((0,2), int)
    two_simplices_indices = np.empty((0,3),int)

    for triangle in two_simplices_coords:  
    
        point_coords = [300,300,300]
    
        for point_number in range(0,3):
        
            [point_coords[point_number]], = np.where(np.all(zero_simplex_coords==triangle[point_number,:],axis = 1))
        
        one_simplices_indices = np.append(one_simplices_indices,[np.sort([point_coords[0],point_coords[1]])],axis = 0)
        one_simplices_indices = np.append(one_simplices_indices,[np.sort([point_coords[1],point_coords[2]])],axis = 0)
        one_simplices_indices = np.append(one_simplices_indices,[np.sort([point_coords[0],point_coords[2]])],axis = 0)
                                                                                                                      
        two_simplices_indices = np.append(two_simplices_indices,[np.sort([point_coords[0],point_coords[1],point_coords[2]])],axis = 0)

    one_simplices_indices = np.unique(one_simplices_indices,axis = 0)
    
    return one_simplices_indices, two_simplices_indices


def getFiltrationValues(witness_coords, zero_simplex_coords, one_simplices_indices, two_simplices_indices):

    dist_vessels = scipy.spatial.distance.cdist(witness_coords,zero_simplex_coords)


    ## We look at filtration values 

    min_values = dist_vessels.min(axis=1)

    landmark_indices = np.empty((0,1), int)
    edge_endpoint_indices = np.empty((0,1), int)
    third_triangle_point_indices = np.empty((0,1), int)
 
    for cell_index in range(np.shape(witness_coords)[0]):
    
        new_landmark, = np.where(dist_vessels[cell_index,] == min_values[cell_index])
        one_removed_minimum = np.delete(dist_vessels[cell_index,], new_landmark)
        #dist_vessels_one_removed = np.append(dist_vessels_one_removed,[one_removed_minimum],axis = 0)

        second_min_value = min(one_removed_minimum)
        new_edge_endpoint, = np.where(dist_vessels[cell_index,] == second_min_value)

        new_edge_endpoint_new_numbering, = np.where(one_removed_minimum == second_min_value)
        two_removed_minima = np.delete(one_removed_minimum, new_edge_endpoint_new_numbering)
    
        third_min_value = min(two_removed_minima)
        new_triangle_endpoint, = np.where(dist_vessels[cell_index,] == third_min_value)
        
        if len(new_landmark) > 1:
            
            if len(new_edge_endpoint) > 1 or len(new_triangle_endpoint)>1:
                
                print('Help we have multiple witnesses for simplces of different dimensions')
            
            for landmark in new_landmark:
                
                landmark_indices = np.append(landmark_indices,landmark)
                edge_endpoint_indices = np.append(edge_endpoint_indices,new_edge_endpoint)
                third_triangle_point_indices = np.append(third_triangle_point_indices,new_triangle_endpoint)
                    
        if len(new_edge_endpoint) > 1:
            
            if len(new_triangle_endpoint)>1:
                
                print('Help we have multiple witnesses for simplces of different dimensions')
            
            for endpoint in new_edge_endpoint:
                
                landmark_indices = np.append(landmark_indices,new_landmark)
                edge_endpoint_indices = np.append(edge_endpoint_indices,endpoint)
                third_triangle_point_indices = np.append(third_triangle_point_indices,new_triangle_endpoint)
        
        if len(new_triangle_endpoint) > 1:
            
            for triangle_endpoint in new_triangle_endpoint:
                
                landmark_indices = np.append(landmark_indices,new_landmark)
                edge_endpoint_indices = np.append(edge_endpoint_indices,new_edge_endpoint)
                third_triangle_point_indices = np.append(third_triangle_point_indices,triangle_endpoint)
        
        if (len(new_landmark) == 1 and len(new_edge_endpoint) == 1 and len(new_triangle_endpoint) == 1):
            
            landmark_indices = np.append(landmark_indices,new_landmark)
            edge_endpoint_indices = np.append(edge_endpoint_indices,new_edge_endpoint)
            third_triangle_point_indices = np.append(third_triangle_point_indices,new_triangle_endpoint)





    witness_triangles = np.column_stack((landmark_indices,edge_endpoint_indices,third_triangle_point_indices))
    witness_triangles = np.sort(witness_triangles,axis = 1)
    filtration_values_two_simplices = np.empty((0,1),int)

    witness_edges = np.column_stack((landmark_indices,edge_endpoint_indices))
    witness_edges = np.sort(witness_edges,axis = 1)


    filtration_values_one_simplices = np.zeros((1,one_simplices_indices.shape[0]),int)
    filtration_values_vertices_prelim = np.zeros((1,zero_simplex_coords.shape[0]),int)


    # We look at filtration values for 2- simplices

    for two_simplex in two_simplices_indices:
    
        witnesses_for_triangle, = np.where(np.all(witness_triangles==two_simplex,axis=1))
        filtration_values_two_simplices = np.append(filtration_values_two_simplices,len(witnesses_for_triangle))
        
        if len(witnesses_for_triangle) != 0:
    
            edge1 = [two_simplex[0],two_simplex[1]]
            edge2 = [two_simplex[1],two_simplex[2]]
            edge3 = [two_simplex[0],two_simplex[2]]
            edge1_index, = np.where(np.all(one_simplices_indices==edge1,axis=1))
            edge2_index, = np.where(np.all(one_simplices_indices==edge2,axis=1))
            edge3_index, = np.where(np.all(one_simplices_indices==edge3,axis=1))
    
            filtration_values_one_simplices[0,edge1_index] = max([len(witnesses_for_triangle)],filtration_values_one_simplices[0,edge1_index])
            filtration_values_one_simplices[0,edge2_index] = max([len(witnesses_for_triangle)],filtration_values_one_simplices[0,edge2_index])
            filtration_values_one_simplices[0,edge3_index] = max([len(witnesses_for_triangle)],filtration_values_one_simplices[0,edge3_index])


    # We look at filtration values for 1- simplices

    i = 0

    for one_simplex in one_simplices_indices:
        
        witnesses_for_edge, = np.where(np.all(witness_edges==one_simplex,axis=1))
    
        if len(witnesses_for_edge) >  filtration_values_one_simplices[0,i]:
            
            filtration_values_one_simplices[0,i] = len(witnesses_for_edge)
      
        filtration_values_vertices_prelim[0,one_simplex[0]] = max(filtration_values_one_simplices[0,i],filtration_values_vertices_prelim[0,one_simplex[0]])
        filtration_values_vertices_prelim[0,one_simplex[1]] = max(filtration_values_one_simplices[0,i],filtration_values_vertices_prelim[0,one_simplex[1]])
     
        i = i + 1

    #We count the number of witnesses per landmark

    filtration_values_vertices_witness,bins = np.histogram(landmark_indices, bins= np.arange(0,zero_simplex_coords.shape[0]+1,1))

    filtration_values_vertices = np.amax(np.column_stack((filtration_values_vertices_witness,filtration_values_vertices_prelim[0,:])),axis = 1)

    return filtration_values_vertices, filtration_values_one_simplices, filtration_values_two_simplices



def get_rel_barcode(filtration_values_vertices,filtration_values_one_simplices,filtration_values_two_simplices,one_simplices_indices,two_simplices_indices):

    st = gudhi.SimplexTree()
    maximum_number_of_witnesses = max([max(filtration_values_vertices),max(filtration_values_one_simplices[0,:]),max(filtration_values_two_simplices)])


    #filtration_values = np.empty((0,1), int)
    
    if maximum_number_of_witnesses != 0:

        for vertex in range(filtration_values_vertices.shape[0]):

            st.insert([vertex], filtration=(maximum_number_of_witnesses-filtration_values_vertices[vertex])/maximum_number_of_witnesses)
            #print([vertex], maximum_number_of_witnesses-filtration_values_vertices[vertex])
            #filtration_values = np.append(filtration_values,maximum_number_of_witnesses-filtration_values_vertices[vertex])
 
        i = 0
        for one_simplex in one_simplices_indices:

            st.insert(one_simplex, filtration=(maximum_number_of_witnesses-filtration_values_one_simplices[0,i])/maximum_number_of_witnesses)
            #print(one_simplex, maximum_number_of_witnesses-filtration_values_one_simplices[0,i])
            #filtration_values = np.append(filtration_values,maximum_number_of_witnesses-filtration_values_one_simplices[0,i])
            i = i + 1

        i = 0
        for two_simplex in two_simplices_indices:

            st.insert(two_simplex, filtration=(maximum_number_of_witnesses-filtration_values_two_simplices[i])/maximum_number_of_witnesses)
            #print(two_simplex, maximum_number_of_witnesses-filtration_values_two_simplices[i])
            #filtration_values = np.append(filtration_values,maximum_number_of_witnesses-filtration_values_two_simplices[i])
            i = i + 1
    else:
        
        for vertex in range(filtration_values_vertices.shape[0]):

            st.insert([vertex], filtration=1)
            #print([vertex], maximum_number_of_witnesses-filtration_values_vertices[vertex])
            #filtration_values = np.append(filtration_values,maximum_number_of_witnesses-filtration_values_vertices[vertex])
 
        i = 0
        for one_simplex in one_simplices_indices:

            st.insert(one_simplex, filtration=1)
            #print(one_simplex, maximum_number_of_witnesses-filtration_values_one_simplices[0,i])
            #filtration_values = np.append(filtration_values,maximum_number_of_witnesses-filtration_values_one_simplices[0,i])
            i = i + 1

        i = 0
        for two_simplex in two_simplices_indices:

            st.insert(two_simplex, filtration=1)
            #print(two_simplex, maximum_number_of_witnesses-filtration_values_two_simplices[i])
            #filtration_values = np.append(filtration_values,maximum_number_of_witnesses-filtration_values_two_simplices[i])
            i = i + 1

    barcode = st.persistence()
     
    return st,barcode


def get_rel_distance_vector(data):
    #define vertices of simplicial complex

    vertices = np.empty((0,2), int)

    vertices_x = data.loc[data["celltypes"] == "Vessel", "points_x"]
    vertices_y = data.loc[data["celltypes"] == "Vessel", "points_y"]

    for vertex in range(len(vertices_x)):
        
        vertices = np.append(vertices,[[vertices_x.iloc[vertex],vertices_y.iloc[vertex]]],axis = 0) 

    # get Delanay triangulation

    one_simplices_indices, two_simplices_indices = getDelaunayTriangulation(vertices)


    #get coordinates of other cell types and distances to blood vessels, assign filtration values, run filtration


    ### 2. Tumour

    Tumour_coords = np.empty((0,2), int)

    Tumour_coords_x = data.loc[data["celltypes"] == "Tumour", "points_x"]
    Tumour_coords_y = data.loc[data["celltypes"] == "Tumour", "points_y"]

    for vertex in range(len(Tumour_coords_x)):
        
        Tumour_coords = np.append(Tumour_coords,[[Tumour_coords_x.iloc[vertex],Tumour_coords_y.iloc[vertex]]],axis = 0) 


    #Get filtration values

    Tumour_filtration_values_vertices, Tumour_filtration_values_one_simplices, Tumour_filtration_values_two_simplices = getFiltrationValues(Tumour_coords, vertices, one_simplices_indices, two_simplices_indices)

    # We feed the filtration values into gudhi and get barcode

    st_Tumour, barcode_Tumour = get_rel_barcode(Tumour_filtration_values_vertices,Tumour_filtration_values_one_simplices,Tumour_filtration_values_two_simplices,one_simplices_indices,two_simplices_indices)
    #gudhi.plot_persistence_diagram(barcode_Tumour)
    #gudhi.plot_persistence_barcode(barcode_Tumour)
    #plt.show()



    #### 3. Necrotic

    Necrotic_coords = np.empty((0,2), int)

    Necrotic_coords_x = data.loc[data["celltypes"] == "Necrotic", "points_x"]
    Necrotic_coords_y = data.loc[data["celltypes"] == "Necrotic", "points_y"]

    for vertex in range(len(Necrotic_coords_x)):
        
        Necrotic_coords = np.append(Necrotic_coords,[[Necrotic_coords_x.iloc[vertex],Necrotic_coords_y.iloc[vertex]]],axis = 0) 

    #Get filtration values

    Necrotic_filtration_values_vertices, Necrotic_filtration_values_one_simplices, Necrotic_filtration_values_two_simplices = getFiltrationValues(Necrotic_coords, vertices, one_simplices_indices, two_simplices_indices)

    # We feed the filtration values into gudhi and get barcode

    st_Necrotic, barcode_Necrotic = get_rel_barcode(Necrotic_filtration_values_vertices,Necrotic_filtration_values_one_simplices,Necrotic_filtration_values_two_simplices,one_simplices_indices,two_simplices_indices)
    #gudhi.plot_persistence_diagram(barcode_Necrotic)
    #gudhi.plot_persistence_barcode(barcode_Necrotic)
    #plt.show()


    #### 4. Macrophage M1

    Macrophage_coords = np.empty((0,2), int)

    Macrophage_coords_x = data.loc[(data["celltypes"] == "Macrophage"), "points_x"]
    Macrophage_coords_y = data.loc[(data["celltypes"] == "Macrophage"), "points_y"]

    for vertex in range(len(Macrophage_coords_x)):
         
        Macrophage_coords = np.append(Macrophage_coords,[[Macrophage_coords_x.iloc[vertex],Macrophage_coords_y.iloc[vertex]]],axis = 0) 


    #Get filtration values

    Macrophage_filtration_values_vertices, Macrophage_filtration_values_one_simplices, Macrophage_filtration_values_two_simplices = getFiltrationValues(Macrophage_coords, vertices, one_simplices_indices, two_simplices_indices)

    # We feed the filtration values into gudhi and get barcode

    st_Macrophage, barcode_Macrophage = get_rel_barcode(Macrophage_filtration_values_vertices,Macrophage_filtration_values_one_simplices,Macrophage_filtration_values_two_simplices,one_simplices_indices,two_simplices_indices)
    #gudhi.plot_persistence_diagram(barcode_Macrophage)
    #gudhi.plot_persistence_barcode(barcode_Macrophage)
    #plt.show()

    #Create distance vector

    bottleneck_distances = np.empty((6), float)

   
    bottleneck_distances[0] = gudhi.bottleneck_distance(st_Tumour.persistence_intervals_in_dimension(0),st_Necrotic.persistence_intervals_in_dimension(0))
    bottleneck_distances[1] = gudhi.bottleneck_distance(st_Tumour.persistence_intervals_in_dimension(1),st_Necrotic.persistence_intervals_in_dimension(1))
    bottleneck_distances[2] = gudhi.bottleneck_distance(st_Tumour.persistence_intervals_in_dimension(0),st_Macrophage.persistence_intervals_in_dimension(0))
    bottleneck_distances[3] = gudhi.bottleneck_distance(st_Tumour.persistence_intervals_in_dimension(1),st_Macrophage.persistence_intervals_in_dimension(1))
    bottleneck_distances[4] = gudhi.bottleneck_distance(st_Necrotic.persistence_intervals_in_dimension(0),st_Macrophage.persistence_intervals_in_dimension(0))
    bottleneck_distances[5] = gudhi.bottleneck_distance(st_Necrotic.persistence_intervals_in_dimension(1),st_Macrophage.persistence_intervals_in_dimension(1))


    wasserstein_distances = np.empty((6), float)
    
    # #wasserstein_distances[0] = gudhi.representations.metrics.pairwise_persistence_diagram_distances(st_stroma.persistence_intervals_in_dimension(0),st_Tumour.persistence_intervals_in_dimension(0), metric='wasserstein', n_jobs=None,order=1, internal_p=2, mode='hera', delta=0.01)
    # #gudhi.representations.metrics.WassersteinDistance(st_stroma.persistence_intervals_in_dimension(0),st_Tumour.persistence_intervals_in_dimension(0),order=1, internal_p=2, mode='hera', delta=0.01)

    #gudhi wasserstein produces an error if there is only an infinite bar in the dim 0 barcode

    if len(st_Tumour.persistence_intervals_in_dimension(0)) == 1 & len(st_Necrotic.persistence_intervals_in_dimension(0))== 1:
        wasserstein_distances[0] = abs(st_Tumour.persistence_intervals_in_dimension(0)[0,0]-st_Necrotic.persistence_intervals_in_dimension(0)[0,0])
    else:   
        wasserstein_distances[0] = gudhi.wasserstein.wasserstein_distance(st_Tumour.persistence_intervals_in_dimension(0),st_Necrotic.persistence_intervals_in_dimension(0), order=1., internal_p=2.)
    
    
    wasserstein_distances[1] = gudhi.wasserstein.wasserstein_distance(st_Tumour.persistence_intervals_in_dimension(1),st_Necrotic.persistence_intervals_in_dimension(1), order=1., internal_p=2.)
    
    if len(st_Tumour.persistence_intervals_in_dimension(0)) == 1 & len(st_Macrophage.persistence_intervals_in_dimension(0)) == 1:
        wasserstein_distances[2] = abs(st_Tumour.persistence_intervals_in_dimension(0)[0,0] - st_Macrophage.persistence_intervals_in_dimension(0)[0,0])
    else:
        wasserstein_distances[2] = gudhi.wasserstein.wasserstein_distance(st_Tumour.persistence_intervals_in_dimension(0),st_Macrophage.persistence_intervals_in_dimension(0), order=1., internal_p=2.)
    
    wasserstein_distances[3] = gudhi.wasserstein.wasserstein_distance(st_Tumour.persistence_intervals_in_dimension(1),st_Macrophage.persistence_intervals_in_dimension(1), order=1., internal_p=2.)
    
    if len(st_Necrotic.persistence_intervals_in_dimension(0)) == 1 & len(st_Macrophage.persistence_intervals_in_dimension(0)) == 1:
        wasserstein_distances[4] = abs(st_Necrotic.persistence_intervals_in_dimension(0)[0,0] - st_Macrophage.persistence_intervals_in_dimension(0)[0,0])
    else:
        wasserstein_distances[4] = gudhi.wasserstein.wasserstein_distance(st_Necrotic.persistence_intervals_in_dimension(0),st_Macrophage.persistence_intervals_in_dimension(0), order=1., internal_p=2.)
    
    wasserstein_distances[5] = gudhi.wasserstein.wasserstein_distance(st_Necrotic.persistence_intervals_in_dimension(1),st_Macrophage.persistence_intervals_in_dimension(1), order=1., internal_p=2.)


    return bottleneck_distances, wasserstein_distances


# Script

path = '/Users/stolz/Desktop/Multiplex Project/JoshData/17082022_all2Params_t500/'


param_info = pd.read_csv('/Users/stolz/Desktop/Multiplex Project/params_file.csv') 


    

all_bottleneck_distance_vectors = np.empty((0,6),float)
all_wasserstein_distance_vectors = np.empty((0,6),float)
both_distance_vectors = np.empty((0,12),float)

chi_macrophageToCSF= np.empty((0),float) 
halfMaximalExtravasationCsf1Conc = np.empty((0),float)	
Timepoint= np.empty((0),float)


for filenr in range(1,1621):
    
    file = "ID-{}_time-500_From2ParamSweep_Data.csv".format(filenr)
    
    if (os.path.exists(path+file) == 0) or (filenr == 1223):
        
        continue
    
    data = pd.read_csv(path+file) 
    bottleneck_distances, wasserstein_distances = get_rel_distance_vector(data)
    
    chi_macrophageToCSF = np.append(chi_macrophageToCSF,[param_info.iloc[filenr-1]['chi_macrophageToCSF'].round(decimals=1)],axis = 0)
    halfMaximalExtravasationCsf1Conc= np.append(halfMaximalExtravasationCsf1Conc,[param_info.iloc[filenr-1]['halfMaximalExtravasationCsf1Conc'].round(decimals=1)],axis = 0)
    #Timepoint= np.append(Timepoint,[param_info.iloc[filenr-1]['Timepoint']],axis = 0)
    
    #bottleneck_distances = get_distance_vector(data)
    
    print(file)
    print(bottleneck_distances)
    print(wasserstein_distances)

    all_bottleneck_distance_vectors = np.append(all_bottleneck_distance_vectors,[bottleneck_distances],axis = 0)
    all_wasserstein_distance_vectors = np.append(all_wasserstein_distance_vectors,[wasserstein_distances],axis = 0)
    both_distance_vectors = np.append(both_distance_vectors,[np.append(bottleneck_distances,wasserstein_distances)],axis = 0)


kmeans_bottleneck = KMeans(n_clusters=3, random_state=0).fit(all_bottleneck_distance_vectors)
labels_bottleneck = kmeans_bottleneck.labels_


kmeans_wasserstein = KMeans(n_clusters=3, random_state=0).fit(all_wasserstein_distance_vectors)
labels_wasserstein = kmeans_wasserstein.labels_

kmeans_both = KMeans(n_clusters=3, random_state=0).fit(both_distance_vectors)
labels_both = kmeans_both.labels_


#bottleneck

show_labels = np.empty((0,3),float)
show_purity = np.empty((0,3),float)

for chi in [0.5,1,1.5,2,2.5,3,3.5,4,4.5]:
    
    for Csf1 in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
    
        indices = np.where((chi_macrophageToCSF == chi) & (halfMaximalExtravasationCsf1Conc == Csf1))
        
        if np.array(indices).size == 0:
            continue
        
        z = np.bincount(labels_bottleneck[indices]).argmax()
        purity = len(np.where(labels_bottleneck[indices] == z)[0])/len(labels_bottleneck[indices])


        show_labels = np.append(show_labels,[[Csf1, chi, z]], axis = 0)
        show_purity = np.append(show_purity,[[Csf1, chi, purity]], axis = 0)


mapping= {0: ("blue", "o"), 1: ("red", "o"), 2: ("yellow","o")}#, 3: ("white","o")}

for c in np.unique(show_labels[:,2]):
    d = show_labels[show_labels[:,2] == c]
    plt.scatter(d[:,0], d[:,1], s=400, c=mapping[c][0], marker=mapping[c][1])

plt.title('Bottleneck clusters')
plt.xlabel('Csf1')
plt.ylabel('chi')
plt.show()

plt.scatter(show_purity[:,0], show_purity[:,1], s=400, c=show_purity[:,2], marker="o", vmin=0, vmax=1)
plt.title('Purity Bottleneck clusters')
plt.xlabel('Csf1')
plt.ylabel('chi')
plt.colorbar()
plt.show()








show_labels = np.empty((0,3),float)
show_purity = np.empty((0,3),float)

for chi in [0.5,1,1.5,2,2.5,3,3.5,4,4.5]:
    
    for Csf1 in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
    
        indices = np.where((chi_macrophageToCSF == chi) & (halfMaximalExtravasationCsf1Conc == Csf1))
        
        if np.array(indices).size == 0:
            continue
        
        z = np.bincount(labels_wasserstein[indices]).argmax()
        
        purity = len(np.where(labels_wasserstein[indices] == z)[0])/len(labels_wasserstein[indices])


        show_labels = np.append(show_labels,[[Csf1, chi, z]], axis = 0)
        show_purity = np.append(show_purity,[[Csf1, chi, purity]], axis = 0)


mapping= {0: ("blue", "o"), 1: ("red", "o"), 2: ("yellow","o")}

for c in np.unique(show_labels[:,2]):
    d = show_labels[show_labels[:,2] == c]
    plt.scatter(d[:,0], d[:,1], s=400, c=mapping[c][0], marker=mapping[c][1])

plt.title('Wasserstein clusters')
plt.xlabel('Csf1')
plt.ylabel('chi')
plt.show()

plt.scatter(show_purity[:,0], show_purity[:,1], s=400, c=show_purity[:,2], marker="o", vmin=0, vmax=1)
plt.title('Purity Wasserstein clusters')
plt.xlabel('Csf1')
plt.ylabel('chi')
plt.colorbar()
plt.show()




show_labels = np.empty((0,3),float)
show_purity = np.empty((0,3),float)


for chi in [0.5,1,1.5,2,2.5,3,3.5,4,4.5]:
    
    for Csf1 in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
    
        indices = np.where((chi_macrophageToCSF == chi) & (halfMaximalExtravasationCsf1Conc == Csf1))
        
        if np.array(indices).size == 0:
            continue
        
        z = np.bincount(labels_both[indices]).argmax()
        purity = len(np.where(labels_both[indices] == z)[0])/len(labels_both[indices])


        show_labels = np.append(show_labels,[[Csf1, chi, z]], axis = 0)
        show_purity = np.append(show_purity,[[Csf1, chi, purity]], axis = 0)


mapping= {0: ("blue", "o"), 1: ("yellow", "o"), 2: ("red","o")}

for c in np.unique(show_labels[:,2]):
    d = show_labels[show_labels[:,2] == c]
    plt.scatter(d[:,0], d[:,1], s=400, c=mapping[c][0], marker=mapping[c][1])

plt.title('Bottleneck + Wasserstein')
plt.xlabel('Csf1')
plt.ylabel('chi')
plt.show()

plt.scatter(show_purity[:,0], show_purity[:,1], s=400, c=show_purity[:,2], marker="o", vmin=0, vmax=1)
plt.title('Purity Bottleneck + Wasserstein clusters')
plt.xlabel('Csf1')
plt.ylabel('chi')
plt.colorbar()
plt.show()
