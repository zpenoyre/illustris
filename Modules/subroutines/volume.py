import numpy as np
import math

# mask for any purpose, index is the first of 3 indices to use for cuts
def mask(array, index, cut0, cut1, cut2):
	mask = (
	(array[:,index] > 0) & (array[:,index] < cut0) & 
	(array[:,index+1] > 0) & (array[:,index+1] < cut1) & 
	(array[:,index+2] > 0) & (array[:,index+2] < cut2)
	)

	return array[mask]

# output: array: absolute distance, if below 2000kpc otherwise just -1
def distance_function(positions_mat, fil_pos_mat):
    
    #initialize number of columns in input file and number of steps to print
    positions_mat_size = positions_mat.shape[0]
    steps = int(max((positions_mat_size)/10,1))
    
    # create distance matrix - 2 columns for data 
    distance_mat = np.zeros(((positions_mat_size), 1))
    
    # use the cutoff of 2000 for distance from filaments/clusters (this is a square, will get smaller with spherical distance)
    cut = 2000
    # enter for-loop which goes over the size of the rows of the position matrix 
    for i in range(0, positions_mat_size):
        
        filaments_mask = (
        	(fil_pos_mat[:,0]>(positions_mat[i,0]-cut)) & (fil_pos_mat[:,0]<(positions_mat[i,0]+cut)) & 
        	(fil_pos_mat[:,1]>(positions_mat[i,1]-cut)) & (fil_pos_mat[:,1]<(positions_mat[i,1]+cut)) & 
        	(fil_pos_mat[:,2]>(positions_mat[i,2]-cut)) & (fil_pos_mat[:,2]<(positions_mat[i,2]+cut)) 
        	)
        
        relevant_filaments = fil_pos_mat[filaments_mask][:,0:3]
        
        if (relevant_filaments.size != 0):
            distance_mat[i,0] = np.min(
                np.linalg.norm(
                    relevant_filaments - positions_mat[i,0:3]
                    ,axis=1)
            )

        else:
            distance_mat[i,0] = -1
            
        if (i%steps == 0):
            print("distance_function done with", (i/steps) *10,"%")

    return distance_mat

# output: array: first entry critical volume fraction, second entry filament volume fraction
def volume_fraction(critical_distance_mat, segments_distance_mat, segments_cut):

    segs_count = 0
    crits_count = 0

    # cut distances off at cutoff
    segs_mask = segments_distance_mat < segments_cut

    # count instances of distance below cutoff
    for i in range (0, critical_distance_mat.shape[0]):
        if(critical_distance_mat[i] == 1):
            crits_count = crits_count +1
        elif(segs_mask[i] == True):
            segs_count = segs_count +1

    # compute fractions            
    vol_crits = crits_count/critical_distance_mat.shape[0]
    vol_segs = segs_count/critical_distance_mat.shape[0]

    return np.array([vol_crits,vol_segs])
