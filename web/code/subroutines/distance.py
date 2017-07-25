# Import libaries and illustris library
import numpy as np
import math

# calulate distance to a point:

# input: subhalo_mask: position of particles in column matrix x,y,z position
#        fil_pos_mat: position of filaments in column matrix x,y,z position

# output: array: 1st column is index of filament row that is nearest to particle. 2nd column is absolute distance

def distance_function(subhalo_mask, fil_pos_mat):
    
    #initialize number of columns in input file and number of steps to print
    subhalo_mask_size = subhalo_mask.shape[0]
    steps = int(max((subhalo_mask_size)/10,1))
    
    # create distance matrix - 2 columns for data 
    distance_mat = np.zeros(((subhalo_mask_size), 2))
    
    # enter for-loop which goes over the size of the rows of the position matrix 
    for i in range(0, subhalo_mask_size):
        
        #calculate actual distance
        distances = np.linalg.norm(fil_pos_mat - subhalo_mask[i,0:3] ,axis=1)
        
        #find minimum - i.e. nearest filament to particle
        min_index = np.argmin(distances, axis=0)

        #for each particle give index of filament point that is nearest and it´s absolute distance
        distance_mat[i,0] = min_index     
        distance_mat[i,1] = distances[min_index]
        
        if (i%steps == 0):
            print ("distance_function done with", (i/steps) *10,"%")

    # return the distance matrix
    return distance_mat


# calulate aligment of spin vector:

# input ang_mom_mat: x,y,z of angular momentum vector
#       distance_mat: index of nearest segment point, distance
#       fil_pos_mat: column matrix x,y,z position

# output: array: 1st column is index of filament row that is nearest to particle. 2nd column is absolute distance, 
#                                              3rd is cos(theta)
#                                              4th is sin(theta)
#                                              5th is angle

def align_function(ang_mom_mat, distance_mat, fil_pos_mat):
    
    #initialize number of columns in input file and number of steps to print
    ang_mom_mat_size = ang_mom_mat.shape[0]
    steps = int(max((ang_mom_mat_size)/4,1))
    
    #create align matrix - 5 columns for data 
    align_mat = np.zeros(((ang_mom_mat_size), 5))
    
    # enter for-loop which goes over the size of the rows of the position matrix 
    for i in range(0, ang_mom_mat_size):
        
        #calculate vector segment is pointing in (0 from index in distance_mat to -1)
        fil_vec = fil_pos_mat[(int(distance_mat[i,0])),0:3] - fil_pos_mat[(int(distance_mat[i,0]) -1),0:3]
        
        #get angular mom. vector of galaxy
        gal_ang_mom = ang_mom_mat[i,:]
        
        #calculate size
        fil_vec_size = np.linalg.norm(fil_vec)
        gal_ang_mom_size = np.linalg.norm(gal_ang_mom)
        
        align_mat[i,0] = distance_mat[i,0]
        align_mat[i,1] = distance_mat[i,1]
        
        align_mat[i,2] = np.dot(fil_vec, gal_ang_mom)/(fil_vec_size * gal_ang_mom_size)
        align_mat[i,3] = np.linalg.norm( np.cross(fil_vec, gal_ang_mom) )/(fil_vec_size * gal_ang_mom_size)
        align_mat[i,4] = np.angle( align_mat[i,2] + align_mat[i,3]*1.j)
        
        if (i%steps == 0):
            print ("align_function done with", (i/steps) *25,"%")

    # return the distance matrix 
    return align_mat