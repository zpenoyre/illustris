import numpy as np
#from numba import jit

# input: start AND end number (we have z=0 == 0 and then go up as we go up in redshift), AND index of subfind
# output: index of subfind at end number
#@jit
def find_galaxy(data, start_no, end_no, start_index):
    # check if we want z=0, the jump right there
    if end_no == 0:
        if data[start_no][start_index,21] != -1:
            if np.argwhere(data[0][:,1] == data[start_no][start_index,21]).size != 0:
                index = np.argwhere( data[0][:,1] == data[start_no][start_index,21] )[0][0];
            else:
                index = -1;
        else:
            index = -1;
    else:
        index = start_index;
        
        # prepare to walk forward in time
        for i in range(start_no, end_no, -1):
            #select next subhalo number
            next_sub = data[i][index,22]
            
            # check if next subhalo available
            if next_sub != -1:
                # select index of next subhalo
                index = np.argwhere( data[i-1][:,1] == next_sub )[0][0];
            else:
                # check wether successor subhalo exists at z=0, if so walk from z=0 backward in time
                if data[i][index,21] != -1:
                    for j in range (0,end_no):
                        prev_sub = data[j][index,23]
                        # check if previous subhalo available
                        if prev_sub != -1:
                            # select index of previous subhalo
                            index = np.argwhere( data[i+1][:,1] == prev_sub )[0][0];
                        else:
                            index = -1
                else:
                    index = -1
                    
    return index

# input :0-3 for z begin and end, z=0:0
# output: column 0: location in web at start, column 1: location in web at end 2: index of subfind at z_end 3: index of subfind at z_start
def web_evolution(data, webdata, z_begin, z_end):
        
    number_map = [0,19,32,50]; start_no = number_map[z_begin]; end_no =  number_map[z_end]  
    a_size= data[start_no][:,0].shape[0];
    output = np.zeros((a_size,4));
    
    for i in range(0, a_size):
        index = find_galaxy(data, start_no, end_no, i);
        output[i, 0] = webdata[z_begin][i,17];
        if index != -1:
            output[i, 1] = webdata[z_end][index,17];
            output[i, 2] = index
            output[i, 3] = i
        else:
            output[i, 1] = np.NAN
            output[i, 2] = np.NAN
    return output

