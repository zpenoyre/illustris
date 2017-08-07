import numpy as np


def time_integration(data, start_no, back_steps_no, integration_ind):

    a_size= data[start_no][:,0].shape[0]
    output = np.zeros((a_size,2))
    times = np.loadtxt('../../Data/time_conversion.txt')

    # for every entry in the gal data
    for i in range (0, a_size):

        # define matrix to hold data of all back steps
        temp_a = np.zeros((2,back_steps_no))

        # go back all steps
        for j in range (0, back_steps_no):

            # if we reach 50, just take the next 6, not previous 6
            if (start_no+back_steps_no)<=50:
                index = find_galaxy_backward(data, start_no, start_no+(j+1),i)
                
                # record data in temp array
                if index != -1:
                    temp_a[0,j] = data[start_no+j+1][index, integration_ind]
                    temp = 0

                # if we encounter -1  record where we did and break out of backsteps
                else:
                    temp = j+1
                    continue

            else:
                index = find_galaxy(data, start_no, start_no-(j+1),i)

                # record data in temp array
                if index != -1:
                    temp_a[0,j] = data[start_no-j-1][index, integration_ind]
                    temp = 0

                # if we encounter -1  record where we did and break out of backsteps
                else:
                    temp = j+1
                    continue
            

            # integration
            # if we could walk all steps to trapezoid integration with gigayears and the datapoint
            # if not, then do the same integration up to the step we could walk back and then scale it acccordingly
            # if we could only step back one give out nan 
        if j == 0:
            times_ind = 133-start_no+1 

            int_value = sum((times[times_ind-back_steps_no+1:times_ind,2]-times[times_ind-back_steps_no:times_ind-1,2]) * (temp_a[0,0:back_steps_no-1]+temp_a[0,1:back_steps_no])/2)

            output[i, 0] = int_value
            output[i, 1] = j
        
        elif j == 1:
            output[i, 0] = np.NAN
        
        else:
            times_ind = 133-start_no+1
            int_value = sum( ( (times[times_ind-j+1:times_ind:,2]-times[times_ind-j:times_ind-1:,2]) * (temp_a[0,0:j-1]+temp_a[0,1:j])/2) * back_steps_no/j)


            output[i, 0] = int_value
            output[i, 1] = j
        
        if (i%int(a_size/5) == 0):
            print ("Integration done with", (i/int(a_size/5))*20,"%")
    return output





# input: start AND end number (we have z=0 == 0 and then go up as we go up in redshift), AND index of subfind
# output: index of subfind at end number
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
                index = np.argwhere( data[i-1][:,1] == next_sub )
                if index.size != 0:
                    index = index[0][0];
                else:
                    index = -1
            else:
                index = -1
                # # check wether successor subhalo exists at z=0, if so walk from z=0 backward in time
                # if data[i][index,21] != -1:
                #     for j in range (0,end_no):
                #         prev_sub = data[j][index,23]
                #         # check if previous subhalo available
                #         if prev_sub != -1:
                #             # select index of previous subhalo
                #             index = np.argwhere( data[i+1][:,1] == prev_sub )[0][0];
                #         else:
                #             index = -1
                # else:
                #     index = -1
                    
    return index


def find_all_data_evolution(data, webdata, start_no, end_no, start_index):
    number_map = [0,19,32,50];
    output = np.zeros((end_no-start_no+1,25))
    output[0,0:21] = data[start_no][start_index,0:21]
    output[0,21] = webdata[0][start_index,17]

    index = start_index
    for i in range(start_no, end_no):
        index  = find_galaxy_backward(data, i, i+1, index);
        if index == -1:
            output[i+1,:] = np.NAN
            continue
        else: 
            output[i+1,0:21] = data[i+1][index,0:21]
            
            if i+1 == 19:
                output[0,22] = webdata[1][index,17]
            if i+1 == 32:
                output[0,23] = webdata[2][index,17]
            if i+1 == 50:
                output[0,24] = webdata[3][index,17]

    return output


# input: start AND end number (we have z=0 == 0 and then go up as we go up in redshift), AND index of subfind
# output: index of subfind at end number
def find_galaxy_backward(data, start_no, end_no, start_index):
    
    index = start_index;
        
    # prepare to walk backward in time
    for i in range(start_no, end_no):
        
        if index == -1:
            continue
        #select next subhalo number
        prev_sub = data[i][index,23]
        
        # check if next subhalo available
        if prev_sub != -1:
            # select index of next subhalo
            index = np.argwhere( data[i+1][:,1] == prev_sub)
            if index.size != 0:
                index = index[0][0];
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

