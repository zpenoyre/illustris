import numpy as np


# input : array with mass data , bins: no of bins in log space
# output: array: 0: median of mass in interval, 1: mass interval step size (dm) 2: dn/dm
def mass_function(data_input, bins):

    binned_mat = np.zeros((bins, 3))
    data = data_input

    max_exp = np.log10(np.amax(data))
    min_exp = np.log10(np.amin(data))
    step_array = np.logspace(min_exp, max_exp, bins)
         
    for i in range (0, bins):
        
        step = step_array[i]
        if (i == bins-1): 
            step1 = np.amax(data)+1
        else:
            step1 = step_array[i+1]
            
        mask = ( (data > step) & (data < step1) )
        
        data_mask = data[mask]

        binned_mat[i,0] = np.median(data_mask)
        binned_mat[i,1] = (step1 + step)/2
        binned_mat[i,2] = np.count_nonzero(data_mask)/((step1-step)*10**10)
        
    return binned_mat

# Logarithmic base 10 Binning
# input: data matrix, index1-2: indices of data to be binned (binned w/ respect to index 1), no of bins, if bmedian =1 median 
# output: 0-1: binned data from index, 2-3: error from indicies, 
#         4: count of data points in bin devided by step size, 5: count of data points in bin
def log_bin_data(bin_data, index1, index2, bins, bmedian):
    
    binned_mat = np.zeros((bins, 6))
    data = bin_data[:, [index1, index2]]

    max_exp = np.log10(np.amax(data[:,0]))
    min_exp = np.log10(np.amin(data[:,0]))
    step_array = np.logspace(min_exp, max_exp, bins)
         
    for i in range (0, bins):
        
        step = step_array[i]
        if (i == bins-1): 
            step1 = np.amax(data[:,0])+1
        else:
            step1 = step_array[i+1]
            
        mask = ( (data[:,0] > step) & (data[:,0] < step1) )
        
        data_mask = data[mask]
        
        if (bmedian == 1):
            binned_mat[i,0] = np.median(data_mask[:,0])
            binned_mat[i,1] = np.median(data_mask[:,1])
            
        else:
            binned_mat[i,0] = np.mean(data_mask[:,0])
            binned_mat[i,1] = np.mean(data_mask[:,1])
        binned_mat[i,2] = (np.std(data_mask[:,0]))/(np.sqrt(np.count_nonzero(data_mask[:,0])))
        binned_mat[i,3] = (np.std(data_mask[:,1]))/(np.sqrt(np.count_nonzero(data_mask[:,1])))
            
        binned_mat[i,4] = np.count_nonzero(data_mask[:,0])/(step1-step)
        binned_mat[i,5] = np.count_nonzero(data_mask[:,1])
        
    return binned_mat

# linear Binning
# input: data matrix, index1-2: indices of data to be binned (binned w/ respect to index 1), no of bins, if bmedian =1 median 
# output: 0-1: binned data from index, 2-3: error from indicies, 
#         4: count of data points in bin devided by step size, 5: count of data points in bin
def bin_data(bin_data, index1, index2, bins, bmedian):
    
    binned_mat = np.zeros((bins, 6))
    data = bin_data[:, [index1, index2]]

    maxb = np.amax(data[:,0])
    minb = np.amin(data[:,0])
    step_array = np.linspace(minb, maxb, bins)
         
    for i in range (0, bins):
        
        step = step_array[i]
        if (i == bins-1): 
            step1 = np.amax(data[:,0])+1
        else:
            step1 = step_array[i+1]
            
        mask = ( (data[:,0] > step) & (data[:,0] < step1) )
        
        data_mask = data[mask]
        
        if (bmedian == 1):
            binned_mat[i,0] = np.median(data_mask[:,0])
            binned_mat[i,1] = np.median(data_mask[:,1])
            
        else:
            binned_mat[i,0] = np.mean(data_mask[:,0])
            binned_mat[i,1] = np.mean(data_mask[:,1])
        binned_mat[i,2] = (np.std(data_mask[:,0]))/(np.sqrt(np.count_nonzero(data_mask[:,0])))
        binned_mat[i,3] = (np.std(data_mask[:,1]))/(np.sqrt(np.count_nonzero(data_mask[:,1])))
            
        binned_mat[i,4] = np.count_nonzero(data_mask[:,0])/(step1-step)
        binned_mat[i,5] = np.count_nonzero(data_mask[:,1])
        
    return binned_mat

