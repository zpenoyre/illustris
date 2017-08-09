import numpy as np

# input: data matrix (1 column), data for x axis (binned with respect to), list of bins, the value for mask to be binned (in this case 0 1 or -1 for cluster/fil), the max for the linspace
# output: linspace, the fraction and the error term
def bin_frac(data, data_bin, bins, fraction_value):
    
    output= np.zeros((len(bins)-1, 4))

    for i in range (0, len(bins)-1):

        mask = data[(data_bin[:]>bins[i]) & (data_bin[:]<bins[i+1])]
        
        if mask.size != 0:
            fraction = (mask[:] == fraction_value).sum()/mask.size
            error = 1/np.sqrt( mask.size)
        else:
            fraction = np.NaN
            error = 0
        
        output[i,0] = np.mean(bins[i:i+2])
        output[i,1] = fraction
        output[i,2] = error
        output[i,3] = (mask[:] == fraction_value).sum()

    return output

def find_frac(data, fraction_value):

    output= np.zeros((2))

    fraction = (data[:]==fraction_value).sum()/data.size
    error = 1/np.sqrt((data[:]==fraction_value).sum())
    
    output[0] = fraction
    output[1] = error

    return output

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

#  Binning
# input: data matrix 1 and 2, bins is the bin list, if bmedian =1 median 
# output: 0-1: binned data from index, 2-3: error from indicies, 
#         4: count of data points in bin devided by step size, 5: count of data points in bin
def bin_data(bin_data_1, bin_data_2, bins, bmedian = 0):
    
    binned_mat = np.zeros((len(bins)-1, 7))
     
    for i in range (0, len(bins)-1):
            
        mask = ( (bin_data_1[:] > bins[i]) & (bin_data_1[:] < bins[i+1]) )
        
        data_mask_1 = bin_data_1[mask]
        data_mask_2 = bin_data_2[mask]

        if data_mask_1.size != 0:

            if (bmedian == 1):
                binned_mat[i,0] = np.median(data_mask_1)
                binned_mat[i,1] = np.median(data_mask_2)
                binned_mat[i,6] = np.median(bins[i:i+2])
                
            else:
                binned_mat[i,0] = np.mean(data_mask_1)
                binned_mat[i,1] = np.mean(data_mask_2)
                binned_mat[i,6] = np.mean(bins[i:i+2])
            binned_mat[i,2] = (np.std(data_mask_1))/(np.sqrt(data_mask_1.size))
            binned_mat[i,3] = (np.std(data_mask_2))/(np.sqrt(data_mask_2.size))
                
            binned_mat[i,4] = (data_mask_1.size)/(bins[i+1]-bins[i])
            binned_mat[i,5] = data_mask_2.size
        
    return binned_mat