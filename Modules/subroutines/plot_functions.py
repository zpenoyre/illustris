# Import libaries and illustris library
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import random
from matplotlib.patches import Polygon
from matplotlib.collections import PolyCollection

def pixelMap(r,m,nPixels,boxSize):
	
	pixels=np.zeros((nPixels,nPixels))
	pixSize=boxSize/nPixels
	
	for i in range(nPixels):
		for j in range(nPixels):
			
			inPix=(r[:,1]>pixSize*i) & (r[:,1]<pixSize*(i+1)) & (r[:,0]>pixSize*j) & (r[:,0]<pixSize*(j+1))
			
			if inPix.size==0:
				continue
			
			pixels[i,j]=np.sum(m[inPix])
	
	return pixels
