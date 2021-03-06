import requests
import numpy as np
import h5py

baseUrl = 'http://www.illustris-project.org/api/'
headers = {"api-key":"cc4ff6392e79c9e08c158e5ae5493718"}
    
# Routine to pull data from online
def get(path, params=None, fName='temp'): # gets data from url, saves to file
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)
    
    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    dataFile=fName+'.hdf5'
    # Saves to file, currently disabled
    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(dataFile, 'wb') as f:
            f.write(r.content)
        return dataFile # return the filename string

    return r
    
# For a chosen galaxy pulls out all the particle data for a set of fields
particleTypeNames=['gas','dm','error','tracers','stars','bhs'] 
def getGalaxy(whichGalaxy, fields, # index of a galaxy and the 2d list of fields (particle type and name of fields)
        simulation='Illustris-1', snapshot=135, # which simulation and snapshot
        fileName='temp',rewriteFile=1): # name of the file where .hdf5 data is stored and whether to rewrite or just read
    
    fields=np.array(fields) # converts to array
    order=np.argsort(fields[:,0])
    disorder=np.argsort(order) # needed to unsort the fields later...
    fields=fields[order,:] # orders by particle type
    nFields=order.size
    
    if rewriteFile==1: # redownloads file from the internet
        url='http://www.illustris-project.org/api/'+simulation+'/snapshots/'+str(snapshot)+'/subhalos/'+str(whichGalaxy)+'/cutout.hdf5?'
    
        thisParticle=0
        thisEntry=0
        firstParticle=1
        while thisParticle<6: # cycles through all particle type

            if (int(fields[thisEntry,0])!=thisParticle): # checks there is at least one field for this particle
                thisParticle+=1
                continue
            if firstParticle==1: # first particle requires no ampersand
                firstParticle=0
            else: # all later particles do
                url+='&' 
            url+=particleTypeNames[thisParticle]+'=' # adds the name of the particle type
        
            firstEntry=1
            while int(fields[thisEntry,0])==thisParticle:
                if firstEntry==1: #first entry requires no comma
                    firstEntry=0
                else: # all later entries do
                    url+=','
                url+=fields[thisEntry,1] # adds every associated field
                thisEntry+=1
                if thisEntry==nFields:
                    break
            if thisEntry==nFields:
                break
            thisParticle+=1
        dataFile=get(url,fName=fileName)
    # end of "if rewriteFile==1:"
    if rewriteFile == 0: # if we're not redownloading need to set path to the file
        dataFile=fileName+'.hdf5'
    
    # actually get the data (saved to .hdf5 file)
    data=[] # initially empty list that we will fill up with the data
    with h5py.File(dataFile,'r') as f:
        for i in range(disorder.size):
            thisField=fields[disorder[i],:] # ensures data returned in original order of fields
            data.append(np.array(f['PartType'+thisField[0]][thisField[1]]))
            # returns all particle data of each field as a numpy array
    return data # returns all the particle fields as a list of numpy arrays in the same order as initial fields

def getSubhaloField(field,simulation='Illustris-1',snapshot=135,fileName='temp',rewriteFile=1):
    if rewriteFile==1: # redownloads file from the internet
        url='http://www.illustris-project.org/api/'+simulation+'/files/groupcat-'+str(snapshot)+'/?Subhalo='+field
        dataFile=get(url,fName=fileName)
    if rewriteFile == 0: # if we're not redownloading need to set path to the file
        dataFile=fileName+'.hdf5'
        
    with h5py.File(dataFile,'r') as f:
        data=np.array(f['Subhalo'][field])
    return data
    
def getHaloField(field,simulation='Illustris-1',snapshot=135,fileName='temp',rewriteFile=1):
    if rewriteFile==1: # redownloads file from the internet
        url='http://www.illustris-project.org/api/'+simulation+'/files/groupcat-'+str(snapshot)+'/?Group='+field
        dataFile=get(url,fName=fileName)
    if rewriteFile == 0: # if we're not redownloading need to set path to the file
        dataFile=fileName+'.hdf5'
        
    with h5py.File(dataFile,'r') as f:
        data=np.array(f['Group'][field])
    return data

def getSim(simName):
    r = get(baseUrl)
    names = [sim['name'] for sim in r['simulations']]
    i = names.index(simName)
    sim = get( r['simulations'][i]['url'] )
    return sim
    
def getSnap(sim,whichSnap):
    snaps = get( sim['snapshots'] )
    snap = get( snaps[whichSnap]['url'] )
    return snap
    
def getSub(subs,whichSub):
    sub_url=subs['results'][whichSub]['url']
    sub=get(sub_url)
    return sub
