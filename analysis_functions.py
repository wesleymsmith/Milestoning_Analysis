import numpy as np
import pandas as pd
import itertools

def analyze_milestone1D_data(dataTable,windowMins,windowMaxs,verbose=False):
    winMins=np.array(windowMins)
    winMaxs=np.array(windowMaxs)
    windowCenters=(windowMins+windowMaxs)/2.0
    simData=dataTable[['Window','Time','X']]
    binEdges=winMins[1:] #np.concatenate([winMins,[winMaxs[-1]]])
    digitize_kwds={"bins":binEdges}
    simData['X_Index']=simData.X.apply(np.digitize,**digitize_kwds)
    
    windows=np.sort(simData.Window.unique())
    xbins=np.sort(simData.X_Index.unique())
    nBins=len(xbins)
    escapeMat=np.zeros((nBins,nBins))
    rVec=np.zeros(nBins-1)
    tSum=0
    #iVal->escape matrix row index
    #xbin->window
    #cVal->place holder for bin index with indexing starting at 1
    for iVal,xbin in enumerate(xbins):
        if xbin in windows:
            tempDat=simData[simData.Window==xbin]
            cVal=xbin+1
            #binVec tracks the index of the current bin. Bin indices
            #are numbered starting at 1 so that we can use values of 0
            #in boolean calculations
            binVec=np.array(tempDat.X_Index+1)
            #binC is equal to 1 if the coordinate is in the bin of the current window
            #and zero if not
            binC=(binVec==cVal)
            #binT tracks which frames a transition has occured in
            #for frames where the coordinate is not in the current window's bin (1-binC[1:])
            #but was in the previous frame (binC[:-1]) list the index of the 
            #current coordinate bin (binVec[1:]) all other frames are 0
            binT=(1-binC[1:])*binC[:-1]*binVec[1:]
            #tCounts returns the 2D array. The first column is a sorted list of all
            #unique values observed in binT. The second column is how many times
            #each of those values occured
            tCounts=np.unique(binT,return_counts=True)
            transInds=tCounts[0][1:] #first entry should always be for binT==0
            transCounts=tCounts[1][1:]
            cCount=np.sum(binC)
            #generate escape matrix
            for iInd,Ind in enumerate(transInds):
                escapeMat[iVal,Ind-1]=1.*transCounts[iInd]/cCount
            #need to exclude frames where coordinate has not just transistioned and is not
            #in window bin from binT
            runVec=binT[np.nonzero(binC[1:]+binT)] 
            #use itertools to get 'runs' in runVec listing the value of the repeated number in the run
            #and the length of the run as a 2D array
            #I.e. runList[:,0] gives the values of the repeated numbers themselves
            #runList[:,1] gives the length of the corresponding runs
            runList=np.array([[int(j[0]),len(j)] for j in \
                              [list(g) for k,g in itertools.groupby(runVec)]])
            #generate R vector
            for iRun,run in enumerate(runList[1:]):
                if run[0]==0:
                    #the bin edge index between bins i and j is min(i,j)
                    rVec[np.min([runList[iRun,0]-1,xbin])]+=run[1]
            tSum=tSum+np.sum(binC)
            if verbose:
                print "--- --- ---"
                print "escapeMatrix entry for window %g:"%xbin
                print '['+', '.join(map(lambda x: '%.5f'%x,escapeMat[iVal,:]))+']'
    if verbose:
        print "--- --- ---"
        
    Rmat=np.matrix(escapeMat)
    Dmat=np.matrix(np.diag(1-np.sum(escapeMat,axis=1)))

    RDmat=Rmat+Dmat

    outEig=np.linalg.eig(RDmat.T)
    si=np.argsort(1-outEig[0])
    if verbose:
        print 'Eigenvalues:',
        print outEig[0][si]
    outVec=np.array(outEig[1])[:,si[0]]
    outVec=outVec/np.sum(outVec)
    
    return (escapeMat,outVec,rVec/tSum)
