import numpy as np
import scipy as sp
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
    nEdges=nBins-1
    escapeMat=np.zeros((nBins,nBins))
    rMat=np.zeros([nBins,nEdges])
    crossArray=np.zeros([nBins,2])
    tSum=0
    countsVec=np.zeros(nBins)
    #iVal->escape matrix row index
    #xbin->window
    #cVal->place holder for bin index with indexing starting at 1
    for iVal,xbin in enumerate(xbins):
        if xbin in windows:
            tempDat=simData[simData.Window==xbin]
            cVal=xbin+1
            binVec=np.array(tempDat.X_Index+1)
            binC=(binVec==cVal)
            binT=(1-binC[1:])*binC[:-1]*binVec[1:]
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
            runList=np.array([[int(j[0]),len(j)] for j in \
                              [list(g) for k,g in itertools.groupby(runVec)]])
            #generate R vector
            for iRun,run in enumerate(runList[1:]):
                if run[0]==0:
                    #the bin edge index between bins i and j is min(i,j)
                    rMat[iVal,np.min([runList[iRun,0]-1,xbin])]+=run[1]
            #tSum=tSum+np.sum(binC)
            countsVec[iVal]=np.sum(binC)
            #count number of times central cell was traversed.
            #we again make use of itertools.
            #First we remove all zero entries from binT
            #Then we compute runs of the non-zero binT. Here, however, we don't care
            #about the counts themselves. I.e. we only care about the binIDs
            #To count the transitions from the left edge to right edge or vice versa
            #we iterate over this list of binIDs and update an nBinEdges x 2 array (crossArray)
            #by adding 1 to the first column of the row for the current bin if the
            #if the current edgeID is the right edge and the previous was the right
            #edge or adding 1 to the second column if the current edgeID is left
            #edge and the previous was the right edge.
            crossings=(np.array([[int(j[0]),len(j)] for j in \
                              [list(g) for k,g in itertools.groupby(binT[np.nonzero(binT)])]]
                              )[:,0]).flatten()
            #rather than needing a for loop, we can use vector arithmetic.
            #we subract entries 0:n-1 from entries 1:n and divide by 2.
            #this will yield a -1 when a traversal from the left edge to the right edge occured
            #or a +1 when a traversal from right edge to left edge occured.
            #we then count the number of -1's and this to column 0 of the row in crossArray for
            #the current window then count the number of 1's and this to column 1 of that row.
            crossings=(crossings[1:]-crossings[:-1])/2
            crossArray[iVal,0]=crossArray[iVal,0]+np.sum(crossings==1)
            crossArray[iVal,1]=crossArray[iVal,1]+np.sum(crossings==-1)
            if verbose:
                print "--- --- ---"
                print "escapeMatrix entry for window %g:"%xbin
                print '['+', '.join(map(lambda x: '%.5f'%x,escapeMat[iVal,:]))+']'
                print "Number of crossings (left-to-right,right-to-left):",
                print "(%g,%g)"%(crossArray[iVal,0],crossArray[iVal,1])
    if verbose:
        print "--- --- ---"
        
    tSum=np.sum(countsVec)
    Emat=np.matrix(escapeMat)
    Dmat=np.matrix(np.diag(1-np.sum(escapeMat,axis=1)))

    Amat=Emat+Dmat

    outEig=np.linalg.eig(Amat.T)
    si=np.argsort(1-outEig[0])
    if verbose:
        print 'Eigenvalues:',
        print outEig[0][si]
    piVec=np.array(outEig[1])[:,si[0]]
    piVec=piVec/np.sum(piVec)
    
    #since we stored Rij in matrix form (nBins x nEdges) we can
    #use matrix and array arithmetic to compute Ri
    #note that rMat.T will yield an nEdges x nBins matrix,
    #so summing over columns (i.e. computing row sums) will
    #give us the needed values for each edge
    #E.g. crossArray[ii,0] = transition from edge ii-1 -> ii
    #     crossArray[ii,1] = transition from edge ii -> ii-1
    Ri=1.*np.sum(rMat.T*piVec/countsVec,axis=1)
    NijMat=np.zeros([nBins-1,nBins-1])
    for ii in np.arange(1,nEdges):
        #note that indexing in python starts at 0 (not 1)
        #ii loops upward over bins omitting first and last bin
        #left to right transitions...
        NijMat[ii-1,ii]=piVec[ii]*crossArray[ii,0]/countsVec[ii]
        #right to left transitions
        #jj loops downward over bins ommiting last and first bin
        jj=nEdges-ii
        NijMat[jj,jj-1]=piVec[jj]*crossArray[jj,1]/countsVec[jj]
    Qmat=np.zeros([nBins,nBins])
    for iRow in np.arange(0,nEdges-1):
        Qmat[iRow,iRow+1]=NijMat[iRow,iRow+1]/Ri[iRow]
        Qmat[iRow+1,iRow]=NijMat[iRow+1,iRow]/Ri[iRow+1]
        
    Qrows=np.nonzero(np.sum(Qmat,axis=1)>0)[0]
    Qmat=Qmat[Qrows[:,None],Qrows]
    
    for iRow,row in enumerate(Qmat):
        Qmat[iRow,iRow]=-np.sum(row)
    
    bVec=np.zeros(nEdges)-1
    tauVec=sp.linalg.lstsq(Qmat,bVec)
    
    return (escapeMat,piVec,rMat,crossArray,countsVec,Ri,NijMat,Qmat,Qrows,tauVec)