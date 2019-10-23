import numpy as np
import pandas as pd

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
            for iInd,Ind in enumerate(transInds):
                escapeMat[iVal,Ind-1]=1.*transCounts[iInd]/cCount
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
    
    return (escapeMat,outVec)