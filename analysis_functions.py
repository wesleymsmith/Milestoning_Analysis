from __future__ import print_function
import numpy as np
import scipy as sp
import pandas as pd
import itertools
import copy
import gc
import sklearn as skl
from sklearn.utils import resample

def parse_multi_indexed_milestone_data(dataTable,
                                   window_index_columns,coordinate_index_columns,
                                   replica_index_columns=None,join_str='_',
                                   windowIndOneDname='wi_1D',
                                   coordinateIndOneDname='ci_1D',
                                   repIndName='Rep',
                                   verbose=False):
    #We first check to make sure that there are the same number of coordinate columns
    #as data columns. As a word of warning, the window and coordinate index columns
    #must also occur in the same order! I.e. the first column listed in the window
    #index columns should correspond to the first column listed in the coordinate index
    #columns.
    if verbose:
        print('Validating input data structure and extracting specified columns')
    if len(window_index_columns) != len(coordinate_index_columns):
        raise ValueError(
            '''ERROR: the number of window index columns {lw} 
               does not match the number of coordinate index columns {cw}'''.format(
            lw=len(window_index_columns),cw=len(coordinate_index_columns)))
    
    datCols=np.concatenate([window_index_columns,coordinate_index_columns])
    
    #Next we extract the needed data subset from the data table provided, a preliminary
    #check is made to ensure that the required columns are present.
    #In the case that no replica index columns are indicated, they will be added to the
    #data table under the assumption that all data is from a single replica
    if replica_index_columns is None:
        for colName in datCols:
            if (not (colName in dataTable.columns)):
                raise ValueError('''ERROR: The column {cn} is not present in the data!'''.format(
                        colName))
        tempData=dataTable[datCols]
        tempData[repIndName]=1
        replica_index_columns=[repIndName]
        dataCols=simData.columns
                                
    else:
        datCols=np.concatenate([datCols,replica_index_columns])
        for colName in datCols:
            if (not (colName in dataTable.columns)):
                raise ValueError('''ERROR: The column {cn} is not present in the data!'''.format(
                        colName))
        tempData=dataTable[datCols]
        if verbose:
            print('Creating one dimensional replica index column')
        tempData[repIndName]=tempData[replica_index_columns].apply(
            lambda x: join_str.join(map(str,x)),axis=1)
    
    if verbose:
        print('Creating stringified window index column')
    tempData[windowIndOneDname+'_str']=tempData[window_index_columns].apply(
        lambda x: join_str.join(map(str,x)),axis=1)
    if verbose:
        print('Creating stringified coordinate index column')
    tempData[coordinateIndOneDname+'_str']=tempData[coordinate_index_columns].apply(
        lambda x: join_str.join(map(str,x)),axis=1)
    if verbose:
        print('pruning data table')
    tempData=tempData[[windowIndOneDname+'_str',repIndName,coordinateIndOneDname+'_str']]
    if verbose:
        print(tempData.head())
            
    if verbose:
        print ('Generating indexing maps:')
    #Next, window and coordinate index columns will be converted to a single index by concatenating
    #their stringified values from each row and creating and taking the sorted array of unique
    #single indices as a mapping from the one dimensional index (single index) to the original
    #N-D index (where N is the number of window / coordinate columns)
    #This array is then used to generate the N-D to one-D mapping in dictionary form where
    #Doing this requires a joining string (sepecified using the 'join_str' argument). The default
    #is an underscore... make sure this join string is not found in the entries of any coordinate
    #or window indexing columns.
    
    #First we collect a list all unique observed coordinate and window indices for each column
    if verbose:
        print('-one dimensional index to N dimension indices map')
    OneD_to_ND_map=tempData[windowIndOneDname+'_str'].sort_values().unique()
    OneD_to_ND_map=np.unique(np.concatenate([
        OneD_to_ND_map,
        tempData[coordinateIndOneDname+'_str'].sort_values().unique()]))
    
    if verbose:
        print('-N dimension indices to one dimensional index map')
    ND_to_OneD_map={
        entry:iEntry for iEntry,entry in enumerate(OneD_to_ND_map)
    }
    
    if verbose:
        print("mapping stringified window index to one dimensional integer index")
    tempData[windowIndOneDname]=tempData[windowIndOneDname+'_str'].map(
        lambda x: ND_to_OneD_map[x])
    
    if verbose:
        print("mapping stringified coordinate index to one dimensional integer index")
    tempData[coordinateIndOneDname]=tempData[coordinateIndOneDname+'_str'].map(
        lambda x: ND_to_OneD_map[x])
    
    return {'OneD_to_ND_map':OneD_to_ND_map,
            'ND_to_OneD_map':ND_to_OneD_map,
            'parsed_data_frame':tempData,
           }

def compute_bin_vector(binInd,xIndSeries,
                       binSet=None,giveBins=False,
                       giveDeltaVal=False):
    if binSet is None:
        #transfrom binInd and xIndSeries such that the minimum observed value is 1
        seriesMin=np.min(xIndSeries)
        seriesMax=np.max(xIndSeries)
        minVal=np.min([binInd,seriesMin])
        maxVal=np.max([binInd,seriesMax])
        bins=np.arange(minVal,maxVal+1)
    else:
        minVal=np.min(binSet)
        bins=binSet
    
    deltaVal=1-minVal
    
    binVal=binInd+deltaVal
    xVals=np.array(xIndSeries+deltaVal)
    
    #binC <- boolean array, true when xVals matches binVal
    binC=np.array(xVals==binVal)
    
    #binT <- contains the value from xVals when xVals has just escaped from
    #        the central bin or 0 otherwise
    binT=(1-binC[1:])*binC[:-1]*xVals[1:]
    
    binVec=np.array(binC,dtype=int)*binVal
    binVec[1:]=binT+binVal*binC[1:]
    
    if giveBins | giveDeltaVal:
        outList=[binVec]
        outData={}
        if giveBins:
            outData['bins']=bins
        if giveDeltaVal:
            outData['deltaVal']=deltaVal
        outList.append(outData)
        return(outList)
    else:
        return(binVec)

def compute_reentry_vector(binInd,xIndSeries,
                           binSet=None,giveBins=False,
                           giveDeltaVal=False,verbose=False):
    if binSet is None:
        #transfrom binInd and xIndSeries such that the minimum observed value is 1
        seriesMin=np.min(xIndSeries)
        seriesMax=np.max(xIndSeries)
        minVal=np.min([binInd,seriesMin])
        maxVal=np.max([binInd,seriesMax])
        bins=np.arange(minVal,maxVal+1)
    else:
        minVal=np.min(binSet)
        bins=binSet
    
    deltaVal=1-minVal
    
    binVal=binInd+deltaVal
    xVals=np.array(xIndSeries+deltaVal)
    
    binC=np.array(xVals==binVal)
    
    reentries=(1-binC[:-1])*binC[1:]*xVals[:-1]
    binR=np.array(binC,dtype=int)*binVal
    binR[:-1]=binR[:-1]+reentries
    if verbose:
        print('binR:',binR)
    
    binRuns=[list(g) for k,g in itertools.groupby(binR)]

    runPairs=[[run[0],len(run)] for run in binRuns]
    if verbose:
        print('runPairs:',runPairs)
    
    lastEscape=0
    reentriesList=[]
    for runPair in runPairs:
        if runPair[0]==0:
            reentriesList.append([0]*runPair[1])
        elif runPair[0]==binVal:
            reentriesList.append([lastEscape]*runPair[1])
        else:
            lastEscape=runPair[0]
            reentriesList.append([runPair[0]]*runPair[1])
    reentryVec=np.concatenate(reentriesList)
    
    if giveBins | giveDeltaVal:
        outList=[reentryVec]
        outData={}
        if giveBins:
            outData['bins']=bins
        if giveDeltaVal:
            outData['deltaVal']=deltaVal
        outList.append(outData)
        return(outList)
    else:
        return(reentryVec)
    
def compute_bin_escape_counts(binInd,xIndSeries,
                              binSet=None,
                              giveBins=False,
                              giveDeltaVal=False,
                              giveBinVec=False,
                              giveBinT=False):
    '''
        Takes as input a given bin index and a time series of x indices
        (can be any iterable array).
        Returns a count of how many times the x index was equal to the bin index
        given along with each bin into which the x index was observed to escape to
        and the number of times that escape was seen. This is returned as dictionary object
        outputDict={
            count:##totalCount##  <-- number of times x was in the given bin plus total numer of escape events
            escapes:(observedEscapeBins,numberOfEscapesToObservedEscapeBin)
        }
        Note: binInd and xIndSeries should be integer valued. During processing they are transformed such that the
        minimum observed value (over binInd and all xIndSeries data) is unity. This simplifies the algorithm somewhat.
        However, the returned 'escapes' entry in the output will match the input series.
    '''
    
    binVec,binningInfo=compute_bin_vector(
        binInd=binInd,xIndSeries=xIndSeries,binSet=binSet,giveBins=True,giveDeltaVal=True)
    deltaVal=binningInfo['deltaVal']
    bins=binningInfo['bins']
    
    #tCounts <- the observed escape bins along with the number of times each of those bins was escaped into
    centerBin=binInd+deltaVal
    binT=(binVec*(1-(binVec==centerBin)))[1:]
    tCounts=np.unique(binT,return_counts=True) #tuple (observedValues,numberOfTimesObserved)
    
    outDict={'count':np.sum(binVec>0),
            'escapes':(tCounts[0][1:]-deltaVal,tCounts[1][1:])}
    if giveBins:
        outDict['binSet']=bins
    if giveBinVec:
        outDict['binVec']=binVec
    if giveBinT:
        outDict['binT']=binT
    if giveDeltaVal:
        outDict['deltaVal']=deltaVal
    
    return(outDict)

def compute_reentry_counts(binInd,xIndSeries,
                              binSet=None,
                              giveBins=False,
                              giveDeltaVal=False,
                              giveReentryVec=False):
    reentryVec,binningInfo=compute_reentry_vector(
        binInd=binInd,xIndSeries=xIndSeries,binSet=binSet,giveBins=True,giveDeltaVal=True)
    deltaVal=binningInfo['deltaVal']
    bins=binningInfo['bins']
    
    rCounts=np.unique(reentryVec,return_counts=True)
    
    outDict={'count':np.sum(reentryVec>0),
             'reentries':(rCounts[0][1:]-deltaVal,rCounts[1][1:])}
    
    if giveBins:
        outDict['binSet']=bins
    if giveReentryVec:
        outDict['reentryVec']=reentryVec
    if giveDeltaVal:
        outDict['deltaVal']=deltaVal
        
    return(outDict)

def build_bin_mappings(binSet):
    binSortArr=np.argsort(binSet)
    bins=binSet[binSortArr]
    binMap={}
    binSetMap={}
    for iBin,binName in enumerate(bins):
        binMap[binName]=iBin
        binSetMap[binName]=binSortArr[iBin]
    return({'bins':bins,
            'binSortArr':binSortArr,
            'binMap':binMap,
            'binSetMap':binSetMap})

def build_edge_mappings(nBins):
    '''maps between pairs of bin indices and edge indices
         the first bin index in the pair must be the
         smallest of the two. Wrapper function are used 
         to enforce this convention'''
    ePairToIndMap=np.zeros((nBins,nBins),dtype=int)-1
    for ii in np.arange(nBins-1):
        for jj in np.arange(ii+1,nBins):
            ePairToIndMap[ii,jj]=int((ii)*((nBins)+(nBins-ii-1))/2)+jj-ii-1
    tempInds=np.nonzero(ePairToIndMap>-1)
    eIndToPairMap=np.array([
        [iPair[0],iPair[1]] for iPair in \
        zip(tempInds[0],tempInds[1])
    ])
    nEdges=len(tempInds[0])
    return({'edgeIndToPair':eIndToPairMap,
            'edgePairToInd':ePairToIndMap,
            'nEdges':nEdges})

def compute_bin_edge_transitions(binInd,escapeVec,reentryVec,binSet,
                             edgeMaps=None,giveBins=False,giveBinMap=False,
                             giveEdgeMaps=False,verbose=False):
    '''
        Returns Nij_alpha and Ri_alpha based on the escape vector and reentry vector given.
        These are assumed to be given in terms of the bins contained in "binSet"
        and in areas where the index in not either escaping, reentering or inside
        the main bin (binInd) the value will be one less than the minimum value
        inside "binSet".
        Internally, binSet gets mapped to bins indexed from 0 to nBins
        Note that binInd, escapeVec, and reentryVec should be in terms of the bins in "binSet",
        Nij_alpha and Ri_alpha are in terms of the edges between these bins.
        To remain general, a bin may have more than multiple bins adjacent to it
        and thus could have many edges. We thus will create mappings to map between
        bin edges denoted as pairs of bins and bin edges with a single indexing value.
        The latter notation enables storing Nij in a 2D matrix format. This allows us
        to build this matrix by summing Nij matrices from each window.
        Edges between bins can be donoted either using the pair of bins which the edge
        divides or as an absolute 1D index. Since these edges are non-directional
        we here adopt the convention that edge pairs will always list the bin with the
        lowest index as the first value of the pair (see "build_edge_mappings" function).
    '''
    binMappingDict=build_bin_mappings(binSet)
    binSortArr=binMappingDict['binSortArr']
    bins=binMappingDict['bins']
    binMap=binMappingDict['binMap']
    binSetMap=binMappingDict['binSetMap']
    
    nBins=len(bins)
    deltaVal=1-np.min(bins)
    
    binVal=binMap[binInd]
    
    if (edgeMaps is None) | \
       (not ('edgeIndToPair' in edgeMaps)) | \
       (not ('edgePairToInd' in edgeMaps)):
        edgeMappingInfo=build_edge_mappings(nBins)
    else:
        edgeMappingInfo=edgeMaps
        
    eIndToPairMap=edgeMappingInfo['edgeIndToPair']
    ePairToIndMap=edgeMappingInfo['edgePairToInd']
    nEdges=len(eIndToPairMap)
    ePairFun=lambda ii,jj: ePairToIndMap[np.min([ii,jj]),np.max([ii,jj])]
    eIndFun=lambda eInd: eIndToPairMap[eInd] if ((eInd>=0) & (eInd<len(eIndToPairMap))) else [-1,-1]
    
    Nij_mat=sp.sparse.lil_matrix((nEdges,nEdges))
    Ri_mat=sp.sparse.lil_matrix((nBins,nBins))
    
    Rtemp=np.unique(reentryVec,return_counts=True)
    for Rentry in zip(Rtemp[0],Rtemp[1]):
        if Rentry[0] in binMap:
            ri=binMap[Rentry[0]]
            Ri_mat[binVal,ri]+=Rentry[1]
    
    center_count=np.sum(escapeVec==binInd)
    
    #need to get how many times bin i was escaped to after
    #last being in bin j
    #essentially this amounts to comparing entries from
    #escapeVec from 1 and greater with entries from
    #
    eMask=np.array(
        (1-(escapeVec==binInd))*np.isin(escapeVec,binSet),
        dtype=bool)
    eVec=escapeVec[eMask]
    rVec=reentryVec[eMask]
    rMask=np.array(
        (1-(rVec==binInd))*np.isin(rVec,binSet),
        dtype=bool)
    eVec=eVec[rMask]
    rVec=rVec[rMask]
    pMask=np.array(np.abs(eVec-rVec)>0,
                   dtype=bool)
    eVec=eVec[pMask]
    rVec=rVec[pMask]
    
    if verbose:
        print('num transition entries:',np.sum(pMask))
    for tPair in zip(rVec,eVec):
        if (tPair[0] in binMap) & (tPair[1] in binMap):
            ei=ePairFun(binVal,binMap[tPair[0]])
            ej=ePairFun(binVal,binMap[tPair[1]])
            if np.isfinite(ei) & np.isfinite(ej) & \
                (ei > -1) & (ej > -1) & \
                (not (ei==ej)):
                Nij_mat[ei,ej]+=1
    
    outDict={
        'Ri_counts':Ri_mat,
        'Nij_counts':Nij_mat,
        'center_count':center_count
    }
    
    #giveBins=False,giveBinMap=False,giveEdgeMaps=False
    if giveBins:
        outDict['bins']=bins
    if giveBinMap:
        outDict['binMap']=binMap
        outDict['binSetMap']=binSetMap
    if giveEdgeMaps:
        outDict['edgeIndToPair']=eIndToPairMap
        outDict['edgePairToInd']=ePairToIndMap
        
    return(outDict)
    
def extract_analysis_columns(milestoneData,windowColumn='Window',xIndexColumn='X_Index',frameCol='Frame',
                                      repColumn=None,groupingColumn=None):
    extractionCols=[windowColumn,xIndexColumn]
    simData=milestoneData[extractionCols]
    
    if groupingColumn is None:
        groupingCol='Group'
        simData[groupingCol]=0
    else:
        groupingCol=groupingColumn
        simData[groupingCol]=milestoneData[groupingCol]
        
    if repColumn is None:
        repCol='Rep'
        simData[repCol]=0
    else:
        repCol=repColumn
        simData[repCol]=milestoneData[repCol]
    
    if frameCol in milestoneData:
        simData[frameCol]=milestoneData[frameCol]
    else:
        simData[frameCol]=np.arange(len(simData))
    
    return(simData)

def add_indexed_milestoning_analysis_columns(milestoneData,
                                            windowColumn='Window',xIndexColumn='X_Index',frameCol='Frame',
                                            repColumn=None,groupingColumn=None,verbose=False,
                                            verboseLevel=0):
    simData=extract_analysis_columns(milestoneData,windowColumn,xIndexColumn,frameCol,
                                      repColumn,groupingColumn)
    if groupingColumn is None:
        groupingCol='Group'
    else:
        groupingCol=groupingColumn
    if repColumn is None:
        repCol='Rep'
    else:
        repCol=repColumn
    
        
    dataFrameList=[]
    groupingGroups=simData.groupby(groupingCol)
    for groupingGroup in groupingGroups:
        groupingName=groupingGroup[0]
        groupingData=groupingGroup[1]
        if verbose:
            print('--- --- --- Grouping Name:',groupingName,'--- --- ---')
        windowGroups=groupingData.groupby(windowColumn)
        binSet=np.sort(np.unique(np.concatenate([
            groupingData[windowColumn].unique(),
            groupingData[xIndexColumn].unique()
        ])))
        nBins=len(binSet)
        deltaVal=1-np.min(binSet)
        if verbose & (verboseLevel>0):
            print('\tbinSet:',binSet,'; deltaVal:',deltaVal)
        for windowGroup in windowGroups:
            windowName=windowGroup[0]
            windowData=windowGroup[1]
            if verbose:
                print('\t--- --- Window Name:',windowName,'--- ---')
                print('\t\t--- Replica Name:',end=" ")
            repGroups=windowData.groupby(repCol)
            for repGroup in repGroups:
                repName=repGroup[0]
                repData=repGroup[1]
                if verbose:
                    print(repName,end=" ")
                binVec=compute_bin_vector(
                    binInd=windowName,xIndSeries=repData[xIndexColumn],
                    binSet=binSet,giveBins=False,giveDeltaVal=False)
                repData['Escape_Vector']=binVec-deltaVal
                reentryVec=compute_reentry_vector(
                    binInd=windowName,xIndSeries=repData[xIndexColumn],
                    binSet=binSet,giveBins=False,
                    giveDeltaVal=False,verbose=False)
                reentryVec[1:]=reentryVec[:-1] #shift forward to align to escape vector
                reentryVec[0]=0
                repData['Reentry_Vector']=reentryVec-deltaVal
                dataFrameList.append(repData.copy())
                gc.collect()
            if verbose:
                print("---")
                print('\t--- --- ------ --- ---')
        if verbose:
            print('--- --- --- ------ --- --- ---')
    return(pd.concat(dataFrameList))

def bootstrap_analysis_group_pi_vector(groupDataFrame,windowColumn,binSet,
                                       bootSampleSize,nBootSamples,repColumn=None,
                                         giveBins=False,giveBinMap=False,
                                         giveEscapeMat=False,giveCounts=False,
                                         giveCountsMat=False):
    bootResults=[]
    stratifyColumns=[windowColumn]
    if (not (repColumn is None)) & \
        (repColumn in groupDataFrame):
        stratifyColumns.append(repColumn)
    print("Running sample (out of {:g}):".format(nBootSamples),end=" ")
    iSample=0
    lastSample=-1
    while iSample < nBootSamples:
        if iSample>lastSample:
            print(iSample+1,end="")
        else:
            print('',end="")
        lastSample=iSample
        bootData=resample(groupDataFrame,
         n_samples=bootSampleSize,replace=True,
         stratify=groupDataFrame[stratifyColumns])
        piResults=compute_analysis_group_pi_vector(
            bootData,windowColumn,binSet,
            giveBins=False,giveBinMap=False,
            giveEscapeMat=True,giveCounts=False,
            giveCountsMat=False)
        eMat=piResults['escapeMat']
        netBreaks=np.sum(np.abs((eMat>0)-(eMat.T>0)))
        if netBreaks==0:
            iSample+=1
            bootResults.append(piResults['piVec'])
            print(",",end=" ")
        else:
            print('x{:g}'.format(netBreaks),end="")
    print(" ")
    return(bootResults)
        
def compute_analysis_group_pi_vector(groupDataFrame,windowColumn,binSet,
                                     giveBins=False,giveBinMap=False,
                                     giveEscapeMat=False,giveCounts=False,
                                     giveCountsMat=False):
    windowGroups=groupDataFrame.groupby(windowColumn)
    
    binMappingDict=build_bin_mappings(binSet)
    binSortArr=binMappingDict['binSortArr']
    bins=binMappingDict['bins']
    binMap=binMappingDict['binMap']
    binSetMap=binMappingDict['binSetMap']
    
    nBins=len(bins)
    deltaVal=1-np.min(bins)

    nBins=len(bins)
    revBinMap=binSetMap
    
    escapeMat=sp.sparse.lil_matrix((nBins,nBins),dtype=float)
    countArray=np.zeros(nBins,dtype=int)
    
    if giveCountsMat:
        countsMat=sp.sparse.lil_matrix((nBins,nBins),dtype=float)
    
    for windowGroup in windowGroups:
        windowName=windowGroup[0]
        iWin=binMap[windowName]
        windowData=windowGroup[1]
        eVec=np.array(windowData['Escape_Vector'],dtype=int)
        eCounts=np.unique(eVec,return_counts=True)
        for eBinName,eBinCount in zip(eCounts[0],eCounts[1]):
            if eBinName in binMap:
                iBin=binMap[eBinName]
                if not(iBin==iWin):
                    escapeMat[iWin,iBin]+=eBinCount
                else:    
                    countArray[iWin]+=eBinCount
                if giveCountsMat:
                    countsMat[iWin,iBin]+=eBinCount
        escapeMat[iWin,:]=escapeMat[iWin,:]/countArray[iWin]
        escapeMat[iWin,iWin]=1-np.sum(escapeMat[iWin,:])
    escapeEig=np.linalg.eig(escapeMat.todense().T)
    si=np.argsort(1-escapeEig[0])
    piVec=np.array(escapeEig[1])[:,si[0]]
    piVec=piVec/np.sum(piVec)
    
    if giveEscapeMat | giveCounts | giveBins | giveBinMap:
        outDataDict={'piVec':copy.deepcopy(piVec)}
        if giveEscapeMat:
            outDataDict['escapeMat']=copy.deepcopy(escapeMat)
        if giveCounts:
            outDataDict['counts']=copy.deepcopy(countArray)
        if giveBins:
            outDataDict['bins']=copy.deepcopy(bins)
        if giveBinMap:
            outDataDict['binMap']=copy.deepcopy(binMap)
            outDataDict['binSetMap']=copy.deepcopy(binSetMap)
        if giveCountsMat:
            outDataDict['countsMat']=copy.deepcopy(countsMat)
        return(outDataDict)
    else:
        return(copy.deepcopy(piVec))
    
def filter_matRow_entries(mat,giveEntryMap=True,
                             matRowAggTestFun=lambda xMat:np.sum(np.abs(xMat),axis=1)>0,
                             rowIndTestFun=lambda iRow: True):
    targets=np.array([
        iRow for iRow in np.nonzero(matRowAggTestFun(mat))[0] \
        if rowIndTestFun(iRow)
    ])
    
    matRed=mat[targets[:,None],targets]
    if giveEntryMap:
        return((matRed,targets))
    else:
        return(matRed)
    
def get_tau(qhat,entryMap,indToPairMap,
            bVec=None,sinkInds=None,
            verbose=False):
    if bVec is None:
        bvec=np.zeros(len(qhat))-1
    elif len(bVec) != len(qhat):
        bvec=np.zeros(len(qhat))-1
    else:
        bvec=copy.deepcopy(bVec)
    
    Qhat=copy.deepcopy(qhat)
    if not (sinkInds is None):
        for sinkInd in sinkInds:
            Qhat[sinkInd,:]=0
            Qhat[:,sinkInd]=0
            #bvec[sinkInd]=0
        if verbose:
            print('new Qhat:',Qhat)
            #print('new bvec',bvec)
    else:
        #Qhat[-1,:]=0
        #Qhat[:,-1]=0
        Qhat=Qhat[:-1,:-1]
        bvec=np.zeros(len(Qhat))-1
        #bvec[-1]=0
        if verbose:
            print('assuming last window is sink')
            print('new Qhat:',Qhat)
            print('new bvec',bvec)
    tauData=np.linalg.lstsq(Qhat,bvec,rcond=None)
    tauVec=tauData[0]
    
    pairStr='{:g}_{:g}'
    tauDict={}
    for iTau,tau in enumerate(tauVec):
        binPair=indToPairMap[entryMap[iTau]]
        if verbose:
            print(binPair)
        binPairStr=pairStr.format(binPair[0],binPair[1])
        if binPairStr in tauDict:
            tauDict[binPairStr]=np.min([tauDict[binPairStr],tau])
        else:
            tauDict[binPairStr]=tau
    return((tauDict,tauVec,tauData))
    
def compute_analysis_group_Qdata(groupDataFrame,windowColumn,binSet,
                                 sourceWindows=None,sinkWindows=None,
                                     repColumn='Rep',givePiVec=True,
                                     giveBins=False,giveBinMap=False,
                                     giveEscapeMat=False,giveCounts=False,
                                     giveCountsMat=False,giveEdgeMap=False,
                                     giveQtargets=True):
    
    piDataDict=compute_analysis_group_pi_vector(groupDataFrame,windowColumn,binSet,
                                     giveBins=True,giveBinMap=True,
                                     giveEscapeMat=True,giveCounts=True,
                                     giveCountsMat=True)
    bins=piDataDict['bins']
    binMap=piDataDict['binMap']
    binSetMap=piDataDict['binSetMap']
    countsVec=piDataDict['counts']
    escapeMat=piDataDict['countsMat']
    
    nBins=len(bins)
    
    edgeInfo=build_edge_mappings(nBins)
    nEdges=len(edgeInfo['edgeIndToPair'])
    #print(edgeInfo)
    #print(nEdges)
    
    eIndToPair=edgeInfo['edgeIndToPair']
    ePairToInd=edgeInfo['edgePairToInd']
    
    windowGroups=groupDataFrame.groupby(windowColumn)
    
    
        
    piData=compute_analysis_group_pi_vector(
        groupDataFrame,windowColumn,binSet,
        giveBins=False,giveBinMap=False,
        giveEscapeMat=giveEscapeMat,giveCounts=True,
        giveCountsMat=True)
    
    piVec=piData['piVec']
    countVec=piData['counts']
    Tcount=np.sum(piVec/countVec)**-1
   
    
    Rmat=sp.sparse.lil_matrix((nBins,nBins))
    Nmat=sp.sparse.lil_matrix((nEdges,nEdges))
    print('Computing R and N:',end=" ")
    for windowGroup in windowGroups:
        repGroups=windowGroup[1].groupby(repColumn)
        print("(",windowGroup[0],':',end=" ")
        for repGroup in repGroups:
            print(repGroup[0],end=" ")
            transitionData=compute_bin_edge_transitions(
                binInd=windowGroup[0],
                escapeVec=np.array(repGroup[1]['Escape_Vector']),
                reentryVec=np.array(repGroup[1]['Reentry_Vector']),
                binSet=binSet,
                edgeMaps=edgeInfo,giveBins=False,giveBinMap=False,
                giveEdgeMaps=False,verbose=False)
            iBin=binMap[windowGroup[0]]
            #print('Nmat shape:',Nmat.shape)
            #print('Nij_counts shape:',transitionData['Nij_counts'].shape)
            Rmat=Rmat+transitionData['Ri_counts']*piVec[iBin]/countVec[iBin]
            Nmat=Nmat+transitionData['Nij_counts']*piVec[iBin]/countVec[iBin]
        print(")",end=" ")
    print("")
    Rmat=Tcount*Rmat
    Nmat=Tcount*Nmat
    outDict={
        'Rmat':Rmat,
        'Nmat':Nmat}
    
    if giveEscapeMat:
        outDict['escapeMat']=piData['escapeMat']
    if giveCountsMat:
        outDict['countsMat']=piData['countsMat']
    if givePiVec:
        outDict['piVec']=piVec
    if giveEdgeMap:
        outDict['edgeMap']=edgeInfo

    Ri_vec=np.zeros(nEdges)
    Rpairs=np.nonzero(Rmat)
    for Rpair in zip(Rpairs[0],Rpairs[1]):
        Redge=ePairToInd[np.min(Rpair),np.max(Rpair)]
        Ri_vec[Redge]+=Rmat[Rpair[0],Rpair[1]]
    Qmat=sp.sparse.lil_matrix((nEdges,nEdges))
    for Redge in np.nonzero(Ri_vec)[0]:
        Qmat[Redge,:]=Nmat[Redge,:]/Ri_vec[Redge]
    
    outDict['Qmat']=Qmat
    
    Qlap=copy.deepcopy(Qmat)
    for iRow,rowVec in enumerate(Qlap):
        Qlap[iRow,iRow]=-np.sum(rowVec)

    nNodes=len(bins)
    
    rowAggTest=lambda tMat: np.sum(np.abs(tMat),axis=1)>0
    indCheckFun=lambda iRow: True
    #indCheckFun=lambda iRow: not ((nNodes-1) in eIndToPair[iRow])
    QlapRed,Qtargets=filter_matRow_entries(
        Qlap,giveEntryMap=True,
        matRowAggTestFun=rowAggTest,rowIndTestFun=indCheckFun)
    
    QlapRed=QlapRed.todense()
    
    outDict['Qhat']=QlapRed
    
    if giveQtargets:
        outDict['Qtargets']=Qtargets

    bvec=np.zeros(len(QlapRed))-1
    
    tauDict,tauVec,tauVecData=get_tau(
        qhat=QlapRed,entryMap=Qtargets,indToPairMap=eIndToPair)
    
    outDict['tauData']=tauVecData
    outDict['tauVec']=tauVec
    outDict['tauDict']=tauDict
    
    return(outDict)
            
def analyze_indexed_milestoning_escapes(milestoneData,windowColumn='Window',xIndexColumn='X_Index',
                                      repColumn=None,groupingColumn=None,
                                      giveEscapeMats=False,giveCounts=False,
                                      giveBins=False,giveBinMaps=False,giveCountsMat=False):
    extractionCols=[windowColumn,xIndexColumn]
    simData=milestoneData[extractionCols]
    
    if groupingColumn is None:
        groupingCol='Group'
        simData[groupingCol]=0
    else:
        groupingCol=groupingColumn
        simData[groupingCol]=milestoneData[groupingCol]
        
    if repColumn is None:
        repCol='Rep'
        simData[repCol]=0
    else:
        repCol=repColumn
        simData[repCol]=milestoneData[repCol]
        
    simGroups=simData.groupby(groupingCol)
    outDataDict={}
    for simGroup in simGroups:
        groupName=simGroup[0]
        groupData=simGroup[1]
        print("Group:",groupName)
        groupDataDict={}
        binSet=np.unique(np.concatenate([
            groupData[windowColumn].unique(),
            groupData[xIndexColumn].unique()
        ]))
        if giveBins:
            groupDataDict['bins']=copy.deepcopy(binSet)
        
        nBins=len(binSet)
        binMap={}
        for iBin,binName in enumerate(binSet):
            binMap[binName]=iBin
        if giveBinMaps:
            groupDataDict['binMap']=copy.deepcopy(binMap)
        
        escapeMat=sp.sparse.lil_matrix((nBins,nBins),dtype=float)
        countArray=np.zeros(nBins,dtype=int)
        
        if giveCountsMat:
            countsMat=sp.sparse.lil_matrix((nBins,nBins),dtype=float)
        
        groupWindows=groupData.groupby(windowColumn)
        for groupWindow in groupWindows:
            windowName=groupWindow[0]
            windowData=groupWindow[1]
            iWin=binMap[windowName]
            print("\tWindow:",windowName,"(iWin = ",iWin,")")
            windowReps=windowData.groupby(repCol)
            for windowRep in windowReps:
                repName=windowRep[0]
                repData=windowRep[1]
                print("\t\tRep:",repName)
                repEscapeData=compute_bin_escape_counts(
                      windowName,repData['X_Index'],
                      binSet=binSet)
                countArray[iWin]+=repEscapeData['count']
                if giveCountsMat:
                    countsMat[iWin,iWin]+=repEscapeData['count']
                for escapeBinName,escapeCount in \
                    zip(repEscapeData['escapes'][0],
                        repEscapeData['escapes'][1]):
                    iEscapeBin=binMap[escapeBinName]
                    escapeMat[iWin,iEscapeBin]+=escapeCount
                    if giveCountsMat:
                        countsMat[iWin,iEscapeBin]+=escapeCount
                print('\t\tcount:',repEscapeData['count'],
                      'escapes:',repEscapeData['escapes'])
                print("\t\t---")
            escapeMat[iWin,:]=escapeMat[iWin,:]/countArray[iWin]
            escapeMat[iWin,iWin]=1-np.sum(escapeMat[iWin,:])
            print("\t--- ---")
            if giveEscapeMats:
                groupDataDict['escapeMatrix']=copy.deepcopy(escapeMat)
            escapeEig=np.linalg.eig(escapeMat.todense().T)
            si=np.argsort(1-escapeEig[0])
            piVec=np.array(escapeEig[1])[:,si[0]]
            piVec=piVec/np.sum(piVec)
            groupDataDict['piVector']=copy.deepcopy(piVec)
            if giveCounts:
                groupDataDict['counts']=copy.deepcopy(countArray)
            if giveCountsMat:
                groupDataDict['countsMat']=copy.deepcopy(countsMat)
        print("--- --- ---")
        #update output dictionary with current group's data
        outDataDict[groupName]=copy.deepcopy(groupDataDict)
        gc.collect()
    
    return(outDataDict)
    
def gen_edge_ordering_maps(nBins):
    tupleToIndex=np.zeros(shape=(nBins,nBins))
    
def compute_bin_crossing_counts(binInd,xIndSeries,
                                edgeOrderingMaps=None,
                                escapeData=None):
    if (escapeData is None) | \
        (not ('binC' in escapeData)):
        escapeData=compute_bin_escape_counts(
            binInd,xIndSeries,giveWorkingArrays=True)

    deltaVal=escapeData['deltaVal']
    binC=escapeData['binC']
    binT=escapeData['binT']
    escapes=escapeData['escapes']
    bins=escapeData['binSet']
    nBins=len(bins)
    
    xVals=(xIndSeries+deltaVal)
    binVal=binInd+deltaVal
    
    
    '''
        Crossings are transitions from one edge (boundary) of a bin
        to a different edge of that bin. Since the bins are potentially
        N-Dimensional objects, a given bin may be in contact with many
        other bins (potentially all other bins in the case of a very compact
        binning setup). A bin edge may be denoted by the pair bins it lies between.
        Similarly, then, a crossing is denoted by the pair of edges being traversed.
        So just as there could be up to nBins^2 edges, there could be up to nEdges^2
        crossings or equivalently, there could be up to nBins^4 possible crossings in
        an nBin setup. If we were to store this in a dense array format, the space
        requirements could very easily become intractible. Thus, we will store the
        crossings in a sparse array format.
        Likewise, we will need to generate mappings to map between tuple and index
        notation for edges. 
        As a note, these 'edges' are just boundaries between
        pairs of bins so the ordering does not matter. E.g. edge(1,2)==edge(2,1).
        Thus we will base our indexing such that the lower bin index is always the
        first value of the pair.
        On the other hand, the ordering of crossings DOES matter, so we must take
        care to preserve ordering there.
    '''
    
    if edgeOrderingMaps is None:
        edgeTupleToIndex
        

def analyze_indexed_milestone_data(
    dataTable,binCol='Window',xIndCol='X_Index',repCol=None):
    '''
    Warning: this algorithm will scale as N^4 in storage, where
        N is the total number of milestoning bins.
    '''
    
    eStr='%g_%g'
    if repCol is None:
        simData=dataTable[[binCol,xIndCol]].copy()
        simData[repCol]=0
    else:
        simData=dataTable[[binCol,xIndCol,repCol]].copy()
        
    minBin_original=np.min([simData[binCol].min(),simData[xIndCol].min()])
    maxBin_original=np.max([simData[binCol].max(),simData[xIndCol].max()])
    nBins=maxBin-minBin
    
    binDelta=1-minBin
    
    simData[binCol]=simData[binCol]-binDelta
    simData[xIndCol]=simData[xIndCol]-binDelta
    
    binExitMat=np.zeros(shape=(nBins,nBins))
    
    
    
    

def analyze_milestone1D_data(dataTable,windowMins=None,windowMaxs=None,
                             useInds=False,verbose=False,multiReplica=False):
    #if "useInds" is set to true, the 'X_Index' column must be present in the data table
    #otherwise, if either useInds is false or there is no 'X_Index' column,
    #an 'X_Index' column is generated internally.
    dataCols=["Window","Time"]
    if multiReplica:
        dataCols=np.concatenate([dataCols,["Rep"]])
    if not (useInds & ('X_Index' in dataTable)):
        if (winMins is None) | (winMaxs is None) | (not ('X' in simData.columns)):
            raise ValueError(
                "Error:  winMins and winMaxs must be defined and X column must be present when X_Index is not provided!")
        dataCols=np.concatenate([dataCols,['X']])
        #need to generate indexing data internally
        winMins=np.array(windowMins)
        winMaxs=np.array(windowMaxs)
        windowCenters=(windowMins+windowMaxs)/2.0
        simData=dataTable[dataCols]
        binEdges=winMins[1:] #np.concatenate([winMins,[winMaxs[-1]]])
        digitize_kwds={"bins":binEdges}
        simData['X_Index']=simData.X.apply(np.digitize,**digitize_kwds)
    else:
        dataCols=np.concatenate([dataCols,["X_Index"]])
        simData=dataTable[dataCols]
    
    windows=np.sort(simData.Window.unique())
    #minWindow=np.min(windows)
    #print('minWindow: ',minWindow)
    #windows=windows-minWindow
    xbins=np.sort(simData.X_Index.unique())#-minWindow
    nBins=len(xbins)
    nEdges=nBins-1
    escapeMat=np.zeros((nBins,nBins))
    rMat=np.zeros([nBins,nEdges])
    crossArray=np.zeros([nBins,2])
    tSum=0
    countsVec=np.zeros(nBins)
    
    #deltaInd=1-np.min(xbins)
    #print('deltaInd: ',deltaInd)
    #iVal->escape matrix row index
    #xbin->window
    #cVal->place holder for bin index with indexing starting at 1
    for iVal,xbin in enumerate(xbins):
        if xbin in windows:
            if verbose:
                print ("--- --- ---")
            winDat=simData[simData.Window==xbin]
            if not multiReplica:
                tempDat=winDat
                cVal=xbin+1 #deltaInd
                binVec=np.array(tempDat.X_Index+1) #+deltaInd-minWindow)
                binC=(binVec==cVal)
                binT=(1-binC[1:])*binC[:-1]*binVec[1:]
                tCounts=np.unique(binT,return_counts=True)
                transInds=tCounts[0][1:] #first entry should always be for binT==0
                transCounts=tCounts[1][1:]
                cCount=np.sum(binC)
                #generate escape matrix
                for iInd,Ind in enumerate(transInds):
                    escapeMat[iVal,Ind-1]=escapeMat[iVal,Ind-1]+1.*transCounts[iInd]/cCount
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
                countsVec[iVal]=countsVec[iVal]+np.sum(binC)
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
            else:
                if verbose:
                    print ("working on replica:",end=" ")
                for rep in winDat.Rep.unique():
                    if verbose:
                         print (rep,end=" ")
                    tempDat=winDat[winDat.Rep==rep]
                    deltaInd=1-tempDat['X_Index'].min()
                    cVal=xbin+1
                    binVec=np.array(tempDat.X_Index+1)
                    binC=(binVec==cVal)
                    binT=(1-binC[1:])*binC[:-1]*binVec[1:]
                    tCounts=np.unique(binT,return_counts=True)
                    transInds=tCounts[0][1:] 
                    transCounts=tCounts[1][1:]
                    cCount=np.sum(binC)
                    for iInd,Ind in enumerate(transInds):
                        escapeMat[iVal,Ind-1]=escapeMat[iVal,Ind-1]+1.*transCounts[iInd]/cCount
                    runVec=binT[np.nonzero(binC[1:]+binT)] 
                    runList=np.array([[int(j[0]),len(j)] for j in \
                                      [list(g) for k,g in itertools.groupby(runVec)]])
                    for iRun,run in enumerate(runList[1:]):
                        if run[0]==0:
                            rMat[iVal,np.min([runList[iRun,0]-1,xbin])]+=run[1]
                    countsVec[iVal]=countsVec[iVal]+np.sum(binC)
                    crossings=(np.array([[int(j[0]),len(j)] for j in \
                                      [list(g) for k,g in itertools.groupby(binT[np.nonzero(binT)])]]
                                      )[:,0]).flatten()
                    crossings=(crossings[1:]-crossings[:-1])/2
                    crossArray[iVal,0]=crossArray[iVal,0]+np.sum(crossings==1)
                    crossArray[iVal,1]=crossArray[iVal,1]+np.sum(crossings==-1)
                    if verbose:
                        print ("")
            if verbose:
                print ("escapeMatrix entry for window %g:"%xbin)
                print ('['+', '.join(map(lambda x: '%.5f'%x,escapeMat[iVal,:]))+']')
                print ("Number of crossings (left-to-right,right-to-left):",end=" ")
                print ("(%g,%g)"%(crossArray[iVal,0],crossArray[iVal,1]))
    if verbose:
        print ("--- --- ---")
        
    tSum=np.sum(countsVec)
    Emat=np.matrix(escapeMat)
    Dmat=np.matrix(np.diag(1-np.sum(escapeMat,axis=1)))

    Amat=Emat+Dmat

    outEig=np.linalg.eig(Amat.T)
    si=np.argsort(1-outEig[0])
    if verbose:
        print ('Eigenvalues:',end=" ")
        print (outEig[0][si])
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
    print('init: Qmat.shape: ',Qmat.shape)
    print('CountsVec: ',countsVec)
    print('EscapeMat: ',escapeMat)
    print('piVec: ',piVec)
    print('NijMat: ',NijMat)
    print('Ri: ',Ri)
    for iRow in np.arange(0,nEdges-1):
        Qmat[iRow,iRow+1]=NijMat[iRow,iRow+1]/Ri[iRow]
        Qmat[iRow+1,iRow]=NijMat[iRow+1,iRow]/Ri[iRow+1]
    print('addElements: Qmat.shape: ',Qmat.shape)
    print (Qmat)
    Qrows=np.nonzero(np.sum(Qmat,axis=1)>0)[0]
    Qmat=Qmat[Qrows[:,None],Qrows]
    print('remove empty row/col: Qmat.shape: ',Qmat.shape)
    
    for iRow,row in enumerate(Qmat):
        Qmat[iRow,iRow]=-np.sum(row)
    print('update diag.: Qmat.shape: ',Qmat.shape)
    
    bVec=np.zeros(nEdges)-1
    print('Qmat.shape: ',Qmat.shape,', bVec.shape: ',bVec.shape)
    tauVec=sp.linalg.lstsq(Qmat,bVec)
    
    return (escapeMat,piVec,rMat,crossArray,countsVec,Ri,NijMat,Qmat,Qrows,tauVec)

def compute_escape_matrix_row_convergence(windowID,windowIndexData,
                                          multiReplica=False,verbose=False,
                                          noAgg=False):
    indData=windowIndexData[windowIndexData.Window==windowID][['X_Index']]
    if multiReplica:
        indData['Rep']=windowIndexData['Rep']
    else:
        indData['Rep']='rep1'
    #print indData.head()
    outDataTables=[]
    for rep in indData.Rep.unique():
        dInd=0
        indVec=np.array(indData['X_Index'])
        windowC=windowID
        indVec=np.array(indData[indData.Rep==rep].X_Index)
        if np.min(indVec)<=0:
            dInd=1-np.min(indVec)
            windowC=windowC+dInd
            indVec=indVec+dInd

        binSet=np.unique(indVec)
        #print binSet
        
        binC=(indVec==windowC)
        binT=(1-binC[1:])*binC[:-1]*indVec[1:]

        #print windowC
        #print dInd
        for binInd in binSet:
            outTable=pd.DataFrame({
                "Frame":np.arange(1,len(binC)),
                "Rep":[rep]*len(binC[:-1]),
                "i":[(windowC-dInd)]*len(binC[:-1]),
                "di":[(binInd-windowC)]*len(binC[:-1])
            })
            if binInd==windowC:
                outTable["N"]=np.cumsum(binC[1:]+(binT>0))
            else:
                outTable["N"]=np.cumsum(binT==binInd)
            outDataTables.append(outTable.copy())
            if verbose:
                print ('window %g, rep %s, di %g, N %g'%(
                    windowC-dInd, rep, binInd-windowC,outTable.N.max()))
            
    outDataTable=pd.concat(outDataTables)
    tempTable=outDataTable
    tempTable['NetFrames']=tempTable.Frame
    if noAgg:
        return outDataTable
    else:
        if verbose:
            print (tempTable.head())
            print ('--- --- ---')
        testAggTab=tempTable.groupby(
            ['Frame','i','di']).agg(
            {'N':np.sum,'NetFrames':np.sum}).reset_index().sort_values(
            ['di','i','Frame'])
        testAggTab['Rep']='Mean'
        testAggTab=testAggTab[tempTable.columns]
        if verbose:
            print (testAggTab.head())
            print ('--- --- ---')
        tempTable=pd.concat([tempTable,testAggTab])
        tempTable=pd.pivot_table(index=['Frame','Rep','i'],
                               columns='di',values=['N','NetFrames'],
                               data=tempTable)
        tempTable.columns=tempTable.columns.map(lambda x: '_'.join([str(xv) for xv in x]))
        tempTable=tempTable.reset_index()
        if verbose:
            print (tempTable.head())
            print ('--- --- ---')
        niCols=[colName for colName in tempTable.columns if 'N_' in colName]
        for niCol in niCols:
            if niCol !='N_0':
                tempTable[niCol]=tempTable[niCol]/tempTable['N_0']
        tempTable['N_0']=tempTable['N_0']/tempTable['NetFrames_0']
        dropCols=[colName for colName in tempTable.columns if 'NetFrames' in colName]
        tempTable=tempTable.drop(dropCols,axis=1)
        tempTable=pd.melt(frame=tempTable,
                        id_vars=[colName for colName in tempTable.columns if not ('N_' in colName)],
                        var_name='di',value_name='nu')
        tempTable.di=tempTable.di.map(lambda x: int(x.replace('N_','')))
        if verbose:
            print (tempTable)
        outDataTable=tempTable

        return(outDataTable)