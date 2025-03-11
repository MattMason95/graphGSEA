# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import libraries 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import pandas as pd
import numpy as np
import regex as re
import re
import os
import itertools
import networkx as nx
<<<<<<< HEAD
from parsing import networkResults
from typing import List
from collections import Counter

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class analyseGSEA:
  def __init__(self):
    print('~~~~~~ analyseGSEA ~~~~~~')
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def getClusterNES(
        self,
        clusters: List,
        metadata: pd.DataFrame):
    '''
    This function will access all of the normalised enrichment scores for the constituent genesets in each subgraph and compute an overall average for the subgraph. 
    '''
    
    clusterIndex, clusterLists = zip(*clusters)

    averageNES = []
    moduleIndex = []

    for i in range(len(clusterLists)):
        moduleIndex.append(clusterIndex[i])
    
        moduleNES = []
        for j in clusterLists[i]:
            nes = metadata['nes'][metadata['node'] == j].values[0]
            moduleNES.append(nes)
        
        average = np.mean(moduleNES)
        averageNES.append(average)

    clusterNES = list(zip(moduleIndex,averageNES))
    
    return clusterNES

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def extractiveSummarisation(
        self, 
        clusterList: List,
        metadata: pd.DataFrame
        ):
     
    ## Calculate background counts of words
    allGenesets = list(metadata['node'].unique())
    allTokens = [item for sublist in [x.split('_')[1:] for x in allGenesets] for item in sublist]

    ## Append short 2-token substrings  
    allTokens += [f"{x.split('_')[1:][idx]}_{x.split('_')[1:][idx+1]}" 
                  for x in allGenesets 
                  for idx in range(len(x.split('_')[1:])-1)]
    
    backgroundCount = Counter(allTokens)
    ## Convert counts to relative proportion of word frequency in background
    relativeBackground = dict([[k,v/len(allGenesets)] for k,v in backgroundCount.items()])

    ## Prepare output containers for geneset tokenisation 
    commonWords = []
    clusterIndex = []

    ## Iterate through the genesets within each subgraph and quantify tokens 
    for idx,subgraph in clusterList:
        ## Calculate local frequency of words 
        tokens = [item for sublist in [x.split('_')[1:] for x in subgraph] for item in sublist] # Take slice from 1st index to omit database name (KEGG, GOBP, etc.)
        ## Append short 2-token substrings  
        tokens += [f"{x.split('_')[1:][idx]}_{x.split('_')[1:][idx+1]}" 
                  for x in subgraph 
                  for idx in range(len(x.split('_')[1:])-1)]
        counter = Counter(tokens)

        ## Calculate the expected frequency of word - from background ratios given the sample size  
        ## Divide observed frequency by the expected frequency to give a relative enrichment
        relativeCounts = dict([[k,v/(relativeBackground[k]*len(subgraph))] for k,v in counter.items()])
        # highestCounts = relativeCounts.most_common()
        clusterIndex.append(idx)
        commonWords.append(relativeCounts)

    clusterTokens = list(zip(clusterIndex,commonWords))

    return clusterTokens

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def clusterInformation(
        self,
        graphData: networkResults,
        clustersize: int =  5
        ):
    '''
    This function will perform several different manipulations:
    - Filtering the subgraph network to exclude subgraphs below a defined threshold
    - Evaluate annotation terminology with tokenisation
    - Compute average normalised enrichment score from the constituent genesets.
    '''

    from collections import Counter
    network = graphData.dataframe
    metadata = graphData.metadata

    ## Generate instance of graph object
    graph = nx.from_pandas_edgelist(network, source='node1', target='node2', edge_attr='jaccard_index')
    ## Remove edges from genesets with Jaccard < 0.5
    noEdge = list(filter(lambda e: e[2] < 0.5, (e for e in graph.edges.data('jaccard_index'))))
    noEdgePairs = list(e[:2] for e in noEdge)
    graph.remove_edges_from(noEdgePairs)
    
    ## Extract the subgraphs from the network
    subGraphs = nx.connected_components(graph)

    largeClusters = []
    graphIndex = []

    ## Iteratively filter subgraphs that are larger than the user-defined clustersize threshold
    for idx,i in enumerate(subGraphs):
        if len(i) > clustersize:
            largeClusters.append(i)
            graphIndex.append(idx)
        else:
            continue

    outputClusters = list(zip(graphIndex,largeClusters))

    relativeEnrichment = self.extractiveSummarisation(outputClusters,metadata)

    moduleNES = self.getClusterNES(outputClusters,metadata)

    return relativeEnrichment, outputClusters, moduleNES
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  
=======

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def Jaccard(moduleNames,geneLists):
  '''
  This function will be nested within the fileAccessor function to provide the capability to calculate Jaccard indices for comparing genesets. 
  The basic principle of the Jaccard index is a calculation of the percentage of overlapping elements within two lists.
  The biological assumption is that genesets with overlapping membership will be related to the same biological process, and will thus reduce redundancy/ 
  '''
  ## Build master and output containers
  masterList = list(zip(moduleNames,geneLists))
  
  node1List = []
  node2List = []
  jaccard = []
  

  ## Use itertools.combinations to iteratively compare all possible pairs in the masterlist 
  for x,y in itertools.combinations(masterList, 2):
    ## Fetch names of compared genesets (node1 and node2) and append to output containers  
    node1 = str(x[0])
    node2 = str(y[0])
    node1List.append(node1)
    node2List.append(node2)

    ## Fetch genes associated with compared genesets and assess these for overlap 
    s1 = set(x[1])
    s2 = set(y[1])
    ## set.intersetion() = Overlap; set.union() = Total unique 
    statistic = float(len(s1.intersection(s2)) / len(s1.union(s2)))
    jaccard.append(statistic)

  ## Zip together output containers for export
  jaccardOutput = list(zip(node1List,node2List,jaccard))  
  
  return jaccardOutput

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def significantGenesets(zippedList,threshold,unzip=False):
  '''
  This is a simple function to apply a statistical threshold to the input data and, when requested (unzip=True) to parse output the elements of the zipped input list.
  '''
  ## Valid conditions for unzip
  unzips = [True, False]
  if unzip not in unzips:
    raise ValueError('Unzip requires a boolean operator: True or False.')
        
  ## This loop goes through the numerical data for the output Jaccard list and applies a threshold filter - retain only those <= the statistical threshold
  filteredGenesets = [x for x in zippedList if float(x[2]) <= threshold]

  ## In most cases, I will want to have the statistical data unzipped - this function will return the unzipped elements from the filteredGeneset
  if unzip == True: 
      def unzip(iterable):
          return zip(*iterable)
      geneList, nesList, fdrList = unzip(filteredGenesets)
      return geneList, nesList, fdrList

## In the event that I want to keep the elements zipped, unzip=False will return the native filteredGeneset object
  else:
      return filteredGenesets
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def clusterInformation(graph,clustersize):
  '''
  This function will perform several different manipulations:
  - Filtering the subgraph network to exclude subgraphs below a defined threshold
  - Evaluate annotation terminology with tokenisation
  - Compute average normalised enrichment score from the constituent genesets.
  '''

  from collections import Counter

  ## Access all of the subgraphs detected within the NetworkX graph object
  subGraphs = nx.connected_components(graph)

  largeClusters = []
  graphIndex = []

  ## Iteratively filter subgraphs that are larger than the user-defined clustersize threshold
  for idx,i in enumerate(subGraphs):
      if len(i) > clustersize:
          largeClusters.append(i)
          graphIndex.append(idx)
      else:
          continue
  outputClusters = list(zip(graphIndex,largeClusters))

  ## Prepare output containers for geneset tokenisation 
  commonWords = []
  clusterIndex = []

  ## Iterate through the genesets within each subgraph and quantify tokens 
  for idx, i in outputClusters:
      tokens = [item for sublist in [x.split('_') for x in i] for item in sublist]
      counter = Counter(tokens)
      highestCounts = counter.most_common()
      clusterIndex.append(idx)
      commonWords.append(highestCounts)

  clusterTokens = list(zip(clusterIndex,commonWords))

  return clusterTokens, outputClusters

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def getClusterNES(data):
  '''
  This function will access all of the normalised enrichment scores for the constituent genesets in each subgraph and compute an overall average for the subgraph. 
  '''
  clusterIndex, clusterLists = zip(*data)

  averageNES = []
  moduleIndex = []

## Iterate through the genesets in each subgraph, extract the normalised enrichment scores, and compute an average for the whole subgraph
  for i in range(len(clusterLists)):
    moduleIndex.append(clusterIndex[i])

    moduleNES = []
    for j in clusterLists[i]:
        nes = metaDf['nes'][metaDf['node'] == j].values[0]
        moduleNES.append(nes)

    average = np.mean(moduleNES)
    averageNES.append(average)

  clusterNES = list(zip(moduleIndex,averageNES))
  return clusterNES
>>>>>>> 52c229c5b0a12b8af863c15cadf577b5d0eea6f9
