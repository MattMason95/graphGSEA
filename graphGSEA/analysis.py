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