# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import libraries 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import pandas as pd
import numpy as np
import regex as re
import re
import os
import itertools
import time
import networkx as nx
from parsing import networkResults
from typing import List, Dict, Any
from collections import Counter
from dataclasses import dataclass

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@dataclass
class analysisResults:
    ''' 
    Data class for main function. Returns an object containing the outputs of the cluster analysis.
    '''
    label: str
    clusters: Dict[str, Any]
    clusterNES: Dict[str, Any]
    enrichment: Dict[str, Any]

    def __getitem__(self, key):
        if key == 'label':
            return self.label
        elif key == 'clusters':
           return self.clusters
        elif key == 'clusterNES':
            return self.clusterNES
        elif key == 'enrichment':
            return self.enrichment
        else:
            raise KeyError(f'Invalid key: {key}.')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class analyseGSEA:
  def __init__(self):
    print('~~~~~~ analyseGSEA ~~~~~~\n')
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def getClusterNES(
        self,
        clusters: Dict,
        metadata: pd.DataFrame):
    '''
    This function will access all of the normalised enrichment scores for the constituent genesets in each subgraph and compute an overall average for the subgraph. 
    '''
    clusterNES = {}

    for idx, cluster in clusters.items():

        moduleNES = []
        for geneset in cluster:
            nes = metadata['nes'][metadata['node'] == geneset].values[0]
            moduleNES.append(nes)
        
        average = np.mean(moduleNES)
        clusterNES[idx] = average
    
    return clusterNES

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def extractiveSummarisation(
        self, 
        clusterList: Dict,
        metadata: pd.DataFrame
        ):
     
    ## Fetch all genesets from meta data
    allGenesets = list(metadata['node'].unique())

    ## Prepare output containers for geneset tokenisation 
    clusterTokens = {}

    ## Iterate through the genesets within each subgraph and quantify tokens 
    for idx,subgraph in clusterList.items():
        ## Remove current genesets from full list to generate a background
        backgroundGenesets = [ele for ele in allGenesets if ele not in list(subgraph)]

        ## Generate single-word tokens from background genesets
        allTokens = [item for sublist in [x.split('_')[1:] for x in backgroundGenesets] for item in sublist]

        ## Append short 2-token substrings
        allTokens += [f"{x.split('_')[1:][idx]}_{x.split('_')[1:][idx+1]}" 
                    for x in backgroundGenesets 
                    for idx in range(len(x.split('_')[1:])-1)]
        
        ## Count token frequencies
        backgroundCount = Counter(allTokens)

        ## Convert counts to relative proportion of word frequency in background
        relativeBackground = dict([[k,v/len(backgroundGenesets)] for k,v in backgroundCount.items()])

        ## Calculate local frequency of words - in the current subgraph
        tokens = [item for sublist in [x.split('_')[1:] for x in subgraph] for item in sublist] # Take slice from 1st index to omit database name (KEGG, GOBP, etc.)
        
        ## Append short 2-token substrings  
        tokens += [f"{x.split('_')[1:][idx]}_{x.split('_')[1:][idx+1]}" 
                  for x in subgraph 
                  for idx in range(len(x.split('_')[1:])-1)]
        
        ## Count token frequencies
        counter = Counter(tokens)

        ## Calculate the expected frequency of word - from background ratios given the sample size  
        ## Try to divide observed frequency by the expected frequency to give a relative enrichment
        ## If word doesn't exist in the background dataset, use the lowest enrichment score
        try:
           relativeCounts = dict([[k,v/(relativeBackground[k]*len(subgraph))] for k,v in counter.items()])
        except:
           minBackground = min(relativeBackground, key=relativeBackground.get)
           relativeCounts = dict([[k,v/(relativeBackground[minBackground]*len(subgraph))] for k,v in counter.items()])
        
        clusterTokens[idx] = relativeCounts

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
    label = graphData.label
    network = graphData.dataframe
    metadata = graphData.metadata
    filters = graphData.filters

    print('Creating instance of graph network.')
    time.sleep(0.5)

    ## Generate instance of graph object
    graph = nx.from_pandas_edgelist(network, source='node1', target='node2', edge_attr='jaccard_index')

    ## Remove edges from genesets with Jaccard index lower than user-specified jaccardFilter (retrieved from networkResults object)
    noEdge = list(filter(lambda e: e[2] < filters['jaccardFilter'], (e for e in graph.edges.data('jaccard_index'))))
    noEdgePairs = list(e[:2] for e in noEdge)
    graph.remove_edges_from(noEdgePairs) 
    
    ## Extract the subgraphs from the network
    print('Evaluating subgraphs.')
    time.sleep(0.5)
    subGraphs = nx.connected_components(graph)

    print('Accessing subgraph constituents.')
    time.sleep(0.5)

    outputClusters = {}
    ## Iteratively filter subgraphs that are larger than the user-defined clustersize threshold
    for idx,subgraph in enumerate(subGraphs):
        if len(subgraph) > clustersize:
            outputClusters[idx] = subgraph
        else: 
            continue

    print('Performing extractive summarisation of geneset annotations.')
    time.sleep(0.5)
    relativeEnrichment = self.extractiveSummarisation(outputClusters,metadata)
    
    print('Calculating mean normalised enrichment scores.')
    time.sleep(0.5)
    moduleNES = self.getClusterNES(outputClusters,metadata)

    print('Analysis complete.\n')
    time.sleep(0.5)
    print('~~~~~~ </analyseGSEA> ~~~~~~')

    output = analysisResults
    output.label = label
    output.clusters = outputClusters
    output.clusterNES = moduleNES
    output.enrichment = relativeEnrichment

    return output
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~