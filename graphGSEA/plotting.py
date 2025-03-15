# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import libraries 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib as mp
import matplotlib.pyplot as plt 
import palettable as pal 
import itertools
import networkx as nx
import regex as re
from parsing import networkResults
from typing import Any

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class plotGSEA:
  def __init__(self):
    print('~~~~~~ plotGSEA ~~~~~~')

  def iterativePlotting(
      data: networkResults,
      saveName: str = None,
      save=False
      ):
    '''
    This function will iteratively plot the Network graphs for each condition within the data extracted from the previous function using the Python library NetworkX.
    '''
    ## Import specific libraries
    import sklearn  
    from sklearn.preprocessing import minmax_scale

    if save and not saveName:
      raise ValueError('Please provide saveName.')

    ## Option for automatic saving of plots to output path
    saves = [True,False]
    if save not in saves:
      raise ValueError('Save requires a boolean operator: True or False.')

    network = data.dataframe
    meta = data.metadata
    thresholds = data.filters

    ## Iterate through the different conditions within the network analysis dataframe 
    for condition in network['condition'].unique():
      ## Truncate to condition 
      shortNetwork = network[network['condition'] == condition]
      shortMeta = meta[meta['condition'] == condition]

      ## Generate edgelist from dataframe
      graph = nx.from_pandas_edgelist(shortNetwork, source='node1', target='node2', edge_attr='jaccard_index')

      ## Specify that no edges are created between nodes with a Jaccard index of < 0.5 - remove tham from the edge list 
      noEdge = list(filter(lambda e: e[2] < thresholds['jaccardFilter'], (e for e in graph.edges.data('jaccard_index'))))
      noEdgePairs = list(e[:2] for e in noEdge)
      
      graph.remove_edges_from(noEdgePairs)

      ## Generate colour maps from the meta data provided - colour red if upregulated or blue if downregulated
      colorMap = []
      for node in graph:
        if node in list(shortMeta['node'][shortMeta['upregulated?'] == 1]):
          colorMap.append('#A43D40')
        else: 
          colorMap.append('#507D96')

      ## Extract the NES for each of the plotted nodes in the graph 
      nodeMap = []
      for node in graph:
        if node in list(shortMeta['node']):
          nes = shortMeta['nes'][shortMeta['node'] == node].values[0]
          nodeMap.append(nes)
        else: 
          continue     

      ## Normalise the scaling of the extracted NES values to use for node sizes
      normalisedNodeMap = minmax_scale(nodeMap,feature_range=(-2,2))

      ## Generate force-directed network graph using the NetworkX library
      pos = nx.spring_layout(graph, k=2/np.sqrt(len(graph.nodes())), iterations=90, seed=101, weight='jaccard_index',center=(0,0))
      
      ## Plotting
      fig, ax = plt.subplots(1,1,figsize=(12,12))
      fig = nx.draw(graph,pos=pos,node_color=colorMap,edgecolors='0.2',edge_color='0.2',linewidths=2,node_size=((normalisedNodeMap+3)*35),width=2)
      plt.tight_layout()

      ## Save if requested
      if save == True:
        plt.savefig(f'Images/{saveName}.png',dpi=150)
      else:
        continue
