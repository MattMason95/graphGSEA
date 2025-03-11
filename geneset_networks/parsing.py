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
from pathlib import Path 
from typing import Union, Any, Optional, Dict, List
import time
from tqdm import tqdm
from dataclasses import dataclass

@dataclass
class networkResults:
    ''' 
    Data class for main function. Returns an object containing the network data and additional metadata.
    '''
    label: str
    dataframe: pd.DataFrame
    metadata: Dict[str, Any]
    statistics: Dict[str, Any]
    lookup: Dict[str, Any]

    def __getitem__(self, key):
        if key == 'label':
            return self.label
        elif key == 'dataframe':
           return self.dataframe
        elif key == 'metadata':
            return self.metadata
        elif key == 'statistics':
            return self.statistics
        elif key == 'lookup':
            return self.lookup
        else:
            raise KeyError(f'Invalid key: {key}.')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class parseGSEA:
  '''
  Primary parsing function for GSEA data.
  
  Attributes:
  filepath - Directory location of the GSEA output EDB folder.   

  '''
  def __init__(
        self,
        filepath: Union[str,Path]
        ):
    print('~~~~~~ parseGSEA ~~~~~~')
    self.filepath = Path(filepath) if isinstance(filepath, (str, Path)) else None

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # Helper functions
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def edbFileParser(
        self,
        edbData,
        label,
        mode='full'
        ):
    '''
    Helper function to parse required information from edb file output from Windows GSEA software.

    Attributes:
    edbData - edbData extracted from the EDB file. 
    label - Typically the parent folder name, which is the title of the GSEA experiment. 
    mode - (Future functionality) user-defined depth of data extraction.

    Returns:
    edbParserOuptut - Zipped list of all the data extracted from the edbData object

    '''
    modes = ['full'] #'basic'
    if mode not in modes:
      raise ValueError("Invalid mode. Expected one of: %s" % modes)

    if mode == 'full':      
      ## Build outputs containers
      genesetList = []
      nesList = []
      fdrList = []
      labels = []

      ## Iterate through the EDB 
      for i in range(len(edbData)):
        for j in range(len(edbData[i])):
          string = edbData[i][j]
          if 'DTG RANKED_LIST' in string:
            genesetParser = re.search(r'(?<=gene_sets.gmt#).*?(?=" ES=")', string).group(0)
            nesParser = float(re.search(r'(?<=NES=").*?(?=" NP=")', string).group(0))
            fdrParser = float(re.search(r'(?<=FDR=").*?(?=" FWER=")', string).group(0))

            ## If the absolute Normalised Enrichment Score is less than 1, skip the geneset
            if abs(nesParser) < 1:
                continue
            ## Otherwise, append all of the extracted information 
            else:
                genesetList.append(genesetParser)
                nesList.append(nesParser)
                fdrList.append(fdrParser)
                labels.append(label)
          else:
            continue
                    
        edbParserOutput = list(zip(genesetList,nesList,fdrList,labels))
        return edbParserOutput
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def uniqueGenes(
        self,
        geneLists,
        lookup
        ) -> Dict:
    '''
    uniqueGenes takes the list of modules represented within each of the network subgraphs and then accesses the gene lookup table for these modules.
    It returns a dictionary containing modules as keys and unique genes as corresponding values.

    Attributes:
    geneLists - A list of network subgraphs and their associated modules.
    lookup - A lookup table containing all modules and their associated genes.

    Returns:
    geneFull - Dictionary containing modules as keys and uniques as values. 
    
    '''
    geneFull = {}
    for genelist in geneLists:
        cluster = genelist[0]
        modules = genelist[1]
        
        fetch_items = []
        
        for module in modules:
            items = lookup[module]
            fetch_items.append(items)
        
        flatItems = [item for sublist in fetch_items for item in sublist]
        geneUnique = list(set(flatItems))
        
        geneFull[cluster] = geneUnique
    
    return geneFull
  
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def Jaccard(
        self,
        moduleNames,
        geneLists,
        ):
    '''
    Jaccard provides the capability to calculate Jaccard indices for comparing genesets. 
    The basic principle of the Jaccard index is a calculation of the percentage of overlapping elements within two lists.
    The biological assumption is that genesets with overlapping membership will be related to the same biological process, and will thus reduce redundancy. 

    Attributes:
    moduleNames - A list of all modules for comparison.
    geneLists - A list of the unique genes for each of the unique modules.

    Returns:
    jaccardOutput - zipped lists containing the two modules being compared (node1 vs node2) and the calculated Jaccard index. 
    
    '''
    ## Build master and output containers
    masterList = list(zip(moduleNames,geneLists))
    
    node1List = []
    node2List = []
    jaccard = []
  

    ## Use itertools.combinations to iteratively compare all possible pairs in the masterlist 
    for x,y in tqdm(itertools.combinations(masterList, 2), colour='yellow'):
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

  def significantGenesets(
        self,
        zippedList,
        threshold: float
        ):
    '''
    significantGenesets filters modules to only include those that have a false-discovery rate (FDR) below a specified significance threshold (default 0.1).
    
    Attributes:
    zippedList - Zipped output from edbParser function
    threshold - User-defined significance threshold (Default = 0.1)
    '''
    ## This loop goes through the numerical data for the output edbParser list and applies a threshold filter - retain only those <= the statistical threshold
    filteredGenesets = [x for x in zippedList if float(x[2]) <= threshold]
    
    return filteredGenesets

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def fileAccessor(
        self,
        jaccardFilter: float = 0.5,
        statThreshold: float = 0.1
        ):
    ''' 
    Primary fileAccessor function.
    fileAccessor accepts a home directory - specified by the user, in the form of a variable-assigned string - from which the OS module will walk through the available directories and files.
    Specifically, the function searches for the destination folders called << edb >> in which the raw data are kept; all other folders are skipped.
    Once within the << edb >> folder the function extracts the two raw data files (.gmt, .edb)
    The default output path for GSEA follows this same architecture, so the function should be generalisable.
    '''
    print()
    print(f'Retrieving data from {self.filepath}.')
    ## Survey directory trees to find the target folders << edb >>
    destinations = []
    for root, dirs, files in os.walk(self.filepath, topdown=True):
        for folder in dirs:
            if 'edb' not in folder:
                continue
            else:
                destination = str(f'{root}/{folder}')  #os.path.join(root,folder)
                print(f'EDB folder found at:{destination}.')
                destinations.append(destination)

    ## Prepare output object
    output = networkResults
    
    ## Iterate through the target file paths identified above
    for target in destinations:
      print('Parsing data at destinations.')
      time.sleep(2)
      ## At this indentation level, we are within one of the target folder paths, so handling of gmt and edb data unique to one folder - i.e. unique to one pairing of conditions 
      print('Parsing:',target)

      gmtData = []
      edbData = []

      for files in os.listdir(target):
        if '.gmt' in files:
          ## Open and extract data from the .gmt file - append to output container
          fileName = os.path.join(target,files)
          with open(fileName) as gmt_:
            rawGmt = gmt_.read()
            splitGmt = re.split(r'\n', rawGmt)
                    
          for i in range(len(splitGmt)):
            GMT = re.split(r'\t',splitGmt[i])
            gmtData.append(GMT)
                    
        elif '.edb' in files:
          ## Open and extract data from the .edb file - append to output container
          fileName = os.path.join(target,files)
          with open(fileName) as edb_:
            rawEdb = edb_.read()
            EDB = re.findall(r'<(.+?)>',rawEdb)
            edbData.append(EDB)
        else:
            continue
          
      ## First handle the .gmt data
      ## Build output containers for the full list of module names and associated genes from those modules (within .gmt data) and the dictionary for lookup 
      moduleNames = []
      moduleGenes = []
      geneDict = {}
      
      for i in range(len(gmtData)):
        ## Iterate through the .gmt data to extract the geneset/module name and the constituent genes of that module
        name = gmtData[i][0]
        genes = gmtData[i][2:]
        geneDict[name] = genes
        moduleNames.append(name)
        moduleGenes.append(genes)
      
      print('PARSED - .gmt file.')
          
      ## Zip together module names and module genes to avoid mismatching
      modulesZipped = list(zip(moduleNames,moduleGenes))
      
      ## Second handle the .edb data
      ## edbFileParser extracts the statistical data (NES and FDR) from the edb file, with the corresponding geneset name
      statsData = self.edbFileParser(edbData,label=target,mode='full')
      print('PARSED - .edb file.')

      ## significantGenesets receives the zipped statistical data from edbFileParser and filters the full list with an FDR cutoff of <0.1 - it also has the functionality to unzip data
      statsDataThresholded = self.significantGenesets(statsData,statThreshold)
      print('PARSED - Significant genesets.')

      ## Extract geneset names from the significant geneset data and use this to truncate the modulesZipped variable to include only the statistically significant genesets                
      sigGenesets = [x[0] for x in statsDataThresholded]
      sigModules = [x for x in modulesZipped if x[0] in sigGenesets]

      print('Computing Jaccard indices.')
      time.sleep(1)
      ## Calculate the Jaccard indices for the significantly enriched genesets
      jaccardIndices = self.Jaccard([x[0] for x in sigModules],[x[1] for x in sigModules])


      ## Produce output dataframe for data within this loop iteration
      ## Format of this dataframe is a node-node relationship table, which is packaged within the jaccardIndices zipped list
      networkData = pd.DataFrame({'node1':[x[0] for x in jaccardIndices],
                                  'node2':[x[1] for x in jaccardIndices],
                                  'jaccard_index':[x[2] for x in jaccardIndices],
                                  'draw_edge':[1 if x[2] > jaccardFilter else 0 for x in jaccardIndices],
                                  'condition':target})
      ## Format of this dataframe is a statistical information table, which is packaged within the statsDataThresholded zipped list
      metaData = pd.DataFrame({'node':[x[0] for x in statsDataThresholded],
                              'nes':[x[1] for x in statsDataThresholded],
                              'fdr':[x[2] for x in statsDataThresholded],
                              'upregulated?':[1 if x[1] > 0 else 0 for x in statsDataThresholded],
                              'condition':target})

      ## Append data from this loop to the networkResults object 
      output.label = target 
      output.dataframe = networkData
      output.metadata = metaData
      output.statistics = statsData
      output.lookup = geneDict
    
    return output

