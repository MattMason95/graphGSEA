# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import libraries 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pandas as pd 
import numpy as np
import regex as re
import re
import os
import itertools
import networkx as nx

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def edbFileParser(edbData,label,mode='full'):
  '''
  This function will take the full array of EDB data and extract the geneset identifier, normalised enrichment scores, and FDR statistics using Regex pattern recognition
  May build up Full versus Basic modes to modify the extracted specific information. 
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

def fileAccessor(homeDirectory,jaccardFilter):
  ''' 
  This is the primary function that integrates all above functions. 
  fileAccessor accepts a home directory - specified by the user, in the form of a variable-assigned string - from which the OS module will walk through the available directories and files.
  Specifically, the function searches for the destination folders called << edb >> in which the raw data are kept; all other folders are skipped.
  Once within the << edb >> folder the function extracts the two raw data files (.gmt, .edb)
  The default output path for GSEA follows this same architecture, so the function should be generalisable.
  '''
  ## Survey directory trees to find the target folders << edb >>
  destinations = []
  for root, dirs, files in os.walk(homeDirectory, topdown=True):
      for folder in dirs:
          print(folder)
          if 'edb' not in folder:
              continue
          else:
              destination = str(f'{root}/{folder}')  #os.path.join(root,folder)
              print('destination:',destination)
              destinations.append(destination)
  
  ## Print out the identified file paths
  print(destinations)

  ## Build output containers 
  networkTemp = []
  metaTemp = []
  masterStats = []
  masterGeneLookup = {}
  
  ## Iterate through the target file paths identified above
  for target in destinations:
    ## At this indentation level, we are within one of the target folder paths, so handling of gmt and edb data unique to one folder - i.e. unique to one pairing of conditions 
    print('target:',target)
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
    
    ## Update out-of-loop dictionary for gene lookup
    masterGeneLookup.update(geneDict)
        
    ## Zip together module names and module genes to avoid mismatching
    modulesZipped = list(zip(moduleNames,moduleGenes))
    
    ## Second handle the .edb data
    ## edbFileParser extracts the statistical data (NES and FDR) from the edb file, with the corresponding geneset name
    statsData = edbFileParser(edbData,label=target,mode='full')

    ## significantGenesets receives the zipped statistical data from edbFileParser and filters the full list with an FDR cutoff of <0.1 - it also has the functionality to unzip data
    statsDataThresholded = significantGenesets(statsData,0.1,unzip=False)

    ## Extract geneset names from the significant geneset data and use this to truncate the modulesZipped variable to include only the statistically significant genesets                
    sigGenesets = [x[0] for x in statsDataThresholded]
    sigModules = [x for x in modulesZipped if x[0] in sigGenesets]

    ## Calculate the Jaccard indices for the significantly enriched genesets
    jaccardIndices = Jaccard([x[0] for x in sigModules],[x[1] for x in sigModules])


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

    ## Append data from this loop to the output containers
    networkTemp.append(networkData)
    metaTemp.append(metaData)
    masterGeneLookup.append(geneLookup)
    masterStats.append(statsData)

  ## Concatenate output dataframes
  masterNetwork = pd.concat(networkTemp)
  masterMeta = pd.concat(metaTemp)
      
  
  return masterNetwork, masterMeta, masterStats, masterGeneLookup
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def uniqueGenes(geneLists,lookup):
    '''
    This function takes the list of modules represnted within each of the network subgraphs and then accesses the gene lookup table for these modules. 
    All unique genes for a given network cluster are returned as dictionary.
    '''

    geneFull = {}
    for list_ in geneLists:
        cluster = list_[0]
        modules = list_[1]
        
        fetch_items = []
        
        for module in modules:
            items = lookup[module]
            fetch_items.append(items)
        
        flatItems = [item for sublist in fetch_items for item in sublist]
        geneUnique = list(set(flatItems))
        
        geneFull[cluster] = geneUnique
    
    return geneFull
