### graphGSEA

## Table of contents
- [About](#about)
- [Installation](#installation)
- [Tutorial](#tutorial)
- [Docs](#docs)

## About <a name = "about"></a>
Geneset enrichment analysis (GSEA) is a mainstay of RNA-sequencing analysis workflows. 
Functional enrichment analyses afford exceptional utility for both unbiased data exploration, and hypothesis driven approaches. However, due to the redundancy of query genesets utilised from the Gene Ontology database, GSEA can output many enriched genesets that overlap significantly in both semantic annotation and gene composition. 

graphGSEA seeks to overcome this redundancy by building relational networks of genesets and analysing subgraphs of related genesets. In this context, genesets are represented as nodes. Nodes are evaluated in a pair-wise manner to compute Jaccard indices of their overlapping genes. Edges are drawn between nodes that have a Jaccard index of > 0.5 (i.e. more than 50% of their genes are shared), under the assumption that genesets sharing the majority of their genes will be involved in similar functions. The stringency/leniency of this threshold can be selected at will.

This package contains several utilities for the parsing of data from the Windows installation of GSEA, the generation of relational network graphs from significantly enriched genesets, the identification of interconnected subgraphs, and a tokenisation of the geneset annotations that populate subgraphs of interest. 
By evaluating subgraphs of interrelated nodes, instead of individual genesets, the redundancy of the analysis is reduced. 

## Installation <a name = "installation"></a>
```python

pip install graphGSEA

```
## Tutorial <a name = "tutorial"></a>
Import the parseGSEA function from the parsing utility folder. Use the fileAcccessor() function to extract data from the designated GSEA experiment directory containing the EDB folder. 

```python

from parsing import parseGSEA

outputs = parseGSEA(filepath='./Example_data').fileAccessor(jaccardFilter=0.5,statThreshold=0.1)


## ~~~~~~ parseGSEA ~~~~~~

## Retrieving data from Example_data.
## EDB folder found at:Example_data\GSEA_7m/edb.
## Parsing data at destinations.
## Parsing: Example_data\GSEA_7m/edb
## PARSED - .gmt file.
## PARSED - .edb file.
## PARSED - Significant genesets.
## Computing Jaccard indices.
## 2552670it [00:28, 91116.15it/s] 
```
Import the plotGSEA function from the plotting utility folder. Pass the outputs from parseGSEA to the iterativePlotting() function to plot network graphs from the parsed GSEA data. 

```python
from parsing import parseGSEA
from plotting import plotGSEA

outputs = parseGSEA(filepath='./Example_data').fileAccessor(jaccardFilter=0.5,statThreshold=0.1)
plotGSEA.iterativePlotting(outputs,save=False)
```

![alt text](https://github.com/MattMason95/GSEA-Network-Graph/blob/graphGSEA/GSEA_7m.png?raw=true)
