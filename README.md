### GSEA Network Graph

Geneset enrichment analysis (GSEA) is a mainstay of RNA-sequencing analysis workflows. Functional enrichment analyses afford exceptional utility for both unbiased data exploration, and hypothesis driven approaches. However, due to the redundancy of query genesets utilised from the Gene Ontology database, GSEA can output many enriched genesets that overlap significantly in both semantic annotation and gene composition. 

In this context, genesets are represented as nodes. Nodes are evaluated in a pair-wise manner to compute Jaccard indices of their overlapping genes. Edges are drawn between nodes that have a Jaccard index of > 0.5 (i.e. more than 50% of their genes are shared), under the assumption that genesets sharing the majority of their genes will be involved in similar functions. The stringency/leniency of this threshold can be selected at will.

This package contains several utilities for the parsing of data from the Windows installation of GSEA, the generation of relational network graphs from significantly enriched genesets, the identification of interconnected subgraphs, and a tokenisation of the geneset annotations that populate subgraphs of interest. 
By evaluating subgraphs of interrelated nodes, instead of individual genesets the redundancy of the analysis is reduced. 
