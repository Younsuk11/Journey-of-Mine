# Bioinformatcis Intership

> Interval 21.6 ~ 22.3

Programming language : python

## Pathway Analysis Automation

**purpose**
> After DEG analysis, for Gene Set Enrichment Analysis and pathyway analysis
Analyzing "Significant target gene pathway" and "Down or Up -regulated genes" by drug samples!

**Method**
* input data : DEG analysis result, Up and Down file(gene list that sorted by LogFC values)
* Web Crowaling with Chrome Webdriver
* python program : DEG_data_preprocessing.py

**result**
> You can get pathway analysis result tables from Enrichr(https://maayanlab.cloud/Enrichr/)

## CDS Checking

**purpose**
> Mapping CDS regions(with NP, XP id) to RNA(with NM, XM ID)

**method**
* input data : fasta file , annotation files
* python program : parsing_mapping_Final.py

## Scoring siRNA project

**purpose**
> Assessment siRNA sequences by scoring scheme and making candidate efficient siRNA
> Comparing target siRNA and random sequence by scoring

**Method**
* python program : scoring_RNA.py

