# Sequencing Run Quality Control

## Nanopore Data

### Run Report

Things to check:
1. Amount of data generated
2. Software versions
3. Basecalling settings
4. Q score over time
5. Pore activity
6. Translocation speed

#### Amount of data generated
Our sequencing runs have generated between 5-10 giga bases of data. This is an important parameter to check on the run report, as we want to maximise the output per flow cell. 

Methods to maximise output:
1. Flow cell washing and reloading library - this reduces the number of inactive pores increasing the amount of data generated
2. Check trouble shooting sections in the run report


#### Software version
Talks with nanopore representatives made it clear the importance of using up-to-date basecalling software to maximise the read classification. At the time of writing this is "Dorado 7.3.9", and the software should be frequently updated.
#### Basecalling settings
Multiple tests showed Super-accurate basecalling (AKA super high accuracy). This can result in lower amounts of data being classified but accuracy base calling. As read accuracy is crucial to assembly of RAGEs we deamed this an important parameter to use.

### PycoQC

## Read Ratios

