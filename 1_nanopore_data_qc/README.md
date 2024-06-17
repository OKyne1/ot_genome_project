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

Also check all the troubleshooting information provided, this is often quite comprihensive and suggests protential issues.

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
[PycoQC]([https://pages.github.com/](https://github.com/a-slide/pycoQC)) "computes metrics and generates interactive QC plots for Oxford Nanopore technologies sequencing data" and is used as an initial check of data quality from a sequencing run.

#### Usage:
Needs to be run from an environment with pycoQC installed
```
pycoQC --summary_file ./path/sequencing_summary.txt --html_outfile HTML_OUTFILE --min_pass_qual MIN_PASS_QUAL
```
The correct min pass qual needs to be selected (otherwise the results will be wrong and confusing), refer back to the run-report for this value. If basecalling was rerun, ensure that the correct run-report is chosen.

#### Splitting the summary file:
```
Barcode_split --summary_file path/to/file/sequencing_summary.txt --output_dir OUTPUT_DIR --output_unclassified
```
This will produce separate files for each barcode, these then need to be run with pycoQC to generate an output (e.g. [pycoqc_barcode_12.html](./pycoqc12.html))

#### Things to check:
1. Read length distributions (all and individual barcodes), we have been trying to maximise reads >10,000.
2. Quality scores, how much data is lost due? Are some barcodes worse than others?
3. Reads per barcode, how even are these and how many unclassified reads are there? Often over representation of a barcode could suggest higher DNA fragmentation in this sample.
4. Channel activity over time, the more red/yellow the better. Low activity in many of the pores may suggest the presence of bubbles.

## Read Ratios
### Joining files
### Mapping to mouse (+ lambda phage)
### Quantification
