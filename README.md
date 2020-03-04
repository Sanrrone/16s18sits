# 16s18sits (16S 18S and ITS)
------------------------------------------------------------------
Welcome to 16s18sits (and you are welcome to suggest a better name :D ), 16S  pipeline based on [dada2](https://benjjneb.github.io/dada2/tutorial.html) software including analysis, summaries and reports in html and pdf.

Actually the pipeline always use all your cores (depending the step) and only support 16S and we are working to include 18S and ITS, sorry about that.

# Requisites

* FastQC
* MultiQC
* R with the following libraries:
	* dada2
	* phyloseq
	* dplyr
	* rmarkdown
	* knitr
	* Kable
* Centrifuge == v1.0.4(beta) (for decontamination step)

# Output
---------------------------------
	
The pipeline will create 3 folders (one per module step):
* **1-qc**: the folder contains the reads that passes the QC filters according to dada2 QC functions. and the file qc_filt.tsv that is a summary for raw and filtered reads.
* **2-decont**: the folder contains the reads that passes the Decontamination (human) step, it keeps the reads that mapped again human genome and the no mapped (decontaminated) reads. Additionally it have multiQC report for decontaminated reads in `multiqc_report.html`
* **3-taxInsight**: the last folder of the pipeline where you can find the final report for mapped reads against silvaDF through dada2 steps. The file you are looking for is `abundance.tsv` where contains the abundance for each organism mapped by your reads. Also it contains a PDF image with your sample and organism found.

# Usage
---------------------------------

The pipeline have a **basic** usage (only with three parameters):
`bash 16s18sits.sh -p [project folder] -f [fastq foder] -pt [fastq pattern]`
`example: bash 16s18sits.sh -p mynewproject -f 0-raw -fp .fastq.gz`

and the **complete** usage:
`bash 16s18sits.sh -p [project folder] -f [fastq foder] -fp [fastq pattern] -to [tolerance] - -module [steps to run] -force`
`example: bash 16s18sits.sh -p mynewproject -f 0-raw -fp .fastq.gz -to 0.1 -module all -force`

options available:

* -p: Project folder, if no exist, the Pipeline will assume is the actual folder.
* -f: fastq folder, mandatory parameter.
* -fp: fastq pattern, it could be '.fastq' or '.fastq.gz'.
* -to: read Tolerance for errors (default 0.1 in a range of 0 to 1). It means you tolerate a maximum of 10% of errors in your reads (QC < 20).
* -ler: numnber of reads for learnError dada2 function (50000 by default).
* -force: overwrite existing pipeline runs.
* -h: print this help.
* --debug: show all code step while the pipeline runs.

more complete readme soon
