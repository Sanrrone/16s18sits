# 16s18sits (16S 18S and ITS)
------------------------------------------------------------------
Welcome to 16s18sits (and you are welcome to suggest a better name :D ), 16S  pipeline based on [dada2](https://benjjneb.github.io/dada2/tutorial.html) software including analysis, summaries and reports in html and pdf.

Actually the pipeline always use all your cores (depending the step) and only support 16S and we are working to include 18S and ITS, sorry about that.

# Output

# Usage

The pipeline have a **basic** usage (only with three parameters):
`bash 16s18sits.sh -p [project folder] -f [fastq foder] -pt [fastq pattern]`
`example: bash 16s18sits.sh -p mynewproject -f 0-raw -pt .fastq.gz`

and the **complete** usage:
`bash 16s18sits.sh -p [project folder] -f [fastq foder] -pt [fastq pattern]`
`example: bash 16s18sits.sh -p mynewproject -f 0-raw -pt .fastq.gz`

options available:

* -p: Project folder, if no exist, the Pipeline will assume is the actual folder.
* -f: fastq folder, mandatory parameter.
* -pt: fastq pattern, it could be '.fastq' or '.fastq.gz'.
* -to: read Tolerance for errors (default 0.1 in a range of 0 to 1). It means you tolerate a maximum of 10% of errors in your reads (QC < 20)
* -force: overwrite existing pipeline runs.
* 
* -h: print this help.

more complete readme soon