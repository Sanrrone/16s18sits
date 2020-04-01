function qc {
	set -e
	echo "ESCLAVO: QC begin"
	if [ ! -d 1-qc ];then
		mkdir 1-qc
	fi
	
	cd 1-qc

	if [ $FORCE ];then
		rm -rf *
	fi

	sampleR1=$(ls -1 $FASTQFOLDER/*$PATTERN | head -n1)
	total=$(if [[ "$PATTERN" =~ "gz" ]] || [[ "$PATTERN" =~ "zip" ]]; then zcat $sampleR1 ;else cat $sampleR1 ;fi | wc -l |awk '{print int($1/4)}') #get total fastq sequences
	lengthSeq=$(if [[ "$PATTERN" =~ "gz" ]] || [[ "$PATTERN" =~ "zip" ]]; then zcat $sampleR1 ;else cat $sampleR1 ;fi | awk -v total=$total '{if(NR%4==2 && NR<(total/4)) print length($1)}' | sort -n | uniq -c |tail -n1 |awk '{print $2}')
	nfiles=$(ls -1 $FASTQFOLDER/*${PATTERN} |wc -l |awk '{print $1}')

	echo "timeElapsed" > tmp0
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "0:0:0"}' >> tmp0
	echo "inputFiles" > tmp1
	ls -1 $FASTQFOLDER/*${PATTERN} >> tmp1
	echo "stepStatus" > tmp2
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "running"}' >> tmp2
	paste tmp0 tmp1 tmp2 > qc.conf && rm tmp0 tmp1 tmp2
	
	echo "
	rm(list=ls())
	if('BiocManager' %in% rownames(installed.packages()) == FALSE) {install.packages('BiocManager')}
	if('dada2' %in% rownames(installed.packages()) == FALSE) {BiocManager::install('dada2')}
	library(dada2)
	options(scipen=999)
	
        getr1r2<-function(path,fqpattern){
                allfiles<-list.files(path, pattern=paste0('*',fqpattern), full.names = TRUE)

                r1<-sort(list.files(path, pattern=paste0('[\\\\.|\\\\_]R1[\\\\.|\\\\_]?.*',fqpattern), full.names = TRUE))
                r2<-sort(list.files(path, pattern=paste0('[\\\\.|\\\\_]R2[\\\\.|\\\\_]?.*',fqpattern), full.names = TRUE))

                if(length(r1)==length(r2) & length(r1)+length(r2) == length(allfiles)){
                        return(list(r1,r2,'[\\\\.|\\\\_]R1'))

                }

                r1<-sort(list.files(path, pattern=paste0('[\\\\.|\\\\_]1',fqpattern), full.names = TRUE))
                r2<-sort(list.files(path, pattern=paste0('[\\\\.|\\\\_]2',fqpattern), full.names = TRUE))

                if(length(r1)==length(r2) & length(r1)+length(r2) == length(allfiles)){

                        return(list(r1,r2,'[\\\\.|\\\\_]1\\\\.f'))
                }else{

                        r1<-sort(list.files(path, pattern=paste0('[\\\\.|\\\\_]1[\\\\.|\\\\_]?.*',fqpattern), full.names = TRUE))
                        r1<-r1[!grepl(x=r1, pattern = '[\\\\.|\\\\_]R?2[\\\\.|\\\\_]', ignore.case = T)]
                        r2<-sort(list.files(path, pattern=paste0('[\\\\.|\\\\_]2[\\\\.|\\\\_]?.*',fqpattern), full.names = TRUE))
                        r2<-r2[!grepl(x=r2, pattern = '[\\\\.|\\\\_]R?1[\\\\.|\\\\_]', ignore.case = T)]
                        return(list(r1,r2,'[\\\\.|\\\\_]R?1'))



                }
        }


	path<-'$FASTQFOLDER'
	fqpattern<-'$(echo "$PATTERN" |sed 's/\./\\\\./g')' #double \\ to escape .
	tolerance<-$TOLERANCE
	readlength<-$lengthSeq
	projectfolder<-'$PROJECTFOLDER'
	paired<-$PAIRED
	
	if(paired){
		r1r2<-getr1r2(path,fqpattern)
		fnFs<-r1r2[[1]]
		fnRs<-r1r2[[2]]

		sample.names <- sapply(strsplit(basename(fnFs), r1r2[[3]]), \`[\`, 1)

		# Place filtered files in filtered/ subdirectory
		filtFs <- file.path(projectfolder, '1-qc', paste0(sample.names, '_F_filt','$PATTERN'))
		filtRs <- file.path(projectfolder, '1-qc', paste0(sample.names, '_R_filt','$PATTERN'))
		names(filtFs) <- sample.names
		names(filtRs) <- sample.names
        	print(paste('ESCLAVO: working on ',fnFs,' ',fnRs))
		readtolerance<-readlength*tolerance
		maxeeformula<- (0.01*readlength)+(0.012589254*readtolerance)
		print(paste0('ESCLAVO: Doing filtering at maxEE ',maxeeformula))
		out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=(readlength-readtolerance),
	              maxN=0, maxEE=maxeeformula, truncQ=2, rm.phix=TRUE, minLen = 50,
	              compress=FALSE, multithread=TRUE, verbose=TRUE)
	}else{
		fnFs <- sort(list.files(path, pattern=fqpattern, full.names = TRUE)) 
		sample.names <- sapply(strsplit(basename(fnFs), '\\\\.fastq'), \`[\`, 1)
		filtFs <- file.path(projectfolder, '1-qc', paste0(sample.names, '_F_filt','$PATTERN'))
		names(filtFs) <- sample.names 
                print(paste('ESCLAVO: working on ',fnFs))
		readtolerance<-readlength*tolerance
                maxeeformula<- (0.01*readlength)+(0.012589254*readtolerance)
                print(paste0('ESCLAVO: Doing filtering at maxEE ',maxeeformula))
                out <- filterAndTrim(fnFs, filtFs, truncLen=(readlength-readtolerance),
                      maxN=0, maxEE=maxeeformula, truncQ=2, rm.phix=TRUE, minLen = 20,
                      compress=FALSE, multithread=TRUE, verbose=TRUE)


	}
	write.table(out,'qc_filt.tsv',sep='\t')

	" > dada2_filt.R
	SECONDS=0
	Rscript --vanilla dada2_filt.R
	duration=$SECONDS
	echo "timeElapsed" > tmp0
	echo $duration | awk -v nfiles=$nfiles -v duration=$duration '{for(i=1;i<=nfiles;i++){print int($1/60/60/nfiles)":"int($1/60/nfiles)":"($1%60)/nfiles}}' >> tmp0
	echo "inputFiles" > tmp1
	ls -1 $FASTQFOLDER/*${PATTERN} >> tmp1
	echo "stepStatus" > tmp2
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "done"}' >> tmp2
	paste tmp0 tmp1 tmp2 > qc.conf && rm tmp0 tmp1 tmp2
	rm dada2_filt.R

	if [ "$PCONF" != "" ]; then
		echo "ESCLAVO: Updating config file: $PCONF"
		sed -i "s/running/done/g" qc.conf
		sed -i "s/pPercent.*/pPercent\t40/g" $PCONF
		sed -i "s/lastStep.*/lastStep\tQC/g" $PCONF
	fi

	echo "ESCLAVO: QC end"
}
