function humanDecont {
	set -e
	echo "ESCLAVO: humanDecont begin"

	QCFOLDER=$1

	if [ "$(which centrifuge)" == "" ];then
		echo "ESCLAVO: no centrifuge installed, aborting this step"
		exit
	fi

	if [ ! -d 2-decont ];then
		mkdir 2-decont
	fi

	if [ $FORCE ];then
		rm -r 2-decont
		mkdir 2-decont
	fi

	cd 2-decont

	if [ ! -f "$ESCLAVOHOME/DB/humanDB.1.cf" ];then
		actualpath=$(pwd)
		cd $ESCLAVOHOME/DB
		#humanDB was created using centrifuge v1.0.4(beta)
		echo "No humanDB found, donwloading..."
		wget --no-check-certificate -r 'https://docs.google.com/uc?export=download&id=11pvlewFGBJXpFvL6tvmZACkf-g3d8nhW' -O humanDB.1.cf
		wget --no-check-certificate -r 'https://docs.google.com/uc?export=download&id=1AZHcGrWqRTFl48X4H20PvaB7bw-Z0sKT' -O humanDB.2.cf
		wget --no-check-certificate -r 'https://docs.google.com/uc?export=download&id=1nVT1Q7KYd-9gFfjHj_6a0qwN1XwBsthW' -O humanDB.3.cf
		wget --no-check-certificate -r 'https://docs.google.com/uc?export=download&id=1rNJw9hiXDCamoGIG3PWi2tYvFn9VZgrz' -O humanDB.4.cf
		echo "Done"
		cd $actualpath
	fi

	if [ "$PAIRED" == "TRUE" ];then	

		ls ../$QCFOLDER/*_F_filt$PATTERN > r1
		ls ../$QCFOLDER/*_R_filt$PATTERN > r2

		paste r1 r2 > r1r2 && rm -f r1 r2
	else
		ls ../$QCFOLDER/*_F_filt$PATTERN > r1r2
	fi
	#### status step file
	nfiles=$(ls -1 $FASTQFOLDER/*${PATTERN} |wc -l |awk '{print $1}')
	echo "timeElapsed" > tmp0
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "0:0:0"}' >> tmp0
	echo "inputFiles" > tmp1
	ls -1 ../$QCFOLDER/*${PATTERN} >> tmp1
	echo "stepStatus" > tmp2
	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "running"}' >> tmp2
	paste tmp0 tmp1 tmp2 > decont.conf && rm tmp0 tmp1 tmp2
	##################################################################	

	SECONDS=0
	cat r1r2 | while read line
	do
		echo "ESCLAVO: decontaminating $line"
		samplename=$(echo $line | awk '{print $1}' | awk -F"/" '{print $NF}' |awk -F"_F_filt" '{print $1}')
		if [ "$PAIRED" == "TRUE" ];then

			readsf1=$(echo $line | awk '{print $1}')
			readsf2=$(echo $line | awk '{print $2}')
			centrifuge -x $ESCLAVOHOME/DB/humanDB -1 $readsf1 -2 $readsf2 --al-conc ${samplename}_cont --un-conc ${samplename}_decont -p $(nproc) > ${samplename}_centrifuge_report_query.tsv
			mv centrifuge_report.tsv ${samplename}_centrifuge_report_reference.tsv
			## subsampling for learnErrors
			head -n $(echo $nfiles |awk '{print int(4*100000/$1)}') ${samplename}_decont.1 > ${samplename}_subR1.fastq
			head -n $(echo $nfiles |awk '{print int(4*100000/$1)}') ${samplename}_decont.2 > ${samplename}_subR2.fastq
			##compressing
			gzip ${samplename}_cont.[12]
			gzip ${samplename}_decont.[12]
			mv ${samplename}_cont.1.gz ${samplename}_cont.1.fastq.gz
			mv ${samplename}_cont.2.gz ${samplename}_cont.2.fastq.gz
			mv ${samplename}_decont.1.gz ${samplename}_decont.1.fastq.gz
			mv ${samplename}_decont.2.gz ${samplename}_decont.2.fastq.gz
		else
                        readsf1=$(echo $line | awk '{print $1}')
                        centrifuge -x $ESCLAVOHOME/DB/humanDB -U $readsf1 --al ${samplename}_cont --un ${samplename}_decont -p $(nproc) > ${samplename}_centrifuge_report_query.tsv
                        mv centrifuge_report.tsv ${samplename}_centrifuge_report_reference.tsv
                        ## subsampling for learnErrors
                        head -n $(echo $nfiles |awk '{print int(4*100000/$1)}') ${samplename}_decont > ${samplename}_sub.fastq
                        ##compressing
                        gzip ${samplename}_cont
                        gzip ${samplename}_decont
                        mv ${samplename}_cont.gz ${samplename}_cont.fastq.gz 
                        mv ${samplename}_decont.gz ${samplename}_decont.fastq.gz

		fi

	done
	
	export PATTERN=".fastq.gz"
	
	echo "rm(list=ls())
	if('BiocManager' %in% rownames(installed.packages()) == FALSE) {install.packages('BiocManager')}
	if('dplyr' %in% rownames(installed.packages()) == FALSE) {install.packages('dplyr')}
	if('dada2' %in% rownames(installed.packages()) == FALSE) {BiocManager::install('dada2')}

	library(dada2)
	library(dplyr)

	paired<-"$PAIRED"

	print('ESCLAVO: reading fastqs')
	if(paired){
	
		filtFs <- sort(list.files('.', pattern='_decont.1.fastq.gz', full.names = TRUE))
		filtRs <- sort(list.files('.', pattern='_decont.2.fastq.gz', full.names = TRUE))
		subfiltFs<- sort(list.files('.', pattern='_subR1.fastq', full.names = TRUE))
		subfiltRs<- sort(list.files('.', pattern='_subR2.fastq', full.names = TRUE))

		sample.names <- sapply(strsplit(basename(filtFs), '_decont'), \`[\`, 1)
		names(filtFs) <- sample.names
		names(filtRs) <- sample.names
	
		print('learning from errors (subsampled at 400000 reads)')
		errF <- learnErrors(subfiltFs, multithread=TRUE, randomize = TRUE)
		errR <- learnErrors(subfiltRs, multithread=TRUE, randomize = TRUE)

		write.table(errF,'dada2_filt_errF.tsv',sep='\t')
		write.table(errR,'dada2_filt_errR.tsv',sep='\t') 

		print('ESCLAVO: Dereplication')
		derepFs <- derepFastq(filtFs, verbose=TRUE)
		derepRs <- derepFastq(filtRs, verbose=TRUE)

		print('ESCLAVO: Inferring sample composition')
		dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
		dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
	
		print('ESCLAVO: Finally merge reads')
		mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)	
		seqtab <- makeSequenceTable(mergers)

	}else{
                filtFs <- sort(list.files('.', pattern='_decont.fastq.gz', full.names = TRUE))
                subfiltFs<- sort(list.files('.', pattern='_sub.fastq', full.names = TRUE)) 

                sample.names <- sapply(strsplit(basename(filtFs), '_decont'), \`[\`, 1)
                names(filtFs) <- sample.names
        
                print('learning from errors (subsampled at 400000 reads)')
                errF <- learnErrors(subfiltFs, multithread=TRUE, randomize = TRUE)

                write.table(errF,'dada2_filt_errF.tsv',sep='\t')

                print('ESCLAVO: Dereplication')
                derepFs <- derepFastq(filtFs, verbose=TRUE)

                print('ESCLAVO: Inferring sample composition')
                dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
        
                print('ESCLAVO: Finally merge reads')
                mergers <- mergePairs(dadaFs, derepFs, verbose=TRUE)   
                seqtab <- makeSequenceTable(mergers)
	}
	#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)]) # to select specific seq length
	seqtab.nochim <- removeBimeraDenovo(seqtab, method='consensus', multithread=TRUE, verbose=TRUE)
	print(paste0('sequences kept after removing chimera step: ',sum(seqtab.nochim)/sum(seqtab)))

	#make summary table
	getN <- function(x) sum(getUniques(x))
	if(length(filtFs)==1){
	  dadaFs<-list(dadaFs)
	  if(paired){
	  	dadaRs<-list(dadaRs)
	  }
	  mergers<-list(mergers)
	}
	rawdf<-read.table('../$QCFOLDER/qc_filt.tsv',sep = '\t',row.names = 1)
	decontdf<-lapply(names(filtFs), function(x){
	  system(paste0('wc -l ',x,\"_centrifuge_report_query.tsv |awk '{print $1 - 1}'\"),intern = T,wait = T)
	})
	decontdf<-data.frame(decont=unlist(decontdf))
	
	if(paired){
		track <- cbind(rawdf, decontdf, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
	}else{
		track <- cbind(rawdf, decontdf, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

	}
	# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
	if(paired){
		colnames(track) <- c('input','QC','decont','denoisedF', 'denoisedR', 'merged', 'nonchim')
	}else{
		colnames(track) <- c('input','QC','decont','denoised', 'merged', 'nonchim')
	}
	rownames(track) <- sample.names

	write.table(track,'decont_summary.tsv', sep='\t')
	save(seqtab.nochim,file = 'seqtab.nochim.RData')
	" > dada2_centrigugeDecont.R

	echo "ESCLAVO: humanDecont end"
#	echo "ESCLAVO: continuing dada2 mapp preprocess"
#	Rscript --vanilla dada2_centrigugeDecont.R
#	duration=$SECONDS
#	echo "timeElapsed" > tmp0
#	echo $duration | awk -v nfiles=$nfiles -v duration=$duration '{for(i=1;i<=nfiles;i++){print int($1/60/60/nfiles)":"int($1/60/nfiles)":"($1%60)/nfiles}}' >> tmp0
#	echo "inputFiles" > tmp1
#	ls -1 ../$QCFOLDER/*${PATTERN} >> tmp1
#	echo "stepStatus" > tmp2
#	echo $nfiles |awk '{for(i=1;i<=$1;i++)print "done"}' >> tmp2
#	paste tmp0 tmp1 tmp2 > decont.conf && rm tmp0 tmp1 tmp2
#	rm dada2_centrigugeDecont.R
#
#	if [ "$PCONF" != "" ]; then
#		echo "ESCLAVO: Updating config file: $PCONF"
#		sed -i "s/running/done/g" decont.conf
#		sed -i "s/pPercent.*/pPercent\t60/g" $PCONF
#		sed -i "s/lastStep.*/lastStep\tQC/g" $PCONF
#	fi

}

