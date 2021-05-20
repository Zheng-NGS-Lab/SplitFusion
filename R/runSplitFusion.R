# 
runSplitFusion <- function(){

options(width=204)
source('config.txt')

dir.create(paste0(output,"/",sample_id), recursive = TRUE, showWarnings = FALSE)
setwd(paste0(output,"/",sample_id))

cat("==========\n", date(), ": SplitFusion started...\n==========\n")

##==== Step 1: fastq to consolidated bam (de-duplication based on unique UMI and adatptor ligation site) ====
		cat("==========\n", date(), ": 1_fastq-bam...\n==========\n")
		steps=gsub(" ", "", steps)
	if ("1_fastq-bam" %in% unlist(strsplit(steps, ","))){
		suppressMessages(system(paste0("bash ", SFpath, "/exec/SplitFusion.1_fastq-bam.sh ", "SA")))
		cat("==========\n", date(), ": 1_fastq-bam completed.\n==========\n")
	}

##==== 2: consolidated.bam to breakpoint candidates (no filter) ====
		cat("==========\n", date(), ": 2_bam-breakpoint...\n==========\n")
	if ("2_bam-breakpoint" %in% unlist(strsplit(steps, ","))) {
		suppressMessages(system(paste0("bash ", SFpath, "/exec/SplitFusion.2_bam-breakpoint.sh ")))
		cat("==========\n", date(), ": 2_bam-breakpoint completed.\n==========\n")
	}

##==== 3: filters ====
		cat("==========\n", date(), ": 3_breakpoint-filter...\n==========\n")
	if ("3_breakpoint-filter" %in% unlist(strsplit(steps, ","))) {
		suppressMessages(system(paste0("bash ", SFpath, "/exec/SplitFusion.3_breakpoint-filter.sh ")))
		cat("==========\n", date(), ": 3_breakpoint-filter completed.\n==========\n")
	}

##==== 4: Annotate breakpoint gene, exon, cDNA position
		cat("==========\n", date(), ": 4_breakpoint-anno...\n==========\n")
	if ("4_breakpoint-anno" %in% unlist(strsplit(steps, ","))) {
        	suppressMessages(system(paste0("bash ", SFpath, "/exec/SplitFusion.4_breakpoint-anno.sh ")))
		cat("==========\n", date(), ": 4_breakpoint-anno completed.\n==========\n")
	}

##==== 5: further processing (in-frame status determination, etc)
		cat("==========\n", date(), ": 5_breakpoint-anno-post...\n==========\n")
	if ( "5_breakpoint-anno-post" %in% unlist(strsplit(steps, ","))) {
		source(paste0(SFpath,"/R/SplitFusion.breakpoint.anno.postscript.R"))
		suppressMessages(SplitFusion.breakpoint.anno.postscript())
		cat("==========\n", date(), ": 5_breakpoint-anno-post completed.\n==========\n")
	}

#        system("rm _*")
if (file.exists(paste0(sample_id, ".brief.summary"))){
	cat("==========\n", date(), ": Split-Fusion completed successfully.\n==========\n")
}else{
	cat("==========\n", date(), ": Split-Fusion completed with ERROR.\n==========\n")
}
## END
}
