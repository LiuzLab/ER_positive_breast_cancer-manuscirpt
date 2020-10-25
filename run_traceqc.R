library(TraceQC)
library(fastqcr)
library(tools)

ref_file <- "./data/ref21.txt"
qc_dir <- "./qc_report/"
fastqc("./data/merged_r1/",qc.dir=qc_dir)
fastq_dir <- "./data/merged_r1/"

for (input_file in list.files("./data/merged_r1/",pattern="LIBILCM")) {
  name <- strsplit(input_file,"\\.")[[1]][1]
  out_alignment_file <- paste("./cache/alignment/",name,".txt",sep="")
  if (!file.exists(out_alignment_file)) {
    print(name)
    sequence_alignment(input_file=paste(fastq_dir,name,".fq.gz",sep=""),
                       ref_file=ref_file,output_file=out_alignment_file,ncores=16)}}

for (input_file in list.files("./cache/alignment/",pattern="LIBILCM")) {
  name <- strsplit(input_file,"\\.")[[1]][1]
  out_file <- paste("./cache/traceQC_obj/",name,".rds",sep="")
  if (!file.exists(out_file)) {
    print(name)
    obj <-
      create_TraceQC_object(
        aligned_reads_file = paste("./cache/alignment/", input_file, sep = ""),
        ref_file = ref_file,
        fastqc_file = paste("./qc_report/", name, "_fastqc.zip", sep =
                              ""),
        alignment_score_threshold = 200,
        abundance_threshold = 10,
        ncores = 1
      )
    saveRDS(obj,file=out_file)
  }}

for (input_file in list.files("./cache/traceQC_obj/")) {
  name <- strsplit(input_file,"\\.")[[1]][1]
  out_file <- paste("/home/humble_local_25t/jasper/projects/CRISPR_barcode/Shawn3/traceQC_report/",name,
                    ".html",sep="")
  if (!file.exists(out_file)) {
    print(name)
    generate_qc_report_from_obj(paste("/home/humble_local_25t/jasper/projects/CRISPR_barcode/Shawn3/cache/traceQC_obj/",
                                      name,".rds",sep=""),
                                out_file,
                                alignment_score_threshold=200)}}