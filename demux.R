library(Biostrings)
library(dplyr)
library(purrr)
library(ShortRead)

# Define the data
# make gfp sequence
gfp <- Biostrings::DNAString('GTCGCCACCATGGTGAGCAAGGGCGAGGA')
rc_gfp <- Biostrings::reverseComplement(gfp)
# Define a alignment scoring matrix
match <- 1
mismatch <- -1
type <- 'local'
gapOpening <- 0
gapExtension <- -1
submat <- Biostrings::nucleotideSubstitutionMatrix(match=match,
                                                   mismatch=mismatch,
                                                   baseOnly=TRUE)



  fqFile <- readFastq(FastqFile("/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/1_aligned/gfp_corrected/reverse/reverse_mapped.fastq"))
  fastq <- sread(fqFile)
  print(head(fastq))

  # polyA in cDNA
  # Barcode is in the begining of FastQ as in polyA
  # Barcode is not reversed or complemented

  num_bases_to_check_gfp_in <- 100
  fastq_start <- substr(fastq,
                        start=1,
                        stop=num_bases_to_check_gfp_in)
  fastq_end <- substr(fastq,
                      start=(nchar(fastq)-num_bases_to_check_gfp_in),
                      stop=nchar(fastq))
  read_length <- nchar(fastq)

  #Pairwise alignment iterating through each element within the DNAStringSet, comparing to gfp as a string
  for(j in 1:length(fastq)){
    cat("Iteration #",j)
    as_gfp <- Biostrings::pairwiseAlignment(pattern=toString(gfp),
                                            subject=fastq_start[j],
                                            substitutionMatrix=submat,
                                            type=type,
                                            scoreOnly=FALSE,
                                            gapOpening=gapOpening,
                                            gapExtension=gapExtension)

    as_rc_gfp <- Biostrings::pairwiseAlignment(pattern=toString(rc_gfp),
                                               subject=fastq_end[j],
                                               substitutionMatrix=submat,
                                               type=type,
                                               scoreOnly=FALSE,
                                               gapOpening=gapOpening,
                                               gapExtension=gapExtension)

    #Normalizing forward and reverse complementary alignment scores
    nas_gfp <- as_gfp@score/gfp@length
    nas_rc_gfp <- as_rc_gfp@score/rc_gfp@length
    
    #determine forward, reverse or invalid read from best normalized alignment score
    if (nas_gfp > 0.5 & (nas_gfp > nas_rc_gfp)) {
      read_type <- 'forward'
    } else if (nas_rc_gfp > 0.5 & (nas_rc_gfp > nas_gfp)) {
      read_type <- 'reverse'
      print("found")
    } else {
      read_type <- 'invalid'
      print("invalid read, continuing")
      next()
    }
    #if(j==20){
    #  stop()
    #}
    #next()

    if (read_type == 'forward') {
      barcode_end <- as_gfp@subject@range@start #locate start of gfp sequence corresponding to end of barcode

      #limit barcode search to 20 positions before gfp start
      fastq_to_search_barcode_in <- substr(fastq[j], start=barcode_end-20, stop=barcode_end)

      #Align read to each of five barcodes as they appear in forward reads
      as_all <- Biostrings::pairwiseAlignment(pattern=c('TGCAGATCTCTTGCC',
                                                        'CTCGAAGCATTGTAA',
                                                        'AACGGTAGCCACCAA',
                                                        'TGCACGAGATTGATG',
                                                        'GACACATAGTCATGG'),
                                              subject=fastq_to_search_barcode_in,
                                              substitutionMatrix=submat,
                                              type=type,
                                              scoreOnly=FALSE,
                                              gapOpening=gapOpening,
                                              gapExtension=gapExtension)


    } else if (read_type == 'reverse') {
      barcode_start <- width(fastq[j])-num_bases_to_check_gfp_in+as_rc_gfp@subject@range@start+as_rc_gfp@subject@range@width-1 

      fastq_to_search_barcode_in <- substr(fastq[j], start=barcode_start, stop=nchar(fastq[j]))

      as_all <- Biostrings::pairwiseAlignment(pattern=c('GGCAAGAGATCTGCA',
                                                        'TTACAATGCTTCGAG',
                                                        'TTGGTGGCTACCGTT',
                                                        'CATCAATCTCGTGCA',
                                                        'CCATGACTATGTGTC'),
                                              subject=fastq_to_search_barcode_in,
                                              substitutionMatrix=submat,
                                              type=type,
                                              scoreOnly=FALSE,
                                              gapOpening=gapOpening,
                                              gapExtension=gapExtension)
    }

    nas_bc1 <- as_all[1]@score/15
    nas_bc2 <- as_all[2]@score/15
    nas_bc3 <- as_all[3]@score/15
    nas_bc4 <- as_all[4]@score/15
    nas_bc5 <- as_all[5]@score/15

    nas_all <- c(nas_bc1, nas_bc2, nas_bc3, nas_bc4, nas_bc5)
    bc_threshold <- 0.8

    m <- which(nas_all==max(nas_all))
    barcode_tie <- F
    barcode_2 <- NA
    barcode_passed_threshold <- F
    count <- 1
    for (i in m) {
      if (nas_all[i] != 0) {
        if (count == 1) {
          count = count + 1
          if (nas_all[i] > bc_threshold) barcode_passed_threshold = T
          if (i == 1) barcode <- 10
          else if (i == 2) barcode <- 30
          else if (i == 3) barcode <- 40
          else if (i == 4) barcode <- 60
          else if (i == 5) barcode <- 100
        } else if (count == 2) {
          count <- count + 1
          barcode_tie <- T
          if (i == 1) barcode_2 <- 10
          else if (i == 2) barcode_2 <- 30
          else if (i == 3) barcode_2 <- 40
          else if (i == 4) barcode_2 <- 60
          else if (i == 5) barcode_2 <- 100
        } else break
      } else {
        barcode <- 0
      }

    }
    rtl <- FALSE
    rts <- FALSE
    if (nchar(fastq[j]) > 1100 ) {
      rtl <- TRUE
    } else if (nchar(fastq[j]) < 750 ) {
      rts <- TRUE
    }
    #write read to the file it belongs to
    if(read_type=="reverse"){
      print("writing reverse")
      if(barcode_passed_threshold==TRUE & barcode==10){#Store TGCAGATCTCTTGCC alignments
        writeFastq(fqFile[j], "/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/1_aligned/gfp_corrected/reverse/strictDemux/revBar1.fastq", mode='a', compress=FALSE)
      }
      if(barcode_passed_threshold==TRUE & barcode==30){#Store CTCGAAGCATTGTAA alignments
        writeFastq(fqFile[j], "/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/1_aligned/gfp_corrected/reverse/strictDemux/revBar2.fastq", mode='a', compress=FALSE)
      }
      if(barcode_passed_threshold==TRUE & barcode==40){#Store AACGGTAGCCACCAA alignments
        writeFastq(fqFile[j], "/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/1_aligned/gfp_corrected/reverse/strictDemux/revBar3.fastq", mode='a', compress=FALSE)
      }
      if(barcode_passed_threshold==TRUE & barcode==60){#Store TGCACGAGATTGATG alignments
        writeFastq(fqFile[j], "/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/1_aligned/gfp_corrected/reverse/strictDemux/revBar4.fastq",mode='a', compress=FALSE)
      }
      if(barcode_passed_threshold==TRUE & barcode==100){#Store GACACATAGTCATGG alignments
        writeFastq(fqFile[j], "/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/1_aligned/gfp_corrected/reverse/strictDemux/revBar5.fastq", mode='a', compress=FALSE)
      }
    }
    else{
      next()
      writeFastq(fqFile[j], "/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/1_aligned/gfp_corrected/reverse/strictDemux/revFor.fastq", mode='a', compress=FALSE)
      print("got a rogue forward")
    }
  }

tryCatch({
  df <- demultiplex_reads(read_path=NA,
                          read_id_fast5_file=NA,
                          submat=submat,
                          type=type,
                          gapOpening=gapOpening,
                          gapExtension=gapExtension,
                          gfp=gfp,
                          rc_gfp=rc_gfp,
                          fp=fp,
                          rc_fp=rc_fp,
                          data=data,
                          multifast5 = 'F')
},
error=function(e){
  print(e)
})

