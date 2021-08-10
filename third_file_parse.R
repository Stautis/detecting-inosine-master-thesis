library(hdf5r)
library(ggplot2)
install.packages('e1071')
library(e1071)

#' Extract the read data from a fast5 file. Assumes that you have hdf5r librarr
#' already installed.
#'
#' This function can deal with both multifast5 files and single fast5 files.
#' It can handle files basecalled with standard or the flip flop model.
#'
#' @param file_path a character string [NA]. Path of the single fast5 file.
#' Use it if the file to be read is a single fast file. If the file to be
#' read is multifast5 file, then keep this parameter as NA. Also set
#' multifast5 flag to FALSE
#'
#' @param read_id_fast5_file a list [NA]. A list of 'read_id' and 'fast5_file'
#' path. Use this option when a read from a multifast5 file is to be read. In
#' such a case, you should set file_path to NA, and set multifast5 flag to TRUE.
#'
#' @param plot_debug a logical [FALSE]. Should data for plotting debug info in
#' the plots be computed
#'
#' @param basecalled_with a character string. Specify if the data is from
#''albacore' or 'guppy'
#'
#' @param basecall_group a character string. Name of the level
#' in the Fast5 file hierarchy from which to read the data e.g. "Basecall_1D_000"
#'
#' @param multifast5 a logical. Set it to TRUE if the file to be processed
#' is multifast5. Set it to FALSE if the file to be processed is a single fast5
#' file
#'
#' @param model a string. Set to 'flipflop' if the basecalling model is flipflop.
#' Set to 'standard' if the basecalling model is standard model.
#'
#' @param plotting_library a string
#'
#' @return a list
#'
extract_read_data <- function(file_path = NA,
                              read_id_fast5_file = NA,
                              plot_debug = FALSE,
                              basecalled_with = 'guppy',
                              basecall_group = 'Basecall_1D_000',
                              multifast5 = FALSE,
                              model = 'flipflop',
                              plotting_library = 'rbokeh') {
  if (!multifast5) {
    f5_obj <- hdf5r::H5File$new(file_path, mode='r')
    f5_tree <- f5_obj$ls(recursive=TRUE)
    f5_tree <- f5_tree$name
    # define all the paths
    raw_signal_path <- grep('.*Signal$', f5_tree, perl = TRUE, value = TRUE)
    if (sum(grepl('Raw/Reads/Read_[0-9]+$', f5_tree)) == 0) {
      read_id_path <-  grep('.*Raw$', f5_tree, perl = TRUE, value = TRUE)
    } else {
      read_id_path <- f5_tree[grepl('Raw/Reads/Read_[0-9]+$', f5_tree)]
    }
    # make fastq, Event/Move path
    event_data_fastq_path <- grep(paste0('.*', basecall_group, '/BaseCalled_template$'),
                                  f5_tree, perl = TRUE, value = TRUE)
    channel_path <- grep('channel_id$', f5_tree, perl=TRUE, value=TRUE)
    
    # make segmentation path based on the basecall group
    sp <- strsplit(basecall_group, split = '_')
    seg_group <- sp[[1]][3]
    segmentation_path <- grep(paste0('.*', seg_group, '/Summary/segmentation$'),
                              f5_tree, perl = TRUE, value = TRUE)
    # make basecalled_template path
    basecall_1d_template_path <- grep(paste0('.*', basecall_group, '/Summary/basecall_1d_template$'),
                                      f5_tree, perl = TRUE, value = TRUE)
  }
  else {
    f5_obj <- hdf5r::H5File$new(read_id_fast5_file$fast5_file, mode='r')
    full_read_id <- read_id_fast5_file$read_id
    # define all the paths
    raw_signal_path <- paste('/', full_read_id, '/Raw/Signal', sep='')
    event_data_fastq_path <- paste0('/', full_read_id, '/Analyses/', basecall_group, '/BaseCalled_template/')
    read_id_path <- paste('/', full_read_id, '/Raw', sep='')
    
    # make segmentation path based on the basecall group
    sp <- strsplit(basecall_group, split='_')
    seg_group <- sp[[1]][3]
    segmentation_path <- paste0('/', full_read_id, '/Analyses/Segmentation_', seg_group, '/Summary/segmentation')
    basecall_1d_template_path <- paste0('/', full_read_id, '/Analyses/', basecall_group, '/Summary/basecall_1d_template')
  }
  
  # get the data
  raw_data <- f5_obj[[raw_signal_path]]$read()
  read_id <- f5_obj[[read_id_path]]$attr_open('read_id')$read()
  fastq <- f5_obj[[event_data_fastq_path]]$open('Fastq')$read()
  start <- f5_obj[[segmentation_path]]$attr_open('first_sample_template')$read()
  digi <- f5_obj[[channel_path]]$attr_open('digitisation')$read()
  offset <- f5_obj[[channel_path]]$attr_open('offset')$read()
  range <- f5_obj[[channel_path]]$attr_open('range')$read()
  called_events <- f5_obj[[basecall_1d_template_path]]$attr_open('called_events')$read()
  
  # compute called_event if it has been set to W by ONT in the FAST5 file (Wierd! I know)
  compute_called_events <- ifelse(called_events == 'W', TRUE, FALSE)
  # get the event data, if present -- or make it, if not present
  make_event_data = FALSE
  if ('Events' %in% names(f5_obj[[event_data_fastq_path]])){
    event_data <- f5_obj[[event_data_fastq_path]]$open('Events')$read()
    if (!('length' %in% colnames(event_data))) {
      stride <- f5_obj[[basecall_1d_template_path]]$attr_open('block_stride')$read()
    } else {
      stride <- event_data$length[1]
    }
    
    # Albacore latest version does not output start column;
    # this is a dirty fix for it
    if (!("start" %in% colnames(event_data))) {
      make_event_data = TRUE
      move <- event_data$move
    }
    
  } else {
    move <- f5_obj[[event_data_fastq_path]]$open('Move')$read()
    stride <- f5_obj[[basecall_1d_template_path]]$attr_open('block_stride')$read()
    make_event_data = TRUE
  }
  
  # extract just the fastq removing quality scores
  fastq <- strsplit(fastq, split = "\n")
  fastq <- fastq[[1]][2]
  
  
  # if event_data wasn't present, make it now
  if (make_event_data) {
    # The line below is there to remove R CMD CHECK
    # "no visible binding for global variable" error
    move_cumsum <- fastq_bases <- NULL
    
    start_col <- seq(from = start,
                     to = start + stride * (length(move) - 1),
                     by = stride)
    
    # in case of RNA, the event data fastQ should be the reverse
    # of actual FASTQ
    # find U in the sequence to determine if it is RNA
    # find if data is dna or rna
    if (!multifast5) {
      context_tags_path <- f5_tree[grepl('.*context_tags$', f5_tree)]
    } else {
      context_tags_path <- paste(first_read_name,
                                 '/context_tags',
                                 sep = '')
    }
    sequencing_kit <- f5_obj[[context_tags_path]]$attr_open('sequencing_kit')$read()
    if (grepl('rna', sequencing_kit)) {
      experiment_type <- 'rna'
    } else {
      experiment_type <- 'dna'
    }
    
    if (experiment_type == 'dna') {
      event_data <- data.frame(move = move,
                               start = start_col,
                               move_cumsum = cumsum(move),
                               fastq_bases = fastq,
                               stringsAsFactors = FALSE)
    } else {
      fastq_rev <- Biostrings::reverse(fastq)
      event_data <- data.frame(move = move,
                               start = start_col,
                               move_cumsum = cumsum(move),
                               fastq_bases = fastq_rev,
                               stringsAsFactors = FALSE)
    }
    
    # The line below is there to remove R CMD CHECK
    # "no visible binding for global variable" error
    model_state <- NULL
    event_data <- dplyr::mutate(event_data,
                                model_state = substr(fastq_bases,
                                                     start=move_cumsum,
                                                     stop=move_cumsum))
    event_data <- dplyr::select(event_data, model_state, start, move)
  }
  
  # compute the number of event if not present already
  if (compute_called_events) {
    called_events <- nrow(event_data)
  }
  
  # make a vector of moves interpolated for every sample i.e., make a sample-wise or per-sample vector of moves
  if (plot_debug & plotting_library == 'ggplot2') {
    if (start != 0) {
      moves_sample_wise_vector <- c(rep(NA, start-1),
                                    rep(event_data$move*0.25+1.5, each=stride),
                                    rep(NA, length(raw_data) - start - stride*called_events + 1))
    } else {
      moves_sample_wise_vector <- c(rep(event_data$move*0.25+1.5, each=stride),
                                    rep(NA, length(raw_data) - start - stride*called_events))
    }
  } else {
    moves_sample_wise_vector <- rep(NA, length(raw_data))
  }
  
  # compute event length vector
  if (model == 'flipflop') {
    # create event length data for tail normalization
    event_length_vector <- rep(NA, called_events)
    count <- 0
    for (i in seq(from=called_events, to=1, by=-1)) {
      if (event_data$move[i] == 1) {
        event_length_vector[i] <- count + 1
        count <- 0
      } else {
        count <- count + 1
      }
    }
    # multiply moves by length of the event (15 for RNA, 5 for DNA)
    event_length_vector <- event_length_vector * stride
    event_data <- cbind(event_data, event_length_vector)
    # remove NAs
    event_length_vector <- event_length_vector[!is.na(event_length_vector)]
    # Normalizer for flip-flop based data
    samples_per_nt <- mean(event_length_vector[event_length_vector <= stats::quantile(event_length_vector, 0.95)])
    # add the start column to the event table for legacy purposes
    start_col <-seq(from=start, to=(start + (nrow(event_data)-1)*stride), by=stride)
    event_data <- dplyr::mutate(event_data, start=start_col)
  } else if (model == 'standard') {
    # create event length data for tail normalization
    l <- stride
    event_length_vector_1 <-  rep(NA, called_events)
    event_length_vector_2 <-  rep(NA, called_events)
    index <- called_events
    length_count <- 1
    divide_by <- 1
    
    # handle the first row of the event data
    if (event_data$move[1]==1) {
      event_length_vector_1[1] <- 1
    } else {
      event_length_vector_1[index] <- 0.5
      event_length_vector_2[index] <- 0.5
    }
    for(i in (called_events):2) {
      if (event_data$move[i]==2 & event_data$move[i-1]==0) {
        # record the index
        index <- i
        length_count <- 1
        divide_by <- 2
      } else if (event_data$move[i]==1 & event_data$move[i-1]==0) {
        index <- i
        length_count <- 1
        divide_by <- 1
      } else if (event_data$move[i]==0 & (event_data$move[i-1]==1 | event_data$move[i-1]==2)) {
        length_count <- length_count + 1
        # put the previous record
        if (divide_by == 1) {
          event_length_vector_1[index] = length_count
        } else {
          event_length_vector_1[index] = length_count/2
          event_length_vector_2[index] = length_count/2
        }
      } else if ((event_data$move[i]==2 & event_data$move[i-1]==1) | (event_data$move[i]==2 & event_data$move[i-1]==2)) {
        event_length_vector_1[i] <- 0.5
        event_length_vector_2[i] <- 0.5
      } else if ((event_data$move[i]==1 & event_data$move[i-1]==1) | (event_data$move[i]==1 & event_data$move[i-1]==2))  {
        event_length_vector_1[i] <- 1
      } else if  (event_data$move[i]==0 & event_data$move[i-1]==0) {
        length_count <- length_count + 1
      }
    }
    # multiply moves by length of the event (15 for RNA, 5 for DNA)
    event_length_vector_1 <- event_length_vector_1 * l
    event_length_vector_2 <- event_length_vector_2 * l
    #event_data <- dplyr::select(event_data, start, length, move, model_state)
    event_data <- cbind(event_data, event_length_vector_1, event_length_vector_2)
    # combine the two vectors
    event_length_vector <- c(event_length_vector_1, event_length_vector_2)
    # reomve NAs
    event_length_vector <- event_length_vector[!is.na(event_length_vector)]
    # compute geometric mean of modfied Albacore events table to get the normalizer
    samples_per_nt <- psych::geometric.mean(event_length_vector)
  }
  f5_obj$close_all()
  
  list(read_id = read_id,
       raw_data = raw_data,
       event_data = event_data,
       moves_sample_wise_vector = moves_sample_wise_vector,
       fastq = fastq,
       start = start,
       stride = stride,
       digi = digi,
       range = range,
       offset = offset,
       samples_per_nt = samples_per_nt)
}

#Read in list of relevant files
files_3 <- read.csv("/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/3_rebasecalled/file_paths.csv", sep="\n", header=FALSE, stringsAsFactors = F)
files_3 <- files_3[[1]]

#prepare pairwise alignment strings
preI <- Biostrings::DNAString('TGAATTCGAATCGATGGGATCCCTAAGA')
rc_preI <- Biostrings::reverseComplement(preI)
postI <- Biostrings::DNAString('TGCGTCCTGTGCCGA')
rc_postI <- Biostrings::reverseComplement(postI)
preT <- Biostrings::DNAString('TGCGTCCTGTGCCGA')
postT <- Biostrings::DNAString('AGCCTGCT')
preG <- Biostrings::DNAString('TCCACTCCTC')
postG <- Biostrings::DNAString('TCAGTGGTG')
preC <- Biostrings::DNAString('TCAGTGGTG')
postC <- Biostrings::DNAString('AGTCAATGCTCGT')
preA <- Biostrings::DNAString('AGCCTGCTG')
postA <- Biostrings::DNAString('TCCACTCCT')
homo_T <- Biostrings::DNAString('TTTTTT')
homo_A <- Biostrings::DNAString('AAAAAA')
homo_G <- Biostrings::DNAString('GGGGGG')
homo_C <- Biostrings::DNAString('CCCCCC')

#define substitution matrix for alignments
match <- 1
mismatch <- -1
type <- 'local'
gapOpening <- 0
gapExtension <- -1
submat <- Biostrings::nucleotideSubstitutionMatrix(match=match,
                                                   mismatch=mismatch,
                                                   baseOnly=TRUE)

#initialize counters, lists and dataframes to save relevant data
third_detectionCounter <- 0
third_counter <- 0
third_avgsum <- 0.0
third_rows <- c()
row_frame <- data.frame(value=double(), row=integer(),eventlength=integer(),read_id=character())
third_raw_list <- c()
third_box_list <- c()
for(row in 1:length(files_3)){
  third_counter=third_counter+1
  print(paste("processing row",row))
  
  #Initialize dataframe and normalize signals
  df <- extract_read_data(file_path=files_3[row])
  read_id <- df$read_id
  raw_signal <- df$raw_data
  event_data <- df$event_data
  seq <- df$fastq
  range <- df$range
  digi <- df$digi
  offset <- df$offset
  scale <- range/digi
  normalized_raw_signal <- scale * (raw_signal + offset)
  
  #Align to reigon flanking inosine homopolymer to the left
  as_preI <- Biostrings::pairwiseAlignment(pattern=toString(preI),
                                           subject=seq,
                                           substitutionMatrix=submat,
                                           type=type,
                                           scoreOnly=FALSE,
                                           gapOpening=gapOpening,
                                           gapExtension=gapExtension)
  
  #Align to the reverse complement of the same region
  as_rc_preI <- Biostrings::pairwiseAlignment(pattern=toString(rc_preI),
                                              subject=seq,
                                              substitutionMatrix=submat,
                                              type=type,
                                              scoreOnly=FALSE,
                                              gapOpening=gapOpening,
                                              gapExtension=gapExtension)
  
  #Normalize alignment scores
  nas_preI <- as_preI@score/preI@length
  nas_rc_preI <- as_rc_preI@score/rc_preI@length
  
  #assess whether alignment threshold has been surpassed
  if (nas_preI > 0.7 | nas_rc_preI > 0.7){
    #write.csv(event_data, "/export/valenfs/data/processed_data/MinION/inosine_project/20181012_max_DNA_Inosine/direct_signal_analysis/output/frame.csv")
    
    #Align to region flanking inosine homopolymer to the right
    as_postI <- Biostrings::pairwiseAlignment(pattern=toString(postI),
                                              subject=seq,
                                              substitutionMatrix=submat,
                                              type=type,
                                              scoreOnly=FALSE,
                                              gapOpening=gapOpening,
                                              gapExtension=gapExtension)
    
    as_rc_postI <- Biostrings::pairwiseAlignment(pattern=toString(rc_postI),
                                                 subject=seq,
                                                 substitutionMatrix=submat,
                                                 type=type,
                                                 scoreOnly=FALSE,
                                                 gapOpening=gapOpening,
                                                 gapExtension=gapExtension)
    
    nas_postI <- as_postI@score/postI@length
    nas_rc_postI <- as_rc_postI@score/rc_postI@length
    
    if (nas_postI > 0.6){
      substr(seq, as_preI@subject@range@start+as_preI@subject@range@width, as_postI@subject@range@start-1) #finds inosine stretch

      #navigate table of raw signals to the signals corresponding to the inosine homopolymer
      if(abs((as_postI@subject@range@start-1) - (as_preI@subject@range@start+as_preI@subject@range@width)) < 20){
        countI <- 0
        string <- ''
        starts <- ''
        lengths <- ''
        first <- 0
        last <- 0
        lastStart <- 0
        eventlength <- 0
        f <- function(x){
          if(x[3]==1){
            countI<<-countI + 1
            print(countI)
          }
          if(countI>=as_preI@subject@range@start+as_preI@subject@range@width & countI<as_postI@subject@range@start){
            if (string==''){
              first<<-as.numeric(x[2])
            }
            string<<-paste(string, x[1], sep="")
            if (x[3]==1){
              last<<-as.numeric(x[4])
              lastStart<<-as.numeric(x[2])
            }
            if (!is.na(x[4])){
              eventlength<<-eventlength+as.numeric(x[4])
            }
            starts<<-paste(starts, x[2])
            lengths<<-paste(lengths, x[4])
          }
        }
        apply(event_data, 1, f)
        if(first==lastStart){
          print("nexting")
          next
        }
	#calculate descriptives
        signalTotal <- 0.0
        count <- 0
        signals <- c()
        for (i in first:(lastStart+last)){
          signals <- c(signals, normalized_raw_signal[i])
          signalTotal <- signalTotal + normalized_raw_signal[i]
          count <- count+1
        }
	#save all data for later calculations
        third_box_list <- c(third_box_list, signals[round(length(signals)*0.25):round(length(signals)*0.75)])
        third_raw_list <- c(third_raw_list, signals)
        row_frame_temp <- data.frame(value=signals[round(length(signals)*0.25):round(length(signals)*0.75)], row_nr=row,length=(last+lastStart)-first,read_id=read_id)
        row_frame <- rbind(row_frame,row_frame_temp)
        avg_signal <- signalTotal/count
        third_avgsum <- third_avgsum+avg_signal
        print(third_avgsum)
      }
      #Distance between right and left flanking regions is too great
      else{
        print("Dismissing candidate! too far apart")
        next
      }
    }
    #reverse complement of inosine homopolymer does not hold relevant information so we skip
    else if (nas_rc_postI > 0.7){
      print("skipping reverse")
      next()
      #third_rows <- c(third_rows, row)
      #print("adding and skipping reverse")
      #next
      if (abs((as_rc_preI@subject@range@start-1) - (as_rc_postI@subject@range@start+as_rc_postI@subject@range@width)) < 20){
        countI <- 0
        string <- ''
        starts <- ''
        lengths <- ''
        first <- 0
        last <- 0
        lastStart <- 0
        f <- function(x){
          if(x[3]==1){
            countI<<-countI + 1
            print(countI)
          }
          if(countI>=as_rc_postI@subject@range@start+as_rc_postI@subject@range@width & countI<as_rc_preI@subject@range@start){
            if (string==''){
              first<<-as.numeric(x[2])
            }
            string<<-paste(string, x[1], sep="")
            if (x[3]==1){
              last<<-as.numeric(x[4])
              lastStart<<-as.numeric(x[2])
            }
            starts<<-paste(starts, x[2])
            lengths<<-paste(lengths, x[4])
          }
        }
        apply(event_data, 1, f)
        signalTotal <- 0.0
        count <- 0
        signals <- c()
        for (i in first:(lastStart+last)){
          signals <- c(signals, normalized_raw_signal[i])
          signalTotal <- signalTotal + normalized_raw_signal[i]
          count <- count+1
        }
        third_box_list <- c(third_box_list, signals[round(length(signals)*0.25):round(length(signals)*0.75)])
        third_raw_list <- c(third_raw_list, signals)
        avg_signal <- signalTotal/count
        third_avgsum <- third_avgsum+avg_signal
        print(third_avgsum)
      }
      else{
        print("Dismissing reverse candidate! too far apart")
        next
      }
    }
    else{
      print("Dismissing candidate! no reverse post section found")
      next
    }
    #sanity check for non-numerical results
    if(is.na(third_avgsum)){
      stop("invalid score")
    }
    third_detectionCounter <- third_detectionCounter + 1
    if (third_detectionCounter > 999){
      print("1000 Candidates detected!")
    }
    else{
      third_rows <- c(third_rows, row)
      print(paste("Candidate detected!", third_detectionCounter))
      next
    }
  }
  else{
    print(paste("got scores", nas_preI,"and",nas_rc_preI))
  }
}

third_detCount2 <- 0
third_avgsum_2 <- 0.0
third_raw_list_2 <- c()
third_box_list_2 <- c()
row_frame_2 <- data.frame(value=double(), row=integer(),eventlength=integer(),read_id=character())
for(j in 1:length(files_3)){
  df2 <- extract_read_data(file_path=files_3[j])
  read_id <- df$read_id
  raw_signal2 <- df2$raw_data
  event_data2 <- df2$event_data
  seq2 <- df2$fastq
  range2 <- df2$range
  digi2 <- df2$digi
  offset2 <- df2$offset
  scale2 <- range2/digi2
  normalized_raw_signal2 <- (raw_signal2 + offset2) * scale2 
  
  
  as_preT <- Biostrings::pairwiseAlignment(pattern=toString(preT),
                                           subject=seq2,
                                           substitutionMatrix=submat,
                                           type=type,
                                           scoreOnly=FALSE,
                                           gapOpening=gapOpening,
                                           gapExtension=gapExtension)
  
  as_postT <- Biostrings::pairwiseAlignment(pattern=toString(postT),
                                            subject=substr(seq2,as_preT@subject@range@start,nchar(seq2)),
                                            substitutionMatrix=submat,
                                            type=type,
                                            scoreOnly=FALSE,
                                            gapOpening=gapOpening,
                                            gapExtension=gapExtension)
  
  nas_preT <- as_preT@score/preT@length
  nas_postT <- as_postT@score/postT@length
  
  if(nas_preT > 0.8 & nas_postT > 0.7){
    if((as_preT@subject@range@start-as_preT@subject@range@start+as_postT@subject@range@start)<(as_preT@subject@range@width+3)){
      print("too close")
      next()
    }
    print("Found T homopolymer flanking regions")
    print(seq2)
    countI <- 0
    string <- ''
    starts <- ''
    lengths <- ''
    first <- 0
    last <- 0
    lastStart <- 0
    raw_list <- c()
    f <- function(x){
      if(x[3]==1){
        countI<<-countI + 1
        print(countI)
      }
      if(countI>=as_preT@subject@range@start+as_preT@subject@range@width & countI<as_preT@subject@range@start+as_postT@subject@range@start){  
        if (string==''){
          first<<-as.numeric(x[2])
        }
        string<<-paste(string, x[1], sep="")
        if (x[3]==1){
          last<<-as.numeric(x[4])
          lastStart<<-as.numeric(x[2])
        }
        starts<<-paste(starts, x[2])
        lengths<<-paste(lengths, x[4])
      }
    }
    apply(event_data2, 1, f)
    signalTotal <- 0.0
    count <- 0
    signals_2 <- c()
    for (i in first:(lastStart+last)){
      signals_2 <- c(signals_2, normalized_raw_signal2[i])
      signalTotal <- signalTotal + normalized_raw_signal2[i]
      count <- count+1
    }
    third_box_list_2 <- c(third_box_list_2, signals_2[round(length(signals_2)*0.25):round(length(signals_2)*0.75)])
    third_raw_list_2 <- c(third_raw_list_2, signals_2)
    row_frame_temp <- data.frame(value=signals[round(length(signals)*0.25):round(length(signals)*0.75)], row_nr=j,length=(last+lastStart)-first,read_id=read_id)
    row_frame_2 <- rbind(row_frame_2,row_frame_temp)
    avg_signal_2 <- signalTotal/count
    third_avgsum_2 <- third_avgsum_2+avg_signal_2
    third_detCount2 <- third_detCount2 +1
  }
  if(length(third_avgsum_2)==0 | is.na(third_avgsum_2)){
    stop()
  }
  print(paste("Got", third_avgsum_2))
  print(paste("done", third_detCount2))
}

third_detCount3 <- 0
third_avgsum_3 <- 0.0
third_raw_list_3 <- c()
third_box_list_3 <- c()
row_frame_3 <- data.frame(value=double(), row=integer(),eventlength=integer(),read_id=character())
for(j in 1:length(files_3)){
  df3 <- extract_read_data(file_path=files_3[j])
  read_id <- df$read_id
  raw_signal3 <- df3$raw_data
  event_data3 <- df3$event_data
  seq3 <- df3$fastq
  range3 <- df3$range
  digi3 <- df3$digi
  offset3 <- df3$offset
  scale3 <- range3/digi3
  normalized_raw_signal3 <- (raw_signal3 + offset3) * scale3 
  
  as_preI <- Biostrings::pairwiseAlignment(pattern=toString(preI),
                                           subject=seq3,
                                           substitutionMatrix=submat,
                                           type=type,
                                           scoreOnly=FALSE,
                                           gapOpening=gapOpening,
                                           gapExtension=gapExtension)
  
  as_preG <- Biostrings::pairwiseAlignment(pattern=toString(preG),
                                            subject=substr(seq3,as_preI@subject@range@start,nchar(seq3)),
                                            substitutionMatrix=submat,
                                            type=type,
                                            scoreOnly=FALSE,
                                            gapOpening=gapOpening,
                                            gapExtension=gapExtension)
  
  as_postG <- Biostrings::pairwiseAlignment(pattern=toString(postG),
                                           subject=substr(seq3,as_preI@subject@range@start,nchar(seq3)),
                                           substitutionMatrix=submat,
                                           type=type,
                                           scoreOnly=FALSE,
                                           gapOpening=gapOpening,
                                           gapExtension=gapExtension)
  
  
  nas_preI <- as_preI@score/preI@length
  nas_preG <- as_preG@score/preG@length
  nas_postG <- as_postG@score/postG@length
  
  if(as_homoG@subject@range@start==1){
    next()
  }
  
  if(nas_preG >0.8 & nas_postG >0.8){
    if(abs(as_postG@subject@range@start-(as_preG@subject@range@start+as_preG@subject@range@width))>15 | abs(as_postG@subject@range@start-(as_preG@subject@range@start+as_preG@subject@range@width))<4){
      print("too close")
      next()
    }
    print("Found G homopolymer")
    print(seq3)
    countI <- 0
    string <- ''
    starts <- ''
    lengths <- ''
    first <- 0
    last <- 0
    lastStart <- 0
    raw_list <- c()
    f <- function(x){
      if(x[3]==1){
        countI<<-countI + 1
        print(countI)
      }
      if(countI>=(as_preI@subject@range@start+as_preG@subject@range@start+as_preG@subject@range@width-1) & countI<(as_preI@subject@range@start+as_postG@subject@range@start-1)){
        if (string==''){
          first<<-as.numeric(x[2])
        }
        string<<-paste(string, x[1], sep="")
        if (x[3]==1){
          last<<-as.numeric(x[4])
          lastStart<<-as.numeric(x[2])
        }
        starts<<-paste(starts, x[2])
        lengths<<-paste(lengths, x[4])
      }
    }
    apply(event_data3, 1, f)
    signalTotal <- 0.0
    count <- 0
    signals_3 <- c()
    for (i in first:(lastStart+last)){
      signals_3 <- c(signals_3, normalized_raw_signal3[i])
      signalTotal <- signalTotal + normalized_raw_signal3[i]
      count <- count+1
    }
    third_box_list_3 <- c(third_box_list_3, signals_3[round(length(signals_3)*0.25):round(length(signals_3)*0.75)])
    third_raw_list_3 <- c(third_raw_list_3, signals_3)
    row_frame_temp <- data.frame(value=signals[round(length(signals)*0.25):round(length(signals)*0.75)], row_nr=j,length=(last+lastStart)-first,read_id=read_id)
    row_frame_3 <- rbind(row_frame_3,row_frame_temp)
    avg_signal_3 <- signalTotal/count
    third_avgsum_3 <- third_avgsum_3+avg_signal_3
    third_detCount3 <- third_detCount3 +1
  }
  if(length(third_avgsum_3)==0 | is.na(third_avgsum_3)){
    stop()
  }
  print(paste("Got", third_avgsum_3))
  print(paste("done", third_detCount3))
}

third_detCount4 <- 0
third_avgsum_4 <- 0.0
third_raw_list_4 <- c()
third_box_list_4 <- c()
row_frame_4 <- data.frame(value=double(), row=integer(),eventlength=integer(),read_id=character())
for(j in 1:length(files_3)){
  df4 <- extract_read_data(file_path=files_3[j])
  read_id <- df$read_id
  raw_signal4 <- df4$raw_data
  event_data4 <- df4$event_data
  seq4 <- df4$fastq
  range4 <- df4$range
  digi4 <- df4$digi
  offset4 <- df4$offset
  scale4 <- range4/digi4
  normalized_raw_signal4 <- (raw_signal4 + offset4) * scale4 
  
  as_preI <- Biostrings::pairwiseAlignment(pattern=toString(preI),
                                           subject=seq4,
                                           substitutionMatrix=submat,
                                           type=type,
                                           scoreOnly=FALSE,
                                           gapOpening=gapOpening,
                                           gapExtension=gapExtension)
  
  as_preC <- Biostrings::pairwiseAlignment(pattern=toString(preC),
                                           subject=substr(seq4,as_preI@subject@range@start,nchar(seq4)),
                                           substitutionMatrix=submat,
                                           type=type,
                                           scoreOnly=FALSE,
                                           gapOpening=gapOpening,
                                           gapExtension=gapExtension)
  
  as_postC <- Biostrings::pairwiseAlignment(pattern=toString(postC),
                                            subject=substr(seq4,as_preI@subject@range@start,nchar(seq4)),
                                            substitutionMatrix=submat,
                                            type=type,
                                            scoreOnly=FALSE,
                                            gapOpening=gapOpening,
                                            gapExtension=gapExtension)
  
  
  nas_preI <- as_preI@score/preI@length
  nas_preC <- as_preC@score/preC@length
  nas_postC <- as_postC@score/postC@length
  
  if(nas_preC >0.8 & nas_postC >0.7){
    if(abs(as_postC@subject@range@start-(as_preC@subject@range@start+as_preC@subject@range@width))>15 | abs(as_postC@subject@range@start-(as_preC@subject@range@start+as_preC@subject@range@width))<4){
      print("too close")
      next()
    }
    print("Found C homopolymer")
    print(seq4)
    countI <- 0
    string <- ''
    starts <- ''
    lengths <- ''
    first <- 0
    last <- 0
    lastStart <- 0
    raw_list <- c()
    f <- function(x){
      if(x[3]==1){
        countI<<-countI + 1
        print(countI)
      }
      if(countI>=(as_preI@subject@range@start+as_preC@subject@range@start+as_preC@subject@range@width-1) & countI<(as_preI@subject@range@start+as_postC@subject@range@start-1)){
        if (string==''){
          first<<-as.numeric(x[2])
        }
        string<<-paste(string, x[1], sep="")
        if (x[3]==1){
          last<<-as.numeric(x[4])
          lastStart<<-as.numeric(x[2])
        }
        starts<<-paste(starts, x[2])
        lengths<<-paste(lengths, x[4])
      }
    }
    apply(event_data4, 1, f)
    signalTotal <- 0.0
    count <- 0
    signals_4 <- c()
    for (i in first:(lastStart+last)){
      signals_4 <- c(signals_4, normalized_raw_signal4[i])
      signalTotal <- signalTotal + normalized_raw_signal4[i]
      count <- count+1
    }
    third_box_list_4 <- c(third_box_list_4, signals_4[round(length(signals_4)*0.25):round(length(signals_4)*0.75)])
    third_raw_list_4 <- c(third_raw_list_4, signals_4)
    row_frame_temp <- data.frame(value=signals[round(length(signals)*0.25):round(length(signals)*0.75)], row_nr=j,length=(last+lastStart)-first,read_id=read_id)
    row_frame_4 <- rbind(row_frame_4,row_frame_temp)
    avg_signal_4 <- signalTotal/count
    third_avgsum_4 <- third_avgsum_4+avg_signal_4
    third_detCount4 <- third_detCount4 +1
  }
  if(length(third_avgsum_4)==0 | is.na(third_avgsum_4)){
    stop()
  }
  print(paste("Got", third_avgsum_4))
  print(paste("done", third_detCount4))
}

third_detCount5 <- 0
third_avgsum_5 <- 0.0
third_raw_list_5 <- c()
third_box_list_5 <- c()
row_frame_5 <- data.frame(value=double(), row=integer(),eventlength=integer(),read_id=character())
miscount <- 0
for(row in 1:length(files_3)){
  print(paste("processing row",row))
  
  df5 <- extract_read_data(file_path=files_3[row])
  read_id <- df$read_id
  raw_signal5 <- df5$raw_data
  event_data5 <- df5$event_data
  seq5 <- df5$fastq
  range5 <- df5$range
  digi5 <- df5$digi
  offset5 <- df5$offset
  scale5 <- range5/digi5
  normalized_raw_signal5 <- scale5 * (raw_signal5 + offset5)
  
  
  as_preI <- Biostrings::pairwiseAlignment(pattern=toString(preI),
                                           subject=seq5,
                                           substitutionMatrix=submat,
                                           type=type,
                                           scoreOnly=FALSE,
                                           gapOpening=gapOpening,
                                           gapExtension=gapExtension)
  
  as_preA <- Biostrings::pairwiseAlignment(pattern=toString(preA),
                                           subject=substr(seq5,as_preI@subject@range@start,nchar(seq5)),
                                           substitutionMatrix=submat,
                                           type=type,
                                           scoreOnly=FALSE,
                                           gapOpening=gapOpening,
                                           gapExtension=gapExtension)
  
  as_postA <- Biostrings::pairwiseAlignment(pattern=toString(postA),
                                            subject=substr(seq5,as_preI@subject@range@start,nchar(seq5)),
                                            substitutionMatrix=submat,
                                            type=type,
                                            scoreOnly=FALSE,
                                            gapOpening=gapOpening,
                                            gapExtension=gapExtension)
  
  nas_preI <- as_preI@score/preI@length
  nas_preA <- as_preA@score/preA@length
  nas_postA <- as_postA@score/postA@length
  
  nas_preI <- as_preI@score/preI@length
  nas_rc_preI <- as_rc_preI@score/rc_preI@length
  
  if(nas_preA >0.8 & nas_postA >0.8){
    if(abs(as_postA@subject@range@start-(as_preA@subject@range@start+as_preA@subject@range@width))>15 | abs(as_postA@subject@range@start-(as_preA@subject@range@start+as_preA@subject@range@width))<4){
      print("too close")
      next()
    }
    print("Found A homopolymer")
    print(seq5)
    print(as_homoA@subject@range@start)
    if(as_homoA@subject@range@start > as_preI@subject@range@start & as_homoA@subject@range@start < as_postI@subject@range@start){
      print("skipping inosine homopolymer mistaken for A's")
      miscount <- miscount+1
      next()
    }
    countI <- 0
    string <- ''
    starts <- ''
    lengths <- ''
    first <- 0
    last <- 0
    lastStart <- 0
    raw_list <- c()
    f <- function(x){
      if(x[3]==1){
        countI<<-countI + 1
        print(countI)
      }
      if(countI>=(as_preI@subject@range@start+as_preA@subject@range@start+as_preA@subject@range@width-1) & countI<(as_preI@subject@range@start+as_postA@subject@range@start-1)){
        if (string==''){
          first<<-as.numeric(x[2])
        }
        string<<-paste(string, x[1], sep="")
        if (x[3]==1){
          last<<-as.numeric(x[4])
          lastStart<<-as.numeric(x[2])
        }
        starts<<-paste(starts, x[2])
        lengths<<-paste(lengths, x[4])
      }
    }
    apply(event_data5, 1, f)
    signalTotal <- 0.0
    count <- 0
    signals_5 <- c()
    for (i in first:(lastStart+last)){
      signals_5 <- c(signals_5, normalized_raw_signal5[i])
      signalTotal <- signalTotal + normalized_raw_signal5[i]
      count <- count+1
    }
    third_box_list_5 <- c(third_box_list_5, signals_5[round(length(signals_5)*0.25):round(length(signals_5)*0.75)])
    third_raw_list_5 <- c(third_raw_list_5, signals_5)
    row_frame_temp <- data.frame(value=signals[round(length(signals)*0.25):round(length(signals)*0.75)], row_nr=row,length=(last+lastStart)-first,read_id=read_id)
    row_frame_5 <- rbind(row_frame_5,row_frame_temp)
    avg_signal_5 <- signalTotal/count
    third_avgsum_5 <- third_avgsum_5+avg_signal_5
    third_detCount5 <- third_detCount5 +1
  }
  else{
    print("didn't find one")
    next()
    print(paste("got scores", nas_preI,"and",nas_rc_preI))
  }
  if(length(third_avgsum_5)==0 | is.na(third_avgsum_5)){
    stop()
  }
  print(paste("Got", third_avgsum_5))
  print(paste("done", third_detCount5))
}

#Make into one dataframe of box values
box_ino <- data.frame(value = row_frame, variable = "I")
box_thy <- data.frame(value = row_frame_2, variable = "T")
box_gua <- data.frame(value = row_frame_3, variable = "G")
box_cyt <- data.frame(value = row_frame_4, variable = "C")
box_ade <- data.frame(value = row_frame_5, variable = "A")
all_dat <- rbind(box_ino, box_thy, box_gua, box_cyt, box_ade)

#Write out data
write.csv(all_dat,"expanded_current_row_variable_dataset.csv",row.names = FALSE)

#make reference boxplot
ont_ref <- data.frame(value = c(90.679010,73.670019,98.890775,86.486336),sd=c(1.513907,1.512557,1.401820,1.517846), variable=c("T","G","C","A"))
#make violinplot
ggplot(ont_ref, aes(x=as.factor(variable))) + geom_boxplot(aes(lower=value-sd, upper=value+sd, middle=value, ymin=value-3*sd, ymax=value+3*sd),stat="identity")

#overlay the boxplot over the violinplot
v <- ggplot(all_dat_pos, aes(x=variable, y=value, yScale="log2")) + geom_violin(scale="width") + ylim(1,150) + stat_summary(fun.y="sample",geom="point", aes(x="T", y=90.67)) + stat_summary(fun.y="sample",geom="point", aes(x="G", y=73.67)) + stat_summary(fun.y="sample",geom="point", aes(x="C", y=98.89)) +  stat_summary(fun.y="sample",geom="point", aes(x="A", y=86.49))
  

