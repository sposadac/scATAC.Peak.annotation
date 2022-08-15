create_data_chunks <- function(data, chunk_size) {
  # internal function
  # data: matrix
  # chunk_size: integer
  data_list <- list()
  # Handle case: nrow(data) < chunk_size
  data_iterat <- nrow(data) %/% chunk_size
  if (data_iterat > 0) {
    for (n in 1:data_iterat) {
      start <- n * chunk_size - (chunk_size - 1)
      end <- min(n * chunk_size, nrow(data))
      data_list[[n]] <- data[start:end, ]
    }
  }
  if ((nrow(data) %% chunk_size) != 0) {
    start <- length(data_list) * chunk_size + 1
    data_list[[length(data_list) + 1]] <- data[start:nrow(data), ]
  }
  return(data_list)
}

stardardize_chr_names <- function(peaks, annotation) {
  # internal function: prefix chromosome names with "chr" (UCSC convention)
  # peaks: matrix
  # annotation: data.frame
  chrom_peaks <- unique(peaks[, 1])
  chromosomes <- unique(as.character(annotation$seqnames))
  if (sum(chromosomes %in% chrom_peaks) == 0) {
    # check if chromosomes in annotations are only numeric + X, Y, MT
    aux <- grepl("^[0-9XYMTxymt]+$", chromosomes, perl=TRUE) 
    if (all(aux)) {
      annotation$seqnames <- paste0("chr", annotation$seqnames)
    } else {
      peaks[, 1] <- paste0("chr", peaks[, 1])
    }
    # NOTE: Ensembl uses "MT" for the mitochondrial chromosome, use Refseq's "chrM"
    annotation$seqnames <- sub("chrMT", "chrM", annotation$seqnames)
    peaks[, 1] <- sub("chrMT", "chrM", peaks[, 1])
    chromosomes <- unique(as.character(annotation$seqnames))
    chrom_peaks <- unique(peaks[, 1])
    if (!all(chrom_peaks %in% chromosomes))
      stop("not all peak chromosome names are found in 'annotation'")
  }
  return(list("peaks"=peaks, "annotation"=annotation,
              "chromosomes"=chromosomes))
}

#' @importFrom future availableCores plan multisession
#' @import future.apply
#' @import Matrix
#' @import Matrix.utils
#' @import GenomicRanges
#
#'
#' @title  get gene annotations from a gtf file
#' @description function to read in gene annotation from a gtf file and filters these annotation based on gene element, coding information, hromosomes and transript support level.
#' @param path_gtf path to gtf file
#' @param skip number of lines to skip in the gtf file until the first gene annotation.
#' @param coding denotes which genes are filtered based on gene biotype and/or transcript biotype information. Available options: "protein_coding", c("processed_transcript", "protein_coding"), etc.
#' @param filter_reg_chr logical. If \code{T}, only genes are included which are on regular chromosomes.
#' @param TSL logical. If \code{T}, retrieve transcript support level(TSL) for given transcripts.
#' @param TSLfilt number. If not \code{NULL} but a number between 1-5 (or a vector containing this values), gene_element is filtered to only contain TSL %in% \code{TSLfilt}.
#' @param filter_transcript_biotype logical. If \code{T}, filtering based on transcript biotype.
#' @return list containing two slots of annotation filtered for transcripts and \code{coding} as well as annotation collased per gene into gene loci.
#'   \item{annotation}{annotation filtered for transcripts, gene_biotype = \code{coding} and/or transcript_biotype = \code{coding}}
#'   \item{gelement_coding}{annotation filtered for transcripts, gene_biotype = \code{coding} and/or transcript_biotype = \code{coding} and optional TSL. In addition, transcipts of the same gene are collpased per gene into gene loci}
#' @examples
#' anno <- get_annotation(gtf)
#' anno_test <- get_annotation(gtf, TSL=T)
#' anno_test <- get_annotation(gtf, TSL=T, TSLfilt=1)
#' anno_test <- get_annotation(gtf, TSL=T, TSLfilt=c(1,2))
#' @export
get_annotation <- function(path_gtf, skip=5, coding="protein_coding", filter_reg_chr=T, TSL=F,TSLfilt=NULL, filter_transcript_biotype=F) {
  start_time2 <- Sys.time()
  cat("read in gtf", "\n")
  annotations <- read.delim(path_gtf, skip=skip, header=F)
  cat("filter_gene_element", "\n")
  annotations <- annotations[annotations$V3 == "transcript",]
  cat("get_coding_information_and_filter", "\n")
  cat("get_gene_coding_information", "\n")
  gene_biotype <- sapply(annotations$V9, function(x){
    y <- strsplit(x, split=" ")[[1]]
    y <- sub(";", "", y)
    if ("gene_biotype" %in% y) {
      z <- y[which( y == "gene_biotype")+1]
    }
    else{
      z <- "no_gene_biotype"
    }
  })
  gene_biotype <- as.character(gene_biotype)
  annotations <- cbind(annotations, gene_biotype=gene_biotype)
  cat("get_transcript_coding_information", "\n")
  transcript_biotype <- sapply(annotations$V9, function(x){
    y <- strsplit(x, split=" ")[[1]]
    y <- sub(";", "", y)
    if ("transcript_biotype" %in% y) {
      z <- y[which( y == "transcript_biotype")+1]
    }
    else{
      z <- "no_transcript_biotype"
    }
  })
  transcript_biotype <- as.character(transcript_biotype)
  annotations <- cbind(annotations, transcript_biotype=transcript_biotype)
  if (filter_transcript_biotype){
    annotations <- annotations[annotations$transcript_biotype %in% coding,]
  }
  else{
    annotations <- annotations[annotations$gene_biotype %in%  coding | annotations$transcript_biotype %in% coding,]}
  cat("retrieve gene name", "\n")
  gene_name <- sapply(annotations$V9, function(x){
    y <- strsplit(x, split = " ")[[1]]
    y <- sub(";", "", y)
    if ("gene_name" %in% y) {
      z <- y[which( y == "gene_name")+1]
    }
    else{
      z <- y[which( y == "gene_id")+1]
    }
  })
  gene_name <- as.character(gene_name)
  annotations <- cbind(annotations, gene_name=gene_name)
  cat("retrieve gene ID", "\n")
  gene_id <- sapply(annotations$V9, function(x){
    y <- strsplit(x, split = " ")[[1]]
    y <- sub(";", "", y)
    z <- y[which( y == "gene_id")+1]
  })
  gene_id <- as.character(gene_id)
  annotations <- cbind(annotations, gene_id=gene_id)
  if (TSL){
    cat("retrieve gene TSL", "\n")
    get_tsl <- sapply(annotations$V9, function(x){
      y <- strsplit(x, split = " ")[[1]]
      y <- sub(";", "", y)
      if ("transcript_support_level" %in% y) {
        z <- y[which( y == "transcript_support_level")+1]
      }
      else{
        z <- "no_TSL"
      }
    })
    annotations <- cbind(annotations, TSL=as.character(get_tsl))
  }
  
  colnames <- c("seqnames", "genome_build", "gene_region", "start", "end", "score", "strand", "frame", "gene_info", "gene_biotype", "transcript_biotype", "gene_name", "gene_id")
  if (TSL) {
    colnames(annotations) <- c(colnames, "TSL")
    gene_element <- annotations[, c("seqnames", "start", "end", "strand", "gene_id", "gene_name", "gene_biotype", "transcript_biotype","TSL")]
    if ( !is.null(TSLfilt)) {
      gene_element <- gene_element[gene_element$TSL %in% TSLfilt,] }
  }
  else{
    colnames(annotations) <- colnames
    gene_element <- annotations[, c("seqnames", "start", "end", "strand", "gene_id", "gene_name", "gene_biotype", "transcript_biotype")]
  }
  TSS_pos <- rep(0, nrow(gene_element))
  TSS_pos[which(gene_element$strand == "+")] <- gene_element$start[which(gene_element$strand == "+")]
  TSS_pos[which(gene_element$strand == "-")] <- gene_element$end[which(gene_element$strand == "-")]
  TSS_pos <- as.character(TSS_pos)
  gene_element <- cbind(gene_element, TSSrange=TSS_pos, TSSset=TSS_pos)
  
  non_unique <- names(table(gene_element$gene_name)[table(gene_element$gene_name) != 1 ])
  cat("collapse genes to locus", "\n")
  
  cores <- as.numeric(availableCores() -2)
  cat("available_cores:", cores, "\n")
  start_time <- Sys.time()
  plan(multisession,workers = cores )
  
  non_unique_collapse <- future.apply::future_lapply(seq_along(non_unique), function(x, non_unique, gene_element, TSL) {
    index <- which(gene_element$gene_name == non_unique[x])
    uni <- gene_element[index[1],]
    gene_id <- gene_element$gene_id[index]
    tss_pos_col <- as.numeric(gene_element$TSSrange[index])
    tss_pos_col <- tss_pos_col[order(tss_pos_col)]
    if (TSL) {
      tsl_uni <- length(unique(gene_element$TSL[index]))
      if (tsl_uni != 1) {
        uni[9] <- paste0(gene_element$TSL[index], collapse = "|")}
      uni[10] <- paste0(tss_pos_col[1],"-", tss_pos_col[length(tss_pos_col)])
      uni[11] <- paste(tss_pos_col, collapse = "|")
    }
    else{
      uni[9] <- paste0(tss_pos_col[1],"-", tss_pos_col[length(tss_pos_col)])
      uni[10] <- paste(tss_pos_col, collapse = "|")
    }
    uni[2] <- min(gene_element[index, "start"])
    uni[3] <- max(gene_element[index, "end"])
    uni[5] <- paste(unique(gene_id), collapse = "|")
    return(uni)
  }, non_unique=non_unique, gene_element=gene_element, TSL=TSL)
  end_time <- Sys.time()
  cat(paste("done", "time", difftime(end_time, start_time, units="secs"), "s", "\n", sep = " "))
  
  non_unique_collapse <- do.call(rbind, non_unique_collapse)
  
  gene_element <- rbind(gene_element[!gene_element$gene_name %in% non_unique,], non_unique_collapse)
  
  if ( filter_reg_chr ) {
    cat("filtering for genes on regular chromosomes", "\n")
    gene_element <- gene_element[grep("^chr|^X$|^Y$|^MT$|^M$|^\\d$|^\\d\\d$", gene_element$seqnames),]
    cat(dim(gene_element), "\n")
  }
  anno_list <- list(annotations[order(annotations$gene_name),], gene_element[order(gene_element$gene_name),])
  names(anno_list) <- c("annotation", paste("transcript", paste(coding, collapse = "_"), sep = "_"))
  end_time2 <- Sys.time()
  cat(paste("overall computing", "time", difftime(end_time2, start_time2, units="secs"), "s", "\n", sep = " "))
  return(anno_list)
}

#' @title filter out genes spanning a too large range.
#' @description this function filters out genes which span a region larger than a given threshold.
#' @param annotations list. Output from \code{get_annotation}.
#' @param gene_element data.frame of gene annotations. Either second slot of \code{get_annotation} or similar annotation with second column representing start and third column representing end of a peak.
#' @param tolarge number. Indicating maximal range of genes to be included.
#' @return data.frame of filtered annotations.
#' @examples
#' anno_prolong <- too_large(gene_element = anno[[2]])
#' @export
too_large <- function(annotations=NULL,gene_element=NULL, tolarge=5e+05 ) {
  if ( is.null(annotations) && is.null(gene_element)) {stop("either annotations or gene_element has to be provided")}
  if ( !is.null(annotations) && !is.null(gene_element)) {stop("either annotations OR gene_ement has to be provided")}
  if (!is.null(annotations)) {
    data <- annotations[[2]]
  }
  else{
    data <- gene_element
  }
  smaller <- apply(data, 1, function(x) {
    length(x[2]:x[3]) < tolarge
  })
  transcript_small <- data[smaller,]
  transcript_small
}


#' @title prolong gene range upstream of TSS.
#' @description this function prolongs the gene ranges in the annotation upstream of TSS by a given number to include the promoter region or other proximal regulatory elements.
#' @param annotations list. Output from get \code{get_annotation}.
#' @param gene_element data.frame of gene annotations. Either second slot of \code{get_annotation}, output of \code{too_large} or similar annotation with second column representing start and third column representing end of a peak.
#' @param prolong number. Indicating the range gene is prolonged upstream of TSS.
#' @return data.frame of upstream prolonged annotations.
#' @examples
#' anno_prolong <- prolong_upstream(gene_element = anno_prolong)
#' @export
prolong_upstream <- function(annotations=NULL,gene_element=NULL, prolong=2000) {
  if ( is.null(annotations) && is.null(gene_element)) {stop("either annotations or gene_element has to be provided")}
  if ( !is.null(annotations) && !is.null(gene_element)) {stop("either annotations OR gene_ement has to be provided")}
  if (!is.null(annotations)) {
    data <- annotations[[2]]
  }
  else{
    data <- gene_element
  }
  prolonged_transcript <- t(apply(data, 1, function(x){
    if(as.character(x[4]) == "-") {
      x[3] <- as.numeric(x[3])+prolong
    }
    else{
      x[2] <- as.numeric(x[2])-prolong
    }
    x
  }))
  data.frame(prolonged_transcript)
}

#' @title Annotate peaks which fall onto genes.
#' @description This function, annotates given peaks with the gene(s) they fall onto. Input is as character vector with chromosome, start- end end position. It can be further annotated whether peaks are overlapping the TSS or gene body. Furthermore the range the peak spans before the TSS is computed if the peak overlaps TSS or the distance from peak start and peak end is computed if gene falls onto the gene body.
#' @param peak_features character vector, e.g. rownames of peak-cell matrix, with chromosome information as well as start and end position separated e.g. by ":" and "-" chr1:1000-2000.
#' @param annotations list. Output from get \code{get_annotation}.
#' @param gene_element data.frame of gene annotations. Either second slot of \code{get_annotation}, output of \code{too_large}, \code{prolong_upstream} or similar annotation with second column representing start and third column representing end of a peak.
#' @param split character vector(or object which can be coerced to such), used for splitting each element of \code{peak_features} into a vector with 3 elements containing chromosome information, start= and end position of the peak.
#' @param TSSmode logical. If \code{T},assigns whether annotated gene overlaps a TSS and the range the peak spans before TSS, i.e. distance peak-start to TSS (+strand) or peak-end to TSS (-strand). Default \code{T}. If \code{F}, required for give_activity, only annotates whether peaks overlap genes and if the peak falls onto the spanning ranges of two genes, gene names are collapse with "|".
#' @return data.frame with a row for each peak, and columns containing chromosome information, start- and end position, gene on which the peak is falling or "nomatch" if peak is falling on no gene, strand information of gene peak falls onto, information whether peak is overlapping TSS, gene body or gene end and distance to TSS.
#' @examples
#' peak_genes <- peaks_on_gene(peak_features = rownames(merged_atac_filt), annotations = anno)
#' @export
peaks_on_gene <- function(peak_features,annotations=NULL, gene_element=NULL, split="[-:]", TSSmode=T) {
  start_time2 <- Sys.time()
  if ( is.null(annotations) && is.null(gene_element)) {stop("either annotations or gene_element has to be provided")}
  if ( !is.null(annotations) && !is.null(gene_element)) {stop("either annotations OR gene_ement has to be provided")}
  if (!is.null(annotations)) {
    data <- annotations[[2]]
  }
  else{
    data <- gene_element
  }
  cat("build_peak_table", "\n")
  peaks <- lapply(peak_features, function(x) {
    strsplit(x, split = "[-:]")[[1]]
  })
  peak_features <- c()
  if( class(as.numeric(peaks[[1]][2])) != "numeric" ) {stop("peak_features need to contain numeric start and end")}
  if( class(as.numeric(peaks[[1]][3])) != "numeric" ) {stop("peak_features need to contain numeric start and end")}
  peaks <- do.call(rbind, peaks)
  ret <- stardardize_chr_names(peaks, data)
  peaks <- ret$peaks
  data <- ret$annotation
  chromosomes <- ret$chromosomes
  ret <- list()
  cat("build_chrom_index", "\n")
  gene_chrom_index <- list()
  for ( i in 1:length(chromosomes)){
    gene_chrom_index[[chromosomes[i]]] <- list()
    index <- which(as.character(data$seqnames) == chromosomes[i])
    gene_chrom_index[[chromosomes[i]]][["starts"]] <- as.numeric(data$start)[index]
    gene_chrom_index[[chromosomes[i]]][["ends"]] <- as.numeric(data$end)[index]
    gene_chrom_index[[chromosomes[i]]][["gene_names"]] <-  data$gene_name[index]
    gene_chrom_index[[chromosomes[i]]][["strand"]] <-  data$strand[index]
    gene_chrom_index[[chromosomes[i]]][["TSSset"]] <-  data$TSSset[index]
  }
  cat("start_overlapping_peaks", "\n")
  cores <- as.numeric(availableCores() -2)
  cat("available_cores:", cores, "\n")
  plan(multisession,workers = cores )
  computing <- cores*1000
  peak_list <- create_data_chunks(peaks, computing)

  peaks <- c()
  cat("processing_", min(nrow(peak_list[[1]]), computing), "_rows_at_a_time", "\n", sep="")
  n_processing <- 0
  for (n in 1:length(peak_list)) {
    n_processing <- n_processing +  nrow(peak_list[[n]])
    cat("processing_rows " , n*computing-(computing-1), " to ", min(n_processing, n*computing), " ")
    start_time <- Sys.time()
    peak_list[[n]] <- future.apply::future_lapply(seq_along(1:nrow(peak_list[[n]])), function(x, peak,chr_index) {
      chr <- peak[x,1]
      gene_starts <- chr_index[[chr]][["starts"]]
      gene_ends <- chr_index[[chr]][["ends"]]
      peak_start <- as.numeric(peak[x,2])
      peak_end <- as.numeric(peak[x,3])
      
      left1 <- gene_starts >= peak_start ### peak overlap left or span the whole gene
      left2 <- gene_starts <= peak_end
      right1 <- gene_ends >= peak_start ### peak overlap right or span the whole gene
      right2 <- gene_ends <= peak_end
      mid1 <- gene_starts <= peak_start ### peak on the gene
      mid2 <- gene_ends >= peak_end
      
      gene_left <- chr_index[[chr]][["gene_names"]][left1 & left2]
      gene_right <- chr_index[[chr]][["gene_names"]][right1 & right2]
      gene_mid <- chr_index[[chr]][["gene_names"]][mid1 & mid2]
      if (TSSmode) { #### ad distance here TSS here
        if ( length(gene_left) > 0  ) {
          index <- which(left1 & left2)
          distances <- chr_index[[chr]][["TSSset"]][index]
          distances <- strsplit(distances, "\\|")
          distances <- lapply(distances, function(x) {as.numeric(x) - peak_start})
          distances <- lapply(distances, function(x) {paste(x, collapse = "|")})
          distances <- unlist(distances)
          gene_left <- paste0(gene_left, "_",chr_index[[chr]][["strand"]][index],"_","TSS.overlap.dist.Pstart","_",as.character(distances),"_", ".")
          gene_left <- sub("\\_-_TSS.overlap.dist.+", "_-_._._.", gene_left)
        } ### do also for right and mid
        if (length(gene_right) > 0) {
          index <- which(right1 & right2)
          distances <- chr_index[[chr]][["TSSset"]][index]
          distances <- strsplit(distances, "\\|")
          distances <- lapply(distances, function(x) { peak_end - as.numeric(x) })
          distances <- lapply(distances, function(x) {paste(x, collapse = "|")})
          distances <- unlist(distances)
          gene_right <- paste0(gene_right, "_",chr_index[[chr]][["strand"]][index],"_" ,"TSS.overlap.dist.Pend.","_",as.character(distances), "_", ".")
          gene_right <- sub("\\_+_TSS.overlap.dist.+", "_+_._._.", gene_right)
        }
        if (length(gene_mid) > 0) {
          ### for alternative TSS
          index <- which(mid1 & mid2)
          distances <- chr_index[[chr]][["TSSset"]][index]
          distances <- strsplit(distances, "\\|")
          distances1 <- lapply(distances, function(x) {as.numeric(x) - peak_start})
          vorzeichenS <- lapply(distances1, function(x) { as.numeric(x>0) })
	  vorzeichenS_0 <- lapply(distances1, function(x) { as.numeric(x==0) })
          distances1 <- lapply(distances1, function(x) {paste(x, collapse = "|")})
          distances1 <- unlist(distances1)
          distances2 <- lapply(distances, function(x) {as.numeric(x) - peak_end})
          vorzeichenE <- lapply(distances2, function(x) { as.numeric(x>0) })
	  vorzeichenE_0 <- lapply(distances2, function(x) { as.numeric(x==0) })
          vorzeichen <- Map(rbind, vorzeichenS, vorzeichenE)
          vorzeichen <- lapply(vorzeichen, function(x) { apply(x, 2, function(y) { as.numeric(as.numeric(y[1] == y[2]) == 0) } )})
          vorzeichen <- Map(rbind, vorzeichen, distances)
	  vorzeichen <- Map(rbind, vorzeichen, vorzeichenS_0)
	  vorzeichen <- Map(rbind, vorzeichen, vorzeichenE_0) 
          vorzeichen <- lapply(vorzeichen, function(x) { paste(x[2,which(x[1,] == "1" | x[3,] == "1" | x[4,] == "1" )], collapse = "|")})
          vorzeichen <- unlist(vorzeichen)
          vorzeichen[which(vorzeichen == "")] <- "." 
          distances2 <- lapply(distances2, function(x) {paste(x, collapse = "|")})
          distances2 <- unlist(distances2)
          gene_mid1 <- paste0(gene_mid, "_",chr_index[[chr]][["strand"]][index],"_","TSS.dist.Pstart.","_",as.character(distances1), "_",vorzeichen)
          gene_mid2 <- paste0(gene_mid, "_", chr_index[[chr]][["strand"]][index],"_","TSS.dist.Pend.","_",as.character(distances2), "_",vorzeichen)
          gene_mid <- c(gene_mid1, gene_mid2)
        }
        gene_names <- c(gene_left, gene_right, gene_mid)
        if (length(gene_names) > 1) { gene_names <- unique(gene_names)}
        if (length(gene_names) == 1) {
          return(c(peak[x,],unlist(strsplit(gene_names, split = "_"))))
        }
        else if (length(gene_names) == 0) {
          return(c(peak[x,], "nomatch", ".", ".", ".","."))
        }
        else {
          to_return <- matrix(rep(peak[x,],length(gene_names)),nrow = length(gene_names), byrow = T)
          to_return <- cbind(to_return, t(sapply(gene_names, function(x) { y <- strsplit(x, split="_"); y <- unlist(y) })))
          return(to_return)
        }
      }
      else{
        gene_names <- c(gene_left, gene_right, gene_mid)
        if (length(gene_names) > 1) { gene_names <- unique(gene_names)}
        if (length(gene_names) == 1) {
          return(c(peak[x,],gene_names))
        }
        else if (length(gene_names) == 0) {
          return(c(peak[x,], "nomatch"))
        }
        else {
          overlap_nam <- paste(gene_names, collapse = "|")
          return(c(peak[x,], overlap_nam))
        }
      }
    }, chr_index=gene_chrom_index, peak=peak_list[[n]])
    peak_list[[n]] <- do.call(rbind, peak_list[[n]])
    end_time <- Sys.time()
    cat(paste("done", "time", difftime(end_time, start_time, units="secs"), "s", "\n", sep = " "))
  }
  peak_ongene <- do.call(rbind, peak_list)
  if (TSSmode==T){
    rownames(peak_ongene) <- paste0(peak_ongene[,1],":", peak_ongene[,2], "-", peak_ongene[,3])
    colnames(peak_ongene) <- c("seqnames", "Pstart", "Pend", "peak_on_gene", "strand", "TSSinfo", "TSSdistance", "overlap.alter.TSS")
    peaks_length <- apply(peak_ongene, 1, function(x) {
      length(as.numeric(x[2]):as.numeric(x[3]))
    })
    peak_ongene<- cbind(peak_ongene, peaks_length=peaks_length)
  }
  else{
    rownames(peak_ongene) <-  paste0(peak_ongene[,1],":", peak_ongene[,2], "-", peak_ongene[,3])}
  end_time2 <- Sys.time()
  cat(paste("overall computing", "time", difftime(end_time2, start_time2, units="secs"), "s", "\n", sep = " "))
  return(peak_ongene)
}


#' @title Aggregate peak counts falling onto the same gene to gene activities
#' @description This function takes the output of \code{peaks_on_gene}, run with TSSmode = \code{F} and the peak-cell count matrix as input and aggregates peaks falling onto the same genes to a gene-activity count matrix.
#' @param gene_peaks data.frame of peak annotations. Output of \code{peak_on_gene} with peak annotated to the genes they fall onto.
#' @param peak_matrix Peak-cell count matrix. With rows representing peaks and columns representing cells.
#' @return sparse matrix of aggregated activity counts with rows as genes and columns as cells.
#' @examples
#' activity <- give_activity(peak_genes, merged_atac_filt)
#' @export
give_activity <- function(gene_peaks, peak_matrix) {
  start_time2 <- Sys.time()
  cat("filter_for_peaks_overlapping_gene", "\n")
  peak_matrix <- peak_matrix[gene_peaks[,4] != "nomatch",]
  gene_peaks <- gene_peaks[gene_peaks[,4] != "nomatch",]
  list_peak <- list()
  cat("filter_for_peaks_which_overlap_more_than_1_gene", "\n")
  index_overlap <- which(grepl("\\|", gene_peaks[,4]))
  activity_overlap <- peak_matrix[index_overlap,]
  gene_peaks_overlap <- gene_peaks[index_overlap,]
  gene_peaks <- gene_peaks[-index_overlap,]
  peak_matrix <- peak_matrix[-index_overlap,]
  cat("split for parallelization", "\n")
  cores <- as.numeric(availableCores() -2)
  cat("available_cores:", cores, "\n")

  computing <- cores*100
  peak_list <- create_data_chunks(gene_peaks_overlap, computing)
  activity_list <- create_data_chunks(activity_overlap, computing)

  plan(multisession,workers = cores )
  activity_overlap <- c()
  gene_peaks_overlap <- c()
  n_processing <- 0
  for (n in 1:length(peak_list)) {
    cat("assign_counts_of_peaks_overlapping_genes_to_individual_genes", "\n")
    n_processing <- n_processing +  nrow(activity_list[[n]])
    cat("processing_rows ", n*computing-(computing-1), " to ", min(n_processing, n*computing), " ")
    start_time <- Sys.time()
    
    peak_list[[n]] <- future_lapply(seq_along(1:nrow(peak_list[[n]])), function(x, activity, peaks) {
      gene_name <- peaks[x, 4]
      gene_name <- as.character(unlist(strsplit(gene_name, split="\\|")))
      mat_activity <- matrix(rep(activity[x,], length(gene_name)), nrow = length(gene_name), byrow = T)
      mat_activity <- list(mat=mat_activity, gene_name=gene_name )
    }, activity=as.matrix(activity_list[[n]]), peaks=as.matrix(peak_list[[n]]))
    
    peak_list_mat <-  lapply(peak_list[[n]], function(x) {
        return(x[["mat"]])
      })
   
    peak_list[[n]] <- lapply(peak_list[[n]], function(x) {
        return(x[["gene_name"]])
    })
    peak_list_mat <- do.call(rbind, peak_list_mat)
    peak_list_mat <- Matrix(as.matrix(peak_list_mat), sparse = T)
    peak_list[[n]] <- Reduce(append, peak_list[[n]])  
    peak_list[[n]] <- list(mat = peak_list_mat, gene_name=peak_list[[n]])
    peak_list_mat <- c()
    end_time <- Sys.time()
    cat(paste("done", "time", difftime(end_time, start_time, units="secs"), "s", "\n", sep = " "))
  }
  cat("append genes and rbind tables", "\n")
  activity_overlap_comb <- lapply(peak_list, function(x) {
      return(x[["mat"]])
  })
  peak_list <- lapply(peak_list, function(x) {
      return(x[["gene_name"]])
  })
    
  activity_overlap_comb <- do.call(rbind, activity_overlap_comb)
  peak_list <- Reduce(append, peak_list)
  
  genes_overlap <- c(as.character(gene_peaks[,4]), peak_list )
  activity_gene_add <- rbind(peak_matrix,activity_overlap_comb)
  rownames(activity_gene_add) <- genes_overlap
  cat("aggregate counts for peaks on the same gene", "\n")
  activity_gene_add_sum <- aggregate.Matrix(activity_gene_add, as.factor(rownames(activity_gene_add)), fun = "sum")
  end_time2 <- Sys.time()
  cat(paste("overall computing", "time", difftime(end_time2, start_time2, units="secs"), "s", "\n", sep = " "))
  return(activity_gene_add_sum)
}


#' @title Annotates peaks to closest gene.
#' @description This function takes the output of \code{peaks_on_gene} and annotates peaks which do not fall onto genes to its closest downstream as well as closest gene in general(if closest gene is not downstream of the peak). The distances to closest (downstream) gene of each peak is also calculated and returned.
#' @param peaks data.frame of peak annotations. Output of \code{peaks_on_gene} with peak annotated to the genes they fall onto.
#' @param annotations list. Output from get \code{get_annotation}.
#' @param gene_element data.frame of gene annotations. Either second slot of \code{get_annotation}, output of \code{too_large}, \code{prolong_upstream} or similar annotation with second column representing start and third column representing end of a peak.
#' @param TSSmode logical. If \code{T} specify if \code{peaks_on_gene} function has been run in TSSmode. ### can be optimized to infer this automatic.
#' @return data.frame with a row for each peak, and columns of peaks_on_gene function and in addition ,closest downstream gene, closest gene in general, distance to gene edge of closest downstream gene and closest gene in general.
#' @examples
#' closest_gene <- peaks_closest_gene(peak_genes,anno)
#' @export
peaks_closest_gene <- function(peaks, annotations=NULL, gene_element=NULL, TSSmode=T, fast=F) {
  start_time2 <- Sys.time()
  if ( is.null(annotations) && is.null(gene_element)) {stop("either annotations or gene_element has to be provided")}
  if ( !is.null(annotations) && !is.null(gene_element)) {stop("either annotations OR gene_ement has to be provided")}
  if (!is.null(annotations)) {
    data <- annotations[[2]]
  }
  else{
    data <- gene_element
  }
  ret <- stardardize_chr_names(peaks, data)
  peaks <- ret$peaks
  data <- ret$annotation
  chromosomes <- ret$chromosomes
  ret <- list()
  cat("build_chrom_index", "\n")
  gene_chrom_index <- list()
  for ( i in 1:length(chromosomes)){
    gene_chrom_index[[chromosomes[i]]] <- list()
    index <- which(as.character(data$seqnames) == chromosomes[i])
    gene_chrom_index[[chromosomes[i]]][["starts"]] <- as.numeric(data$start)[index]
    gene_chrom_index[[chromosomes[i]]][["ends"]] <- as.numeric(data$end)[index]
    gene_chrom_index[[chromosomes[i]]][["gene_names"]] <-  data$gene_name[index]
    gene_chrom_index[[chromosomes[i]]][["strand"]] <- data$strand[index]
  }
  cat("separate into peaks with no match and peaks which showed overlap", "\n")
  if (fast==F) {
  peaks_annotated <- peaks[peaks[,4] != "nomatch",]
  if ( TSSmode==T){
    cat("cbind_peaks_on_gene", "\n")
    peaks_annotated_end <- peaks_annotated[grep("TSS",peaks_annotated[,6], invert = T), , drop=FALSE]
    peaks_annotated_end <- cbind(peaks_annotated_end, closest_downstream_gene=rep("",nrow(peaks_annotated_end)) ,closest_gene=rep("",nrow(peaks_annotated_end)),dist_clos_downstream=rep("",nrow(peaks_annotated_end)), dist_gene.edge_closest=rep("",nrow(peaks_annotated_end)))
    cat("cbind_peaks_on_TSS", "\n")
    peaks_annotated_TSS_ol <- peaks_annotated[grep("TSS.overlap",peaks_annotated[,6]), , drop=FALSE]
    peaks_annotated_TSS_ol <- cbind(peaks_annotated_TSS_ol, closest_downstream_gene=peaks_annotated_TSS_ol[,4],closest_gene=peaks_annotated_TSS_ol[,4], dist_clos_downstream=rep("",nrow(peaks_annotated_TSS_ol)),dist_gene.edge_closest=rep("",nrow(peaks_annotated_TSS_ol)))
    peaks_annotated_TSS_on <- peaks_annotated[grep("TSS.dist",peaks_annotated[,6]), , drop=FALSE]
    peaks_annotated_TSS_on <- cbind(peaks_annotated_TSS_on, closest_downstream_gene=peaks_annotated_TSS_on[,4],closest_gene=peaks_annotated_TSS_on[,4], dist_clos_downstream=rep("",nrow(peaks_annotated_TSS_on)),dist_gene.edge_closest=rep("",nrow(peaks_annotated_TSS_on)))

  }
  else{
    peaks_annotated <- cbind(peaks_annotated, closest_downstream_gene=peaks_annotated[,4] ,closest_gene=peaks_annotated[,4])}
  }
  peaks_not_annotated <- peaks[peaks[,4] == "nomatch",]
  peaks_nam <- rownames(peaks)
  peaks_nam_not_annotated <- rownames(peaks_not_annotated)
  peaks <- c()
  cat("start_looking_for_closest_gene", "\n")
  cores <- as.numeric(availableCores() -2)
  cat("available_cores:", cores, "\n")

  computing <- cores*1000
  peak_list <- create_data_chunks(peaks_not_annotated, computing)

  peaks_not_annotated <- c()
  plan(multisession,workers = cores )
  n_processing <- 0
  for (n in 1:length(peak_list)) {

    n_processing <- n_processing +  nrow(peak_list[[n]])
    cat("processing_rows " , n*computing-(computing-1), " to ", min(n_processing, n*computing), " ")
    start_time <- Sys.time()
    peak_list[[n]] <- future.apply::future_lapply(seq_along(1:nrow(peak_list[[n]])), function(x, peak,chr_index,TSSmode) {

      chr <- peak[x,1]
      gene_starts <- chr_index[[chr]][["starts"]]
      gene_ends <- chr_index[[chr]][["ends"]]

      left <- gene_starts - as.numeric(peak[x,3]) ### upstream +
      right <- as.numeric(peak[x,2]) - gene_ends ### upstream -
      abs_left <- abs(left) ### shortest distance to gene start
      abs_right <- abs(right) ### shortest distance to gene ends
      left[left < 0] <- Inf  ## if < 0 peak starts after genes start
      right[right < 0] <- Inf ## if < 0 peak starts before genes end
      strand <- chr_index[[chr]][["strand"]]
      left[strand == "-"] <- Inf  ### peak before gene start but genes on the - strand => not upstream
      right[strand == "+"] <- Inf ### peak after gene end but genes on the + strand => not upstream
      if (TSSmode) {
        if (min(left) < min(right)) {
          gene_downstream <- paste(chr_index[[chr]][["gene_names"]][which(left == min(left))], collapse = "|")#,min(left), sep = "_")
          distance <- min(left)
        }
        else if (min(left) > min(right)) {
          gene_downstream <- paste(chr_index[[chr]][["gene_names"]][which(right == min(right))], collapse = "|")#, min(right), sep = "_")
          distance <- min(right)
        }
        else{
          gene_downstream <- "no_downstream"
          distance <- ""
        }
        if (min(abs_left) < min(abs_right)) {
          gene_general <- paste(chr_index[[chr]][["gene_names"]][which(abs_left == min(abs_left))], collapse = "|")#, min(abs_left), sep = "_")
          dist_general <- min(abs_left)
        }
        else{
          gene_general <- paste(chr_index[[chr]][["gene_names"]][which(abs_right == min(abs_right))], collapse = "|")#, min(abs_right), sep = "_")
          dist_general <- min(abs_right)
        }
        return(c(peak[x,], gene_downstream, gene_general, distance, dist_general))}
      else{
        if (min(left) < min(right)) {
          gene_downstream <- paste(paste(chr_index[[chr]][["gene_names"]][which(left == min(left))], collapse = "|"),min(left), sep = "_")
        }
        else if (min(left) > min(right)) {
          gene_downstream <- paste(paste(chr_index[[chr]][["gene_names"]][which(right == min(right))], collapse = "|"), min(right), sep = "_")
        }
        else{
          gene_downstream <- "no_downstream"
        }
        if (min(abs_left) < min(abs_right)) {
          gene_general <- paste(paste(chr_index[[chr]][["gene_names"]][which(abs_left == min(abs_left))], collapse = "|"), min(abs_left), sep = "_")
        }
        else{
          gene_general <- paste(paste(chr_index[[chr]][["gene_names"]][which(abs_right == min(abs_right))], collapse = "|"), min(abs_right), sep = "_")
        }
        return(c(peak[x,], gene_downstream, gene_general))
      }

    },chr_index=gene_chrom_index, peak=peak_list[[n]], TSSmode=TSSmode )
    peak_list[[n]] <- do.call(rbind, peak_list[[n]])
    end_time <- Sys.time()
    cat(paste("done", "time", difftime(end_time, start_time, units="secs"), "s", "\n", sep = " "))
  }
  peaks_close <- do.call(rbind,peak_list)
  rownames(peaks_close) <- peaks_nam_not_annotated
  cat(dim(peaks_close), "\n")
  if (fast == F) {
  if (TSSmode) {
    cat(dim(peaks_annotated_TSS_ol), "\n")
    cat(dim(peaks_annotated_TSS_on), "\n")
    peaks_on_gene <- rbind(peaks_annotated_TSS_ol, peaks_annotated_TSS_on)
    peaks_on_gene <- rbind(peaks_on_gene, peaks_annotated_end)
    cat(dim(peaks_on_gene), "\n")
    peaks_close <- rbind(peaks_on_gene, peaks_close)
    cat(dim(peaks_close), "\n")
  }
  else{
    cat(dim(peaks_annotated), "\n")
    cat(dim(peaks_close), "\n")
    peaks_close <- rbind(peaks_annotated, peaks_close)}
  }
  colnames(peaks_close) <- c("seqnames",                "Pstart"                 , "Pend"   , "peak_on_gene",           "strand"        ,          "TSSinfo"                ,"TSSdistance"       ,   "overlap.alter.TSS"   ,  "peaks_length"  ,         "closest_downstream_gene" ,"closest_gene" , "dist_clos_downstream"   , "dist_gene.edge_closest")
  end_time2 <- Sys.time()
  cat(paste("overall computing", "time", difftime(end_time2, start_time2, units="secs"), "s", "\n", sep = " "))
  return(peaks_close)
}


#' @title Get combined overlapping peak set.
#' @description This function takes the peaks which have been called on different samples, and thus could represent overlapping peaks with different start and or end positions to provide an overlapping peak set.
#' @param atac_layer list of peak-cell count matrix of the different samples, with rows as peaks and cells as columns.
#' @param filter_reg_chr logical. If \code{T}, only genes are included which are on regular chromosomes..
#' @return GenomicRanges object of combined peaks of different samples of atac_layer list.
#' @examples
#' combined.peaks <- give_combined_peaks(atac_layer)
#' @export
give_combined_peaks <- function(atac_layer, filter_reg_chr=T) {
  cat("build_peak_table for each sample", "\n")
  atac_layer_peaks <- lapply(atac_layer, function(x) {
    peaks <- rownames(x)
    peaks <- lapply(peaks, function(x) {
      strsplit(x, split = "[-:]")[[1]]
    })
    if( class(as.numeric(peaks[[1]][2])) != "numeric" ) {stop("peak_features need to contain numeric start and end")}
    if( class(as.numeric(peaks[[1]][3])) != "numeric" ) {stop("peak_features need to contain numeric start and end")}
    peaks <- do.call(rbind, peaks)
    colnames(peaks) <- c("chr", "start", "end")
    return(data.frame(peaks))
  })
  if (filter_reg_chr==T){
    cat("filter for peaks on regular chromosomes", "\n")
    atac_layer_peaks_gr <- lapply(atac_layer_peaks, function(x){
      y <- x[grep("chr", x[,1]),]
    })
  }
  cat("make genomic ranges for peaks of each sample", "\n")
  atac_layer_peaks_gr <- lapply(atac_layer_peaks_gr, function(x){
    y <- makeGRangesFromDataFrame(x)
  })
  cat("reduce peaks of different samples to comman overapping peak set", "\n")
  to_reduce <- paste0("c(", paste0("atac_layer_peaks_gr", "[[", 1:length(atac_layer_peaks_gr), "]]", collapse = "," ), ")")
  combined.peaks <- reduce(x = eval(parse(text = to_reduce)))
  return(combined.peaks)
}


#' @title Annotate peaks to overlapping peaks
#' @description This function annotates peaks which have been called on different samples, and thus could represent overlapping peaks with different start and or end positions to provided overlapping peaks. In addition a provided peak-cell count matrix can be aggregated across the common overlapping peak set.
#' @param peak_features character vector, e.g. rownames of peak-cell matrix, with chromosome information as well as start and end position separated by e.g. by "-" and ":".
#' @param combined.peaks GenomicRanges object of combined peaks. Output of \code{give_combined_peaks} function.
#' @param do.aggregate logical. If \code{T}, aggregate provided peak-cell count matrix across the common overlapping peak set.
#' @param peak_matrix Peak-cell count matrix. With rows representing peaks and columns representing cells.
#' @param insert_run1 output of function with do.aggregate = \code{F}, i.e. matrix with a row for each peak, containing chromosome information, start- and end position and annotated peak with combined peak set
#' @return If \code{do.aggregate} = \code{F}: matrix with a row for each peak, containing chromosome information, start- and end position and annotated peak with combined peak set. If \code{do.aggregate} = \code{T} and \code{peak_matrix} provided:
#'   \item{agg.comb.peaks}{aggregated peak-cell count matrix across the common overlapping peak set.}
#'   \item{peaks.combined}{data.frame with a row for each peak, containing chromosome information, start- and end position and annotated peak with combined peak set.}
#' @examples
#' overlapped_peaks <- peak_overlap(peak_features=rownames(merged_atac_filt), combined.peaks=combined.peaks)
#' overlapped_peaks <- peak_overlap(do.aggregate=T, peak_matrix=merged_atac_filt,insert_run1 = overlapped_peaks)
#' @export

peak_overlap <- function(peak_features=NULL, combined.peaks=NULL, do.aggregate=F,peak_matrix=NULL, insert_run1=NULL ) {
  if (do.aggregate==T) {
    if(is.null(peak_matrix)) { stop("if do.aggregate=T, peak_matrix has to be provided")}
  }
  if (is.null(insert_run1)) {
  if(is.null(peak_features ) || is.null(combined.peaks)) {stop("both 'peak_features' and 'combined.peaks' need to be provided")}
  start_time2 <- Sys.time()
  chromosomes <- unique(as.character(combined.peaks@seqnames))
  gene_chrom_index <- list()
  for ( i in 1:length(chromosomes)){
    gene_chrom_index[[chromosomes[i]]] <- list()
    index <- which(as.character(combined.peaks@seqnames) == chromosomes[i])
    gene_chrom_index[[chromosomes[i]]][["starts"]] <- as.numeric(data.frame(combined.peaks@ranges)$start)[index]
    gene_chrom_index[[chromosomes[i]]][["ends"]] <- as.numeric(data.frame(combined.peaks@ranges)$end)[index]
    gene_chrom_index[[chromosomes[i]]][["gene_names"]] <-   paste0(chromosomes[i], ":", gene_chrom_index[[chromosomes[i]]][["starts"]], "-", gene_chrom_index[[chromosomes[i]]][["ends"]])
  }
  peaks <- lapply(peak_features, function(x) {
    strsplit(x, split = "[-:]")[[1]]
  })
  peak_features <- c()
  if( class(as.numeric(peaks[[1]][2])) != "numeric" ) {stop("peak_features need to contain numeric start and end")}
  if( class(as.numeric(peaks[[1]][3])) != "numeric" ) {stop("peak_features need to contain numeric start and end")}
  peaks <- do.call(rbind, peaks)
  
  cores <- as.numeric(availableCores() -2)
  cat("available_cores:", cores, "\n")
  computing <- cores*1000
  peak_list <- create_data_chunks(peaks, computing)
  peaks <- c()
  plan(multisession,workers = cores )
  n_processing <- 0
  for (n in 1:length(peak_list)) {
    cat("get overlapping peaks", "\n")
    n_processing <- n_processing +  nrow(peak_list[[n]])
    cat("processing_rows " , n*computing-(computing-1), " to ", min(n_processing, n*computing), " ")
    start_time <- Sys.time()
    peak_list[[n]] <- future_lapply(seq_along(1:nrow(peak_list[[n]])), function(x, peaks, chr_index) {
      chr <- peaks[x,1]
      peak_start <- as.numeric(peaks[x,2])
      peak_end <- as.numeric(peaks[x,3])
      starts <- chr_index[[chr]][["starts"]]
      ends <- chr_index[[chr]][["ends"]]
      
      mid1 <- starts <= peak_start ### peak on the gene
      mid2 <- ends >= peak_end
      gene_names <- chr_index[[chr]][["gene_names"]][mid1 & mid2]
      return(c(chr, peak_start, peak_end, gene_names))
    },peaks=peak_list[[n]],chr_index=gene_chrom_index)
    peak_list[[n]] <- do.call(rbind, peak_list[[n]])
    end_time <- Sys.time()
    cat(paste("done", "time", difftime(end_time, start_time, units="secs"), "s", "\n", sep = " "))
  }
  peaks_combined <- do.call(rbind, peak_list)
  peak_list <- c()
  rownames(peaks_combined) <- paste0(peaks_combined[,1],":", peaks_combined[,2], "-", peaks_combined[,3])
  colnames(peaks_combined) <- c("seqnames", "start", "end", "combined_peaks")
  end_time2 <- Sys.time()
  cat(paste("overall computing", "time", difftime(end_time2, start_time2, units="secs"), "s", "\n", sep = " "))
  }
  else{
    peaks_combined <- insert_run1
  }
  if (do.aggregate==T) {
    if (class(peak_matrix) != "dgCMatrix" ) {
      stop("please provide 'peak_matrix' as sparce matrix in dgCMatrix format")
    }
    cat("aggregate overlapping peaks", "\n")
    combined_peaks_matrix <- aggregate.Matrix(peak_matrix, as.factor(peaks_combined[,4]), fun = "sum")
    cat("done aggregate overlapping peaks", "\n")
    combined_list <- list(agg.comb.peaks=combined_peaks_matrix, peaks.combined=peaks_combined)
    return(combined_list)
  }
  else{
    return(peaks_combined)}
}

#' @title Merge peak-cell count matrices.
#' @description This function takes the peak-count matrices for different samples, stored in a list with a slot for every sample and merges these matrices.
#' @param list.samples list with a slot for every sample's peak-cell count matrix
#' @return merged peak-cell count matrix
#' @examples
#' atac_layer_merged <- merge_samples(atac_layer)
#' @export
merge_samples <- function(list.samples) {
  start_time <- Sys.time()
  layer_merge <- lapply(seq_along(list.samples), function(x, list.samples) {
    rang <- 1:length(list.samples)
    rang <- rang[-x]
    feat_others <- c()
    for ( i in rang) {
      feat_others <- c(feat_others, rownames(list.samples[[i]]))
    }
    mat_merge <- matrix(rep(0,as.numeric(length(feat_others))*as.numeric(ncol(list.samples[[x]]))), ncol = ncol(list.samples[[x]]))
    mat_merge <- Matrix(mat_merge, sparse = T)
    rownames(mat_merge) <- feat_others
    mat_merge <- rbind(list.samples[[x]],mat_merge)
    mat_merge <- mat_merge[order(rownames(mat_merge)),]
    mat_merge
  }, list.samples=list.samples)
  end_time <- Sys.time()
  cat(paste("done", "time", difftime(end_time, start_time, units="secs"), "s", "\n", sep = " "))
  layer_merge <- do.call(cbind, layer_merge)
  layer_merge
}

