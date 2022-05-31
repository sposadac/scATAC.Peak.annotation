#' @import future.apply
#' @import Matrix
#' @import Matrix.utils
#' @import GenomicRanges
#' @import future
#'
#' @title  get gene annotations from a gtf file
#' @description function to read in gene annotation from a gtf file and filters these annotation based on gene element, coding information and chromosomes.
#' @param path_gtf path to gtf file
#' @param skip number of lines to skip in the gtf file until the first gene annotation.
#' @param gelement denotes which gene element are filtered. Available options: "transcript", "exon", "intron"
#' @param coding denotes which genes are filtered based on coding information. Available options: "protein_coding",
#' @param filter_reg_chr logical. If \code{T}, only genes are included which are on regular chromosomes.
#' @return list containing two slots of full annotation filtered for \code{gelement} and \code{coding} aswell as annotations collpased per gene into gene loci.
#'   \item{annotation}{annotation filtered for \code{gelement} and \code{coding}}
#'   \item{gelement_coding}{annotation filtered for \code{gelement} and \code{coding} and collpased per gene into gene loci}
#' @examples
#' anno <- get_annotation(gtf)
#' @export
get_annotation <- function(path_gtf, skip=5,gelement="transcript", coding="protein_coding", filter_reg_chr=T) {
  cat("read in gtf", "\n")
  annotations <- read.delim(path_gtf, skip=skip, header=F)
  cat("filter_gene_element", "\n")
  annotations <- annotations[annotations$V3 == gelement,]
  cat("filter_coding_information", "\n")
  gene_biotype <- sapply(annotations$V9, function(x){
    y <- strsplit(x, split=" ")[[1]][length( strsplit(x, split=" ")[[1]])]
    y <- sub(";", "", y)
  })
  gene_biotype <- as.character(gene_biotype)
  annotations <- cbind(annotations, gene_biotype=gene_biotype)
  annotations <- annotations[annotations$gene_biotype ==  coding,]
  cat("retrieve gene name and cbind", "\n")
  gene_name <- sapply(annotations$V9, function(x){
    y <- strsplit(x, split = " ")[[1]]
    if ("gene_name" %in% y) {
      z <- y[which( y == "gene_name")+1]
      z <- sub(";", "", z)
    }
    else{
      z <- y[2]
      z <- sub(";", "", z)
    }
  })
  gene_name <- as.character(gene_name)
  annotations <- cbind(annotations, gene_name=gene_name)
  cat("retrieve gene ID and cbind", "\n")
  gene_id <- sapply(annotations$V9, function(x){
    y <- strsplit(x, split=" ")[[1]][2]
    y <- sub(";", "", y)
  })
  gene_id <- as.character(gene_id)
  annotations <- cbind(annotations, gene_id=gene_id)
  colnames(annotations) <- c("seqnames", "genome_build", "gene_region", "start", "end", "dot", "strand", "dot", "gene_info", "gene_biotype", "gene_name", "gene_id")
  gene_element <- annotations[, c("seqnames", "start", "end", "strand", "gene_id", "gene_name", "gene_biotype")]
  non_unique <- names(table(gene_name)[table(gene_name) != 1 ])
  cat("collapse genes to locus", "\n")
  for (i in 1:length(non_unique)) {
    cat(non_unique[i], "\n")
    index <- which(gene_element$gene_name == non_unique[i])
    cat(length(index), "\n")
    mini <- min(gene_element[index, "start"])
    maxi <- max(gene_element[index, "end"])
    cat(dim(gene_element), "\n")
    gene_element[index[1],"start"] <- mini
    gene_element[index[1], "end"] <- maxi
    gene_element <- gene_element[-index[-1],]
    cat(dim(gene_element), "\n")
  }
  if ( filter_reg_chr ) {
    cat("filtering for genes on regular chromosomes", "\n")	  
    gene_element <- gene_element[grep("chr", gene_element$seqnames),]
  }
  anno_list <- list(annotations, gene_element)
  names(anno_list) <- c("annotation", paste(gelement, coding, sep = "_"))
  return(anno_list)
}

#' @title filter out genes spanning a large range.
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


#' @title prolong gene range proximal.
#' @description this function prolongs the gene ranges in the annotation proximal by a given number to include the promoter region or other proximal regulatory elements.
#' @param annotations list. Output from get \code{get_annotation}.
#' @param gene_element data.frame of gene annotations. Either second slot of \code{get_annotation}, output of \code{too_large} or similar annotation with second column representing start and third column representing end of a peak.
#' @param prolong number. Indicating range the range of gene elements is prolonged proximal.
#' @return data.frame of proximal prolonged annotations.
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
#' @description This function, annotates given peaks with the gene(s) they fall onto. Input is as character vector with chromosome, start- end end position. if the peak falls onto the spanning ranges of two genes, gene names are collapse with "|".
#' @param peak_features character vector, e.g. rownames of peak-cell matrix, with chromosome information as well as start and end position separated e.g. by "-" and ":".
#' @param annotations list. Output from get \code{get_annotation}.
#' @param gene_element data.frame of gene annotations. Either second slot of \code{get_annotation}, output of \code{too_large}, \code{prolong_upstream} or similar annotation with second column representing start and third column representing end of a peak.
#' @param split character vector(or object which can be coerced to such), used for splitting each element of \code{peak_features} into a vector with 3 elements containing chromosome information, start= and end position of the peak.
#' @return data.frame with a row for each peak, containing chromosome information, start- and end position as well as the gene on which the peak is falling or "nomatch" if peak is falling on no gene.
#' @examples
#' peak_genes <- peaks_on_gene2(peak_features = rownames(merged_atac_filt), gene_element = anno_prolong)
#' @export
peaks_on_gene <- function(peak_features,annotations=NULL, gene_element=NULL, split="[-:]") {
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
  if( class(as.numeric(peaks[[1]][2])) != "numeric" ) {stop("peak_features need to contain numeric start and end")}
  if( class(as.numeric(peaks[[1]][3])) != "numeric" ) {stop("peak_features need to contain numeric start and end")}
  peaks <- do.call(rbind, peaks)
  chromosomes <- unique(as.character(data$seqnames))
  cat("build_chrom_index", "\n")
  gene_chrom_index <- list()
  for ( i in 1:length(chromosomes)){
    gene_chrom_index[[chromosomes[i]]] <- list()
    index <- which(as.character(data$seqnames) == chromosomes[i])
    gene_chrom_index[[chromosomes[i]]][["starts"]] <- as.numeric(data$start)[index]
    gene_chrom_index[[chromosomes[i]]][["ends"]] <- as.numeric(data$end)[index]
    gene_chrom_index[[chromosomes[i]]][["gene_names"]] <-  data$gene_name[index]
  }
  cat("start_overlapping_peaks", "\n")
  cores <- as.numeric(future::availableCores() -2)
  cat("available_cores:", cores, "\n")
  plan(multisession,workers = cores )
  peak_list <- list()
  computing <- cores*1000
  data_iterat <- nrow(peaks)%/%computing
  for (n in 1:data_iterat) {
    peak_list[[n]] <- peaks[(n*computing-(computing-1)):(n*computing),]
  }
  if ((nrow(peaks) %% computing) != 0) {
    peak_list[[length(peak_list) +1 ]] <- peaks[(length(peak_list)*computing+1):nrow(peaks),]}
  cat("processing_",computing,"_rows_at_a_time", "\n")
  for (n in 1:length(peak_list)) {
    cat("processing_rows " , n*computing-(computing-1), " to ", n*computing, " ")
    start_time <- Sys.time()
    peak_list[[n]] <- future_lapply(seq_along(1:nrow(peak_list[[n]])), function(x, peak,chr_index) {
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

    }, chr_index=gene_chrom_index, peak=peak_list[[n]])
    end_time <- Sys.time()
    cat(paste("done", "time", difftime(end_time, start_time, units="secs"), "s", "\n", sep = " "))
  }
  peak_list2 <- lapply(peak_list, function(x) { y <- do.call(rbind, x) })
  peak_ongene <- do.call(rbind, peak_list2)
  rownames(peak_ongene) <- peak_features
  end_time2 <- Sys.time()
  cat(paste("overall computing", "time", difftime(end_time2, start_time2, units="secs"), "s", "\n", sep = " "))
  return(peak_ongene)
}

#' @title Aggregate peak counts falling onto the same gene to gene activities
#' @description This function takes the output of \code{peak_on_gene} and the peak-cell count matrix as input and aggregates peaks falling onto the same genes to a gene-activity count matrix.
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
  cat("split for parallelization", "\n")
  cores <- as.numeric(future::availableCores() -2)
  cat("available_cores:", cores, "\n")

  peak_list <- list()
  activity_list <- list()
  computing <- cores*100
  data_iterat <- nrow(activity_overlap)%/%computing

  for (n in 1:data_iterat) {
    peak_list[[n]] <- gene_peaks_overlap[(n*computing-(computing-1)):(n*computing),]
  }
  if ((nrow(gene_peaks_overlap) %% computing) != 0 ) {
    peak_list[[length(peak_list) +1 ]] <- gene_peaks_overlap[(length(peak_list)*computing+1):nrow(gene_peaks_overlap),]}

  for (n in 1:data_iterat) {
    activity_list[[n]] <- activity_overlap[(n*computing-(computing-1)):(n*computing),]
  }
  if ((nrow(activity_overlap) %% computing) != 0 ) {
    activity_list[[length(activity_list) +1 ]] <- activity_overlap[(length(activity_list)*computing+1):nrow(activity_overlap),]}
  plan(multisession,workers = cores )
  activity_overlap <- c()
  gene_peaks_overlap <- c()
  for (n in 1:length(peak_list)) {
    cat("assign_counts_of_peaks_overlapping_genes_to_individual_genes", "\n")
    cat("processing_rows " , n*computing-(computing-1), " to ", n*computing, " ")
    start_time <- Sys.time()

    peak_list[[n]] <- future_lapply(seq_along(1:nrow(peak_list[[n]])), function(x, activity, peaks) {
      gene_name <- peaks[x, 4]
      gene_name <- as.character(unlist(strsplit(gene_name, split="\\|")))
      mat_activity <- matrix(rep(activity[x,], length(gene_name)), nrow = length(gene_name), byrow = T)
      mat_activity <- list(mat=mat_activity, gene_name=gene_name )
    }, activity=as.matrix(activity_list[[n]]), peaks=as.matrix(peak_list[[n]]))
    end_time <- Sys.time()
    cat(paste("done", "time", difftime(end_time, start_time, units="secs"), "s", "\n", sep = " "))
  }
  cat("append genes and rbind tables", "\n")
  activity_overlap_comb <- lapply(peak_list, function(x) {
    y <- lapply(x, function(z) {
      return(z[["mat"]])
    })
    y
  })
  activity_overlap_comb <- lapply(activity_overlap_comb, function(x) {
    y <- do.call(rbind, x)
    y
  })
  activity_overlap_comb <- do.call(rbind, activity_overlap_comb)

  activity_overlap_comb_genes <- lapply(peak_list, function(x) {
    y <- lapply(x, function(z) {
      return(z[["gene_name"]])
    })
    y
  })
  activity_overlap_comb_genes <- lapply(activity_overlap_comb_genes, function(x) {
    y <- Reduce(append, x)
    y
  })
  activity_overlap_comb_genes <- Reduce(append, activity_overlap_comb_genes)

  genes_overlap <- c(as.character(gene_peaks[-index_overlap,4]), activity_overlap_comb_genes)
  activity_overlap_comb <- Matrix(as.matrix(activity_overlap_comb), sparse = T)
  activity_gene_add <- rbind(peak_matrix[-index_overlap,],activity_overlap_comb)
  rownames(activity_gene_add) <- genes_overlap
  cat("aggregate counts for peaks on the same gene", "\n")
  activity_gene_add_sum <- aggregate.Matrix(activity_gene_add, as.factor(rownames(activity_gene_add)), fun = "sum")
  end_time2 <- Sys.time()
  cat(paste("overall computing", "time", difftime(end_time2, start_time2, units="secs"), "s", "\n", sep = " "))
  return(activity_gene_add_sum)
}

#' @title Annotates peaks to closest gene.
#' @description This function takes the output of \code{peak_on_gene} and annotates peaks which do not fall onto genes to its closest proximal as well as closest distal genes(if closest gene is not proximal). The distances to closest proximal and distal genes is also calculated.
#' @param gene_peaks data.frame of peak annotations. Output of \code{peak_on_gene} with peak annotated to the genes they fall onto.
#' @param annotations list. Output from get \code{get_annotation}.
#' @param gene_element data.frame of gene annotations. Either second slot of \code{get_annotation}, output of \code{too_large}, \code{prolong_upstream} or similar annotation with second column representing start and third column representing end of a peak.
#' @return data.frame with a row for each peak, containing chromosome information, start- and end position and the closest proximal gene and if the closest gene is not proximal the closest distal gene as well. Distances to closest (proximal) gene is included. And given input is the output of \code{peak_on_gene}, if peak is falling on gene, this gene represents closest gene.
#' @examples
#' closest_gene <- peaks_closest_gene(peak_genes,anno_prolong)
#' @export
peaks_closest_gene <- function(peaks, annotations=NULL, gene_element=NULL) {
  start_time2 <- Sys.time()
  if ( is.null(annotations) && is.null(gene_element)) {stop("either annotations or gene_element has to be provided")}
  if ( !is.null(annotations) && !is.null(gene_element)) {stop("either annotations OR gene_ement has to be provided")}
  if (!is.null(annotations)) {
    data <- annotations[[2]]
  }
  else{
    data <- gene_element
  }
  chromosomes <- unique(as.character(data$seqnames))
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
  peaks_annotated <- peaks[peaks[,4] != "nomatch",]
  peaks_annotated <- cbind(peaks_annotated, closest_proximal_gene=peaks_annotated[,4] ,closest_gene=peaks_annotated[,4])
  peaks_not_annotated <- peaks[peaks[,4] == "nomatch",]

  cat("start_looking_for_closest_gene", "\n")
  cores <- as.numeric(future::availableCores() -2)
  cat("available_cores:", cores, "\n")

  peak_list <- list()
  computing <- cores*1000
  data_iterat <- nrow(peaks_not_annotated)%/%computing
  for (n in 1:data_iterat) {
    peak_list[[n]] <- peaks_not_annotated[(n*computing-(computing-1)):(n*computing),]
  }
  if ((nrow(peaks_not_annotated) %% computing) != 0 ) {
    peak_list[[length(peak_list) +1 ]] <- peaks_not_annotated[(length(peak_list)*computing+1):nrow(peaks_not_annotated),]}

  plan(multisession,workers = cores )
  for (n in 1:length(peak_list)) {

    cat("processing_rows " , n*computing-(computing-1), " to ", n*computing, " ")
    start_time <- Sys.time()
    peak_list[[n]] <- future_lapply(seq_along(1:nrow(peak_list[[n]])), function(x, peak,chr_index) {

      chr <- peak[x,1]
      gene_starts <- chr_index[[chr]][["starts"]]
      gene_ends <- chr_index[[chr]][["ends"]]

      left <- gene_starts - as.numeric(peak[x,3]) ### proximal +
      right <- as.numeric(peak[x,2]) - gene_ends ### proximal -
      abs_left <- abs(left) ### shortest distance to gene start
      abs_right <- abs(right) ### shortest distance to gene ends
      left[left < 0] <- Inf  ## if < 0 peak starts after genes start
      right[right < 0] <- Inf ## if < 0 peak starts before genes end
      strand <- chr_index[[chr]][["strand"]]
      left[strand == "-"] <- Inf  ### peak before gene start but genes on the - strand => not proximal
      right[strand == "+"] <- Inf ### peak after gene end but genes on the + strand => not proximal
      if (min(left) < min(right)) {
        gene_proximal <- paste(paste(chr_index[[chr]][["gene_names"]][which(left == min(left))], collapse = "|"),min(left), sep = "_")
      }
      else if (min(left) > min(right)) {
        gene_proximal <- paste(paste(chr_index[[chr]][["gene_names"]][which(right == min(right))], collapse = "|"), min(right), sep = "_")
      }
      else{
        gene_proximal <- "no_proximal"
      }
      if (min(abs_left) < min(abs_right)) {
        gene_general <- paste(paste(chr_index[[chr]][["gene_names"]][which(abs_left == min(abs_left))], collapse = "|"), min(abs_left), sep = "_")
      }
      else{
        gene_general <- paste(paste(chr_index[[chr]][["gene_names"]][which(abs_right == min(abs_right))], collapse = "|"), min(abs_right), sep = "_")
      }
      return(c(peak[x,], gene_proximal, gene_general))

    },chr_index=gene_chrom_index, peak=peak_list[[n]])
    end_time <- Sys.time()
    cat(paste("done", "time", difftime(end_time, start_time, units="secs"), "s", "\n", sep = " "))
  }

  peak_list2 <- lapply(peak_list, function(x) {
    y <- do.call(rbind,x)
  })
  peaks_close <- do.call(rbind,peak_list2)
  rownames(peaks_close) <- rownames(peaks_not_annotated)
  peaks_close <- rbind(peaks_annotated, peaks_close)
  peaks_close <- peaks_close[rownames(peaks),]
  end_time2 <- Sys.time()
  cat(paste("overall computing", "time", difftime(end_time2, start_time2, units="secs"), "s", "\n", sep = " "))
  return(peaks_close)
}

#' @title Get combined overlapping peak set.
#' @description This function takes the peaks which have been called on different samples, and thus could represent overlapping peaks with different start and or end positions to provide an overlapping peaks set.
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
  cat("reduce peaks of different samples to comman overapping peak set")
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
#' @return If \code{do.aggregate} = \code{F}. data.frame with a row for each peak, containing chromosome information, start- and end position and annotated peak with combined peak set. If \code{do.aggregate} = \code{T} and \code{peak_matrix} provided:
#'   \item{agg.comb.peaks}{aggregated peak-cell count matrix across the common overlapping peak set.}
#'   \item{peaks.combined}{data.frame with a row for each peak, containing chromosome information, start- and end position and annotated peak with combined peak set.}
#' @examples
#' overlapped_peaks <- peak_overlap(rownames(merged_atac_filt), combined.peaks=combined.peaks)
#' @export

peak_overlap <- function(peak_features, combined.peaks, do.aggregate=F,peak_matrix=NULL ) {
  if (do.aggregate==T) {
    if(is.null(peak_matrix)) { stop("if do.aggregate=T, peak_matrix has to be provided")}
  }

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
  if( class(as.numeric(peaks[[1]][2])) != "numeric" ) {stop("peak_features need to contain numeric start and end")}
  if( class(as.numeric(peaks[[1]][3])) != "numeric" ) {stop("peak_features need to contain numeric start and end")}
  peaks <- do.call(rbind, peaks)

  cores <- as.numeric(future::availableCores() -2)
  cat("available_cores:", cores, "\n")
  computing <- cores*1000
  data_iterat <- nrow(peaks)%/%computing
  peak_list <- list()

  for (n in 1:data_iterat) {
    peak_list[[n]] <- peaks[(n*computing-(computing-1)):(n*computing),]
  }
  if ((nrow(peaks) %% computing) != 0 ) {
    peak_list[[length(peak_list) +1 ]] <- peaks[(length(peak_list)*computing+1):nrow(peaks),]}
  plan(multisession,workers = cores )
  for (n in 1:length(peak_list)) {
    cat("get overlapping peaks", "\n")
    cat("processing_rows " , n*computing-(computing-1), " to ", n*computing, " ")
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
    end_time <- Sys.time()
    cat(paste("done", "time", difftime(end_time, start_time, units="secs"), "s", "\n", sep = " "))
  }
  peak_list2 <- lapply(peak_list, function(x) { y <- do.call(rbind, x) })
  peaks_combined <- do.call(rbind, peak_list2)
  rownames(peaks_combined) <- peak_features
  colnames(peaks_combined) <- c("seqnames", "start", "end", "combined_peaks")
  end_time2 <- Sys.time()
  cat(paste("overall computing", "time", difftime(end_time2, start_time2, units="secs"), "s", "\n", sep = " "))
  if (do.aggregate==T) {
    if (class(peak_matrix) != "dgTMatrix" ) {
      peak_matrix <- Matrix(as.matrix(peak_matrix),  sparse = T)
    }
    combined_peaks_matrix <- aggregate.Matrix(peak_matrix, as.factor(peaks_combined[,4]), fun = "sum")
    peaks_list3 <- list(agg.comb.peaks=combined_peaks_matrix, peaks.combined=peaks_combined)
    return(peaks_list3)
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

