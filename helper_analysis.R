#------------------------------------------------------------------------------#
#                              Analysis helper                                 #
#                                 functions                                    #
#------------------------------------------------------------------------------#

library(here)
library(GenomicRanges)
library(tidyverse)

# Load paths to data files
source(here('helper_data_file_locations.R'))
# Load IO code
source(here('helper_io.R'))


#' Add genome name to GRanges
#'
#' Adds genome name to \code{GRanges} object.
#' @param gr \code{GRanges} object. No default.
#' @param name Genome name \code{string}. No default.
#' @return Same \code{GRanges} object with a genome name.
#' @examples
#' \dontrun{
#' x <- add_genome_name_to_GR(WT, name='SK1_S288CYue')
#' }
#' @export
add_genome_name_to_GR <- function(gr, name) {
  number_seqs <- length(levels(gr@seqnames))
  gr@seqinfo@genome <- rep(name, number_seqs)
  
  gr
}


#' Get chromosome coordinates
#'
#' Returns chromosome length and centromere coordinates for an hybrid of the SK1
#' and S288C assemblies from Yue et al. 2017 (here named "SK1Yue"). Chromosome
#' lengths were calculated in the bash shell using commands like the following 
#' (SK1 example): \cr
#' \code{cat SK1.genome.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100)
#' "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'}
#' @param genome Character object specifying the genome version; accepts one of
#' the following options:
#' \enumerate{
#'   \item \code{"SK1"}
#'   \item \code{"S288C"}
#'   \item \code{"SK1_S288C"}
#' }
#' Defaults to \code{"SK1_S288C"}.
#' @return SK1 coordinates as a \code{GRanges} object.
#' @examples
#' \dontrun{
#' x <- get_chr_coordinates(genome="SK1")
#' }
#' @export

get_chr_coordinates <- function(genome='SK1_S288C') {
  SK1_gff <- rtracklayer::import.gff(gff_files[['SK1']])
  S288C_gff <- rtracklayer::import.gff(gff_files[['S288C']])
  
  SK1_start <- SK1_gff[SK1_gff$type == 'centromere']@ranges@start
  SK1_end <- SK1_start + SK1_gff[SK1_gff$type == 'centromere']@ranges@width - 1
  
  S288C_start <- S288C_gff[S288C_gff$type == 'centromere']@ranges@start
  S288C_end <- 
    S288C_start + S288C_gff[S288C_gff$type == 'centromere']@ranges@width - 1
  
  # Coordinate df
  if (genome == 'SK1') {
    coord_table <- data.frame(
      "Chromosome" = paste0('chr', as.roman(1:16)),
      "Start" = SK1_start, "End" = SK1_end,
      "LenChr" =  c(228861, 829469, 340914, 1486921, 589812, 299318,
                    1080440, 542723, 449612, 753937, 690901, 1054145,
                    923535, 791982, 1053869, 946846))
  } else if (genome == 'S288C') {
    coord_table <- data.frame(
      "Chromosome" = paste0('chr', as.roman(1:16)),
      "Start" = S288C_start, "End" = S288C_end,
      "LenChr" = c(219929, 813597, 341580, 1566853, 583092, 271539,
                   1091538, 581049, 440036, 751611, 666862, 1075542,
                   930506, 777615, 1091343, 954457))
  } else if (genome == 'SK1_S288C') {
    coord_table <- data.frame(
      "Chromosome" = c(paste0('chr', as.roman(1:16), '_SK1'),
                       paste0('chr', as.roman(1:16), '_S288C')),
      "Start" = c(SK1_start, S288C_start),
      "End" = c(SK1_end, S288C_end),
      "LenChr" = c(
        c(228861, 829469, 340914, 1486921, 589812, 299318,
          1080440, 542723, 449612, 753937, 690901, 1054145,
          923535, 791982, 1053869, 946846),
        c(219929, 813597, 341580, 1566853, 583092, 271539,
          1091538, 581049, 440036, 751611, 666862, 1075542,
          930506, 777615, 1091343, 954457)))
  } else stop('"genome" argument must be one of ',
              '"SK1", "S288C" or "SK1_S288C".')
  
  # Convert to GRanges
  SK1_S288C <- with(coord_table,
                    GRanges(Chromosome, IRanges(Start + 1, End),
                            seqlengths=setNames(LenChr, Chromosome)))
  # Add genome name and return
  add_genome_name_to_GR(SK1_S288C, name=genome)
}


#' Average signal genome-wide and on each chromosome
#'
#' Computes average signal (in the \code{score} \code{GRanges} metadata column)
#' genome-wide and on each chromosome (each individual sequence determined by
#' \code{seqnames}).
#' @param gr Input signal track data as a \code{GRanges} object (see
#' \code{?"GRanges-class"} for more details). No default.
#' @return List with two elements:
#' \enumerate{
#'   \item \code{seq_avrg} Average \code{score} on each sequence (dataframe)
#'   \item \code{genome_avrg} Average \code{score} genome-wide (named vector)
#' }
#' @examples
#' \dontrun{
#' average_chr_signal(GRanges_object)
#' }
#' @export
average_chr_signal <- function(gr){
  # IO checks
  if (!is(gr, "GRanges")) stop('input must be a GRanges object.')
  if (!"score" %in% names(mcols(gr))) {
    stop(deparse(substitute(gr)), ' does not have a "score" metadata column.')
  }
  
  message('Computing average signal...')
  avrg <- function(x) (sum(width(x) * score(x), na.rm = T) / sum(width(x)))
  genome_avrg <- avrg(gr)
  seq_avrg <- sapply(split(gr, seqnames(gr)), avrg)
  
  # Convert to dataframe
  seq_avrg <- data.frame(chr=names(seq_avrg), avrg_signal=seq_avrg,
                         row.names=NULL, stringsAsFactors = F)
  
  message('Done!')
  return(list("seq_avrg"=seq_avrg, "genome_avrg"=genome_avrg))
}


#' Normalize by SK1 genome average
#'
#' Computes average signal genome-wide on the SK1 chromosomes and divides every
#' signal value by computed average.
#' @param gr \code{GRanges} object. No default.
#' @return Same \code{GRanges} object with average-normalized signal.
#' @examples
#' \dontrun{
#' x <- norm_by_genome_average(gr=WT)
#' }
#' @export
norm_by_SK1_genome_average <- function(gr){
  # Compute genome average
  message('   SK1 chromosome signal average: ', appendLF=FALSE)
  
  suppressMessages(
    SK1_mean <- average_chr_signal(gr[str_detect(gr@seqnames, 'SK1')])[[2]]
  )
  message(round(SK1_mean, 3))
  
  # Divide each score value by genome average signal
  message('   Normalize signal...')
  gr$score <- gr$score / SK1_mean
  
  message('   ---')
  message('   Done!')
  return(gr)
}