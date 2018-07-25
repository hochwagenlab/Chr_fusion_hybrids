#------------------------------------------------------------------------------#
#                                                                              #
#                                   Figure X                                   #
#                                                                              #
#------------------------------------------------------------------------------#

library(here)
library(wesanderson)
library(GenomicRanges)

# Load paths to data files
source(here('helper_data_file_locations.R'))
# Load IO code
source(here('helper_io.R'))
# Load analysis functions
source(here('helper_analysis.R'))
# Load ggplot2
source(here('helper_ggplot2_settings.R'))

# Colors -----------------------------------------------------------------#
# wes_colors <- wes_palette(name = "Darjeeling1", n = 4, type = 'continuous')
wes_colors <- wes_palette(name = "Moonrise1", n = 4, type = 'continuous')
strain_colors <- c('black', 'lightblue', wes_colors[2:4])
strain_colors <- c('black', 'lightblue', 'darkgoldenrod1', 'darkolivegreen3', 'indianred1')
names(strain_colors) <- c('WT', 'chrIV-I(cen1∆)', 'chrIV(cen4∆)-I',
                          'chrIX(cen9∆)-III-I(cen1∆)',
                          'chrIX-III(cen3∆)-I(cen1∆)')
chr_colors <- c(Others='black', 'chr I'='mediumorchid3', 'chr III'='green3',
                'chr IV'='maroon3', 'chr IX'='yellow2')

#------------------------------------------------------------------------------#
#                                   Panel C                                    #
# Average signal per chromosome                                                #
#------------------------------------------------------------------------------#
# Load signal track data
signal_track_data <- lapply(MACS2_pileups, import_bedGraph)

# Normalize by average signal on SK1 chromosomes
normalized_data <- lapply(signal_track_data, norm_by_SK1_genome_average)

# Compute average signal per bp on each chromosome
average_signal_per_chr <- lapply(normalized_data,
                                 function(x) average_chr_signal(x)[[1]])

# Add chromosome sizes
get_chr_len <- function(genome) {
  len <- get_chr_coordinates(genome)
  len <- as.data.frame(GenomeInfoDb::seqlengths(len))
  len <- tibble::rownames_to_column(len)
  colnames(len) <- c('chr', 'len')
  return(len)
}

chr_len <- get_chr_len('SK1_S288C')

average_signal_per_chr <- lapply(average_signal_per_chr,
                                 function(x) dplyr::left_join(x, chr_len))

# Add strain name
for (i in seq_along(average_signal_per_chr)) {
  average_signal_per_chr[[i]]$strain <- names(average_signal_per_chr)[i]
}

# Combine data frames
average_signal_per_chr <- do.call('rbind', average_signal_per_chr)

# Adjust chr lengths
chrIV_I_len <- subset(average_signal_per_chr,
                      chr == 'chrIV_S288C', select=len)[1, ] +
               subset(average_signal_per_chr,
                      chr == 'chrI_S288C', select=len)[1, ]
chrIX_III_I_len <- subset(average_signal_per_chr,
                          chr == 'chrIX_S288C', select=len)[1, ] +
                   subset(average_signal_per_chr,
                          chr == 'chrIII_S288C', select=len)[1, ] +
                   subset(average_signal_per_chr,
                          chr == 'chrI_S288C', select=len)[1, ]

average_signal_per_chr[
  (average_signal_per_chr$strain == 'chrIV-I(cen1∆)'
   | average_signal_per_chr$strain == 'chrIV(cen4∆)-I')
   & average_signal_per_chr$chr == 'chrI_S288C', 'len'] <- chrIV_I_len
average_signal_per_chr[
  (average_signal_per_chr$strain == 'chrIV-I(cen1∆)'
   | average_signal_per_chr$strain == 'chrIV(cen4∆)-I') 
   & average_signal_per_chr$chr == 'chrIV_S288C', 'len'] <- chrIV_I_len
average_signal_per_chr[
  (average_signal_per_chr$strain == 'chrIX(cen9∆)-III-I(cen1∆)' 
   | average_signal_per_chr$strain == 'chrIX-III(cen3∆)-I(cen1∆)') 
   & average_signal_per_chr$chr == 'chrIX_S288C', 'len'] <- chrIX_III_I_len
average_signal_per_chr[
  (average_signal_per_chr$strain == 'chrIX(cen9∆)-III-I(cen1∆)'
   | average_signal_per_chr$strain == 'chrIX-III(cen3∆)-I(cen1∆)')
   & average_signal_per_chr$chr == 'chrIII_S288C', 'len'] <- chrIX_III_I_len
average_signal_per_chr[
  (average_signal_per_chr$strain == 'chrIX(cen9∆)-III-I(cen1∆)' 
   | average_signal_per_chr$strain == 'chrIX-III(cen3∆)-I(cen1∆)')
   & average_signal_per_chr$chr == 'chrI_S288C', 'len'] <- chrIX_III_I_len

# Add color-specifying variable
average_signal_per_chr$chr_color <- 'Others'

average_signal_per_chr[
  average_signal_per_chr$chr == 'chrI_S288C', 'chr_color'] <- 'chr I'
average_signal_per_chr[
  average_signal_per_chr$chr == 'chrIII_S288C', 'chr_color'] <- 'chr III'
average_signal_per_chr[
  average_signal_per_chr$chr == 'chrIV_S288C', 'chr_color'] <- 'chr IV'
average_signal_per_chr[
  average_signal_per_chr$chr == 'chrIX_S288C', 'chr_color'] <- 'chr IX'

# Plot (only S288C chrs)
plot_data <- subset(average_signal_per_chr, str_detect(chr, 'S288C'))

# Add shape-specifying variable
plot_data$chr_shape <- 'Native'

chr_subset <- plot_data$strain %in% c('chrIV(cen4∆)-I', 'chrIV-I(cen1∆)') &
  plot_data$chr %in% c('chrI_S288C', 'chrIV_S288C')
plot_data$chr_shape[chr_subset] <- 'Fused'

chr_subset <- plot_data$strain %in% c('chrIX(cen9∆)-III-I(cen1∆)',
                                      'chrIX-III(cen3∆)-I(cen1∆)') &
  plot_data$chr %in% c('chrI_S288C', 'chrIII_S288C', 'chrIX_S288C')
plot_data$chr_shape[chr_subset] <- 'Fused'

# Order strains
plot_data$strain <- factor(plot_data$strain, levels=c(
  'WT', 'chrIV(cen4∆)-I', 'chrIV-I(cen1∆)',
  'chrIX-III(cen3∆)-I(cen1∆)', 'chrIX(cen9∆)-III-I(cen1∆)'))

ggplot(plot_data,
       aes(len / 10^6, avrg_signal, colour = chr_color, shape = chr_shape)) +
  geom_hline(yintercept = 1, lty = 3) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5)) +
  geom_point(stat = "identity", size = 1.8, alpha = 1, stroke=1) +
  facet_wrap(~ strain, ncol = 5) +
  scale_color_manual('', values = chr_colors) +
  scale_shape_manual('', values = c(19, 1)) +
  scale_y_continuous(limits = c(0.5, 1.75)) +
  labs(title = '', x = 'Chromosome size (Mb)',
       y = 'Average Red1\noccupancy') +
  theme(strip.background = element_blank(),
        # panel.border = element_rect(colour = "black"),
        panel.spacing = unit(0.5, "lines"),
        legend.key.height=unit(0.75, "line"))


#------------------------------------------------------------------------------#
#                             Supplementary figure                             #
#                                   Panel B                                    #
#------------------------------------------------------------------------------#

# Plot (only SK1 chrs)
average_signal_per_chr_SK1 <- subset(average_signal_per_chr,
                                     str_detect(chr, 'SK1'))

# Order strains
average_signal_per_chr_SK1$strain <- factor(
  average_signal_per_chr_SK1$strain, levels=c('WT', 'chrIV(cen4∆)-I',
                                              'chrIV-I(cen1∆)',
                                              'chrIX-III(cen3∆)-I(cen1∆)',
                                              'chrIX(cen9∆)-III-I(cen1∆)'))

ggplot(average_signal_per_chr_SK1, aes(len / 10^6, avrg_signal)) +
  geom_hline(yintercept = 1, lty = 3) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0, 1.5)) +
  geom_point(stat = "identity", colour = 'black', size = 2, alpha = 0.8) +
  facet_wrap(~ strain, ncol = 5) +
  # scale_color_manual('', values = chr_colors) +
  scale_y_continuous(limits = c(0.5, 1.75)) +
  labs(title = '', x = 'Chromosome size (Mb)',
       y = 'Average Red1\noccupancy') +
  theme(strip.background = element_blank(),
        # panel.border = element_rect(colour = "black"),
        panel.spacing = unit(1, "lines"))


#------------------------------------------------------------------------------#
#                                   Panel A                                    #
# Signal along chromosomes                                                     #
#------------------------------------------------------------------------------#

# Add seqlengths to data
# Add seqlengths
add_seqlengths <- function(gr, ref_genome='SK1_S288C') {
  message(deparse(substitute(gr)), ': ', appendLF=F)
  
  if (ref_genome %in% c('SK1', 'S288C', 'SK1_S288C')) {
    gr_lens <- get_chr_coordinates(ref_genome)
  } else stop('"ref_genome" must be a valid input to ',
              'the "genome" argument of "get_chr_coordinates()".')
  # Sort sequences and levels to make sure they match
  gr <- sort(sortSeqlevels(gr))
  gr_lens <- sortSeqlevels(gr_lens)
  
  # Add seqlengths to signal object
  seqlengths(gr) <- seqlengths(gr_lens)
  message('Done!')
  gr
}

normalized_data <- lapply(
  normalized_data, add_seqlengths, ref_genome='SK1_S288C')

get_signal_to_plot <- function(gr, sample_name, chr_name, window=100) {
  # Get selected chr
  gr <- keepSeqlevels(gr, chr_name, pruning.mode="coarse")
  
  # Compute tiling windows
  bins <- tileGenome(GenomeInfoDb::seqlengths(gr), tilewidth=window,
                     cut.last.tile.in.chrom=TRUE)
  
  # Get signal as "RleList";
  score <- coverage(gr, weight="score")
  
  # Compute average signal per tile
  bins <- binnedAverage(bins, score, "binned_score")
  
  # Get positions as the midpoints of the intervals
  positions <- bins@ranges@start + floor(bins@ranges@width / 2)
  
  # Make data frame (convert positions to Kb; signal is the binned score)
  data.frame(position=positions / 1000, signal=bins$binned_score, strain=sample_name)
}


cen_midpoint <- function(chr, ref_genome='SK1_S288C', y_coord=0) {
  cen <- get_chr_coordinates(genome=ref_genome)
  
  if (!missing(chr)) {
    # Keep only required chr
    cen <- cen[cen@seqnames == chr]
    message('   Drop all chromosomes except ', chr, '...')
    gr <- keepSeqlevels(cen, chr, pruning.mode="coarse")
  }
  
  start <- cen@ranges@start
  half_width <- cen@ranges@width / 2
  cen_mid <- round(start + half_width)
  
  data.frame(seqnames = cen@seqnames, cen_mid = cen_mid, y = y_coord)
}


plot_signal_along_chr <- function(chr='chrI_S288C', chr_name='chr I', x1, x2,
                                  ref_genome='SK1_S288C', sample_names,
                                  x2_annotation_coords=c(200, 12),
                                  compress_window=50, colors, plot_cen=TRUE,
                                  annotate_strains=FALSE) {
  
  if (!missing(x2)) {
    if (length(sample_names) != 2) {
      stop('Length of "sample_names" argument input must be 2 ',
           'when both "x1" and "x2" are provided.', call.=FALSE)
    }
    if (length(colors) != 2) {
      stop('Length of "colors" argument input must be 2 ',
           'when both "x1" and "x2" are provided.', call.=FALSE)
    }
  }
  
  # Get signal data
  x1_signal <- get_signal_to_plot(x1, sample_name=sample_names[1],
                                  chr_name=chr, window=compress_window)
  
  if (!missing(x2)) {
    x2_signal <- get_signal_to_plot(x2, sample_name=sample_names[2],
                                    chr_name=chr, window=compress_window)
    plot_title <- paste(sample_names[1], ' + ', sample_names[2])
  } else plot_title <- sample_names[1]
  
  p <- ggplot(x1_signal, aes(position, signal, colour=strain)) +
    scale_color_manual('', values=colors, guide=FALSE) +
    scale_fill_manual('', values=colors, guide=FALSE) +
    geom_area(position='identity', aes(fill=strain, colour=strain),
              alpha=1, size=0.25) +
    scale_y_continuous(limits = c(-0.5, 12)) +
    labs(title = '', x = paste0("Position on ", chr_name, " (Kb)"),
         y = "Red1 occupancy") +
    theme(strip.background = element_blank())
  
  if (!missing(x2)) {
    p <- p +  geom_area(data=x2_signal, position='identity',
                        aes(fill=strain, colour=strain),
                        alpha=1, size=0.25)
  }
  
  if (plot_cen) {
    cen_mid <- cen_midpoint(chr=chr, ref_genome='SK1_S288C', y_coord=0)
    p <- p + geom_point(data = cen_mid, aes(cen_mid/1000, 0.5),
                        size = 1.5, colour = 'black') +
      annotate("text", x=cen_mid$cen_mid / 1000, y=-0.5,
               label = "CEN", size=2, colour='black')
  }
  
  if (annotate_strains) {
    p <- p + annotate('text', x=25, y=12, label=sample_names[1],
                      size=3, colour=colors[1])
    if (!missing(x2)) {
      p <- p + annotate('text', x=x2_annotation_coords[1],
                        y=x2_annotation_coords[2], label=sample_names[2],
                        size=3, colour=colors[2])
    }
  }
  
  return(p)
}

plot_signal_along_chr(
  chr='chrIX_S288C', chr_name='chr IX',
  x1=normalized_data[['WT']], x2=normalized_data[['chrIX(cen9∆)-III-I(cen1∆)']],
  sample_names=c('WT', 'chrIX(cen9∆)-III-I(cen1∆)'), compress_window=100,
  colors=strain_colors[c('WT', 'chrIX(cen9∆)-III-I(cen1∆)')],
  plot_cen=TRUE, annotate_strains=TRUE)

plot_signal_along_chr(
  chr='chrIX_S288C', chr_name='chr IX',
  x1=normalized_data[['WT']], x2=normalized_data[['chrIX-III(cen3∆)-I(cen1∆)']],
  sample_names=c('WT', 'chrIX-III(cen3∆)-I(cen1∆)'), compress_window=100,
  colors=strain_colors[c('WT', 'chrIX-III(cen3∆)-I(cen1∆)')],
  plot_cen=TRUE, annotate_strains=TRUE)


#------------------------------------------------------------------------------#
#                             Supplementary figure                             #
#                                   Panel A                                    #
#------------------------------------------------------------------------------#
plot_signal_along_chr(
  chr='chrI_S288C', chr_name='chr I',
  x1=normalized_data[['WT']], x2=normalized_data[['chrIV(cen4∆)-I']],
  sample_names=c('WT', 'chrIV(cen4∆)-I'), x2_annotation_coords=c(120, 12),
  compress_window=50, colors=strain_colors[c('WT', 'chrIV(cen4∆)-I')],
  plot_cen=TRUE, annotate_strains=TRUE)

plot_signal_along_chr(
  chr='chrI_S288C', chr_name='chr I',
  x1=normalized_data[['WT']], x2=normalized_data[['chrIV-I(cen1∆)']],
  sample_names=c('WT', 'chrIV-I(cen1∆)'), x2_annotation_coords=c(120, 12),
  compress_window=50, colors=strain_colors[c('WT', 'chrIV-I(cen1∆)')],
  plot_cen=TRUE, annotate_strains=TRUE)

plot_signal_along_chr(
  chr='chrI_S288C', chr_name='chr I',
  x1=normalized_data[['WT']], x2=normalized_data[['chrIX(cen9∆)-III-I(cen1∆)']],
  sample_names=c('WT', 'chrIX(cen9∆)-III-I(cen1∆)'),
  x2_annotation_coords=c(120, 13), compress_window=50,
  colors=strain_colors[c('WT', 'chrIX(cen9∆)-III-I(cen1∆)')],
  plot_cen=TRUE, annotate_strains=TRUE)

plot_signal_along_chr(
  chr='chrI_S288C', chr_name='chr I',
  x1=normalized_data[['WT']], x2=normalized_data[['chrIX-III(cen3∆)-I(cen1∆)']],
  sample_names=c('WT', 'chrIX-III(cen3∆)-I(cen1∆)'),
  x2_annotation_coords=c(120, 13), compress_window=50,
  colors=strain_colors[c('WT', 'chrIX-III(cen3∆)-I(cen1∆)')],
  plot_cen=TRUE, annotate_strains=TRUE)

plot_signal_along_chr(
  chr='chrIII_S288C', chr_name='chr III',
  x1=normalized_data[['WT']], x2=normalized_data[['chrIX(cen9∆)-III-I(cen1∆)']],
  sample_names=c('WT', 'chrIX(cen9∆)-III-I(cen1∆)'), compress_window=75,
  colors=strain_colors[c('WT', 'chrIX(cen9∆)-III-I(cen1∆)')],
  plot_cen=TRUE, annotate_strains=TRUE)

plot_signal_along_chr(
  chr='chrIII_S288C', chr_name='chr III',
  x1=normalized_data[['WT']], x2=normalized_data[['chrIX-III(cen3∆)-I(cen1∆)']],
  sample_names=c('WT', 'chrIX-III(cen3∆)-I(cen1∆)'), compress_window=75,
  colors=strain_colors[c('WT', 'chrIX-III(cen3∆)-I(cen1∆)')],
  plot_cen=TRUE, annotate_strains=TRUE)

plot_signal_along_chr(
  chr='chrIV_S288C', chr_name='chr IV',
  x1=normalized_data[['WT']], x2=normalized_data[['chrIV(cen4∆)-I']],
  sample_names=c('WT', 'chrIV(cen4∆)-I'), compress_window=350,
  colors=strain_colors[c('WT', 'chrIV(cen4∆)-I')],
  plot_cen=TRUE, annotate_strains=TRUE)

plot_signal_along_chr(
  chr='chrIV_S288C', chr_name='chr IV',
  x1=normalized_data[['WT']], x2=normalized_data[['chrIV-I(cen1∆)']],
  sample_names=c('WT', 'chrIV-I(cen1∆)'), compress_window=350,
  colors=strain_colors[c('WT', 'chrIV-I(cen1∆)')],
  plot_cen=TRUE, annotate_strains=TRUE)


#------------------------------------------------------------------------------#
#                                   Panel C                                    #
# Signal in peals along chromosomes                                            #
#------------------------------------------------------------------------------#
# Import peak files
peaks <- lapply(MACS2_peaks, import_MACS2_peaks)

# Calculate average score per peak and put in tidy table to plot
mean_score_per_peak <- function(gr, peaks, chr='chrI_S288C',
                                ref_genome='SK1_S288C') {
  message('Running ', deparse(substitute(gr)), '...')
  genome_info <- get_chr_coordinates(genome = ref_genome)
  
  # Sort sequences and levels to make sure they match
  gr <- sort(sortSeqlevels(gr))
  peaks <- sort(sortSeqlevels(peaks))
  genome_info <- sortSeqlevels(genome_info)
  
  # Add info to signal object
  seqlengths(gr) <- seqlengths(genome_info)
  seqlengths(peaks) <- seqlengths(genome_info)
  
  gr_chr <- keepSeqlevels(gr, chr, pruning.mode="coarse")
  peaks_chr <- keepSeqlevels(peaks, chr, pruning.mode="coarse")
  
  # Get signal as "RleList"; the signal is stored in the "SK1_norm_score" metadata column
  score <- coverage(gr_chr, weight="score")
  
  # Compute and return signal per tile
  bins <- binnedAverage(peaks_chr, score, "binned_score")
  message('   Done!')
  
  return(bins$binned_score)
}

get_chr_quantiles <- function(chr='chrI_S288C', peaks) {
  chr_quants <- peaks[seqnames(peaks) == chr]
  
  chr_quants$WT <- mean_score_per_peak(
    normalized_data[['WT']], chr_quants, chr=chr)
  chr_quants$chrIV_I_CEN4 <- mean_score_per_peak(
    normalized_data[['chrIV-I(cen1∆)']], chr_quants, chr=chr)
  chr_quants$chrIX_III_I_CEN3 <- mean_score_per_peak(
    normalized_data[['chrIX(cen9∆)-III-I(cen1∆)']], chr_quants, chr=chr)
  chr_quants$chrIV_I_CEN1 <- mean_score_per_peak(
    normalized_data[['chrIV(cen4∆)-I']], chr_quants, chr=chr)
  chr_quants$chrIX_III_I_CEN9 <- mean_score_per_peak(
    normalized_data[['chrIX-III(cen3∆)-I(cen1∆)']], chr_quants, chr=chr)
  
  chr_quants$quant <- cut(chr_quants$WT, breaks=quantile(chr_quants$WT),
                          labels=1:4, include.lowest=TRUE)
  
  chr_quants$chrIV_I_CEN4_ratio <- chr_quants$chrIV_I_CEN4 / chr_quants$WT
  chr_quants$chrIX_III_I_CEN3_ratio <- chr_quants$chrIX_III_I_CEN3 / chr_quants$WT
  chr_quants$chrIV_I_CEN1_ratio <- chr_quants$chrIV_I_CEN1 / chr_quants$WT
  chr_quants$chrIV_I_CEN1_ratio <- chr_quants$chrIV_I_CEN1 / chr_quants$WT
  chr_quants$chrIX_III_I_CEN9_ratio <- chr_quants$chrIX_III_I_CEN9 / chr_quants$WT
  
  chr_quants
}

chrI_quants <- get_chr_quantiles(chr='chrI_S288C', peaks=peaks[['WT']])
chrIII_quants <- get_chr_quantiles(chr='chrIII_S288C', peaks=peaks[['WT']])
chrIX_quants <- get_chr_quantiles(chr='chrIX_S288C', peaks=peaks[['WT']])
chrIV_quants <- get_chr_quantiles(chr='chrIV_S288C', peaks=peaks[['WT']])


plot_signal_in_peak_ratio <- function(chr_quants_gr, chr='chrI_S288C',
                                      chr_name='chr I', trim_chr=c(0, 0),
                                      y_lim=c(-1, 1), fit_span=0.75,
                                      show_legend=FALSE) {
  position <- round(start(chr_quants_gr) + width(chr_quants_gr) / 2)
  to_plot <- data.frame(
    position=position, score_ratio=chr_quants_gr$chrIV_I_CEN1_ratio,
    quant=chr_quants_gr$quant, strain='chrIV(cen4∆)-I')
  to_plot <- rbind(to_plot, data.frame(
    position=position, score_ratio=chr_quants_gr$chrIV_I_CEN4_ratio,
    quant=chr_quants_gr$quant, strain='chrIV-I(cen1∆)'))
  to_plot <- rbind(to_plot, data.frame(
    position=position, score_ratio=chr_quants_gr$chrIX_III_I_CEN3_ratio,
    quant=chr_quants_gr$quant, strain='chrIX(cen9∆)-III-I(cen1∆)'))
  to_plot <- rbind(to_plot, data.frame(
    position=position, score_ratio=chr_quants_gr$chrIX_III_I_CEN9_ratio,
    quant=chr_quants_gr$quant, strain='chrIX-III(cen3∆)-I(cen1∆)'))
  
  # Trim chromosome
  start <- trim_chr[1]
  end <- tail(end(chr_quants_gr), 1) - trim_chr[2]
  to_plot <- to_plot[to_plot$position >= start & to_plot$position <= end, ]
  
  cen_mid <- cen_midpoint(chr=chr, ref_genome='SK1_S288C', y_coord=0)
  
  p <- ggplot(to_plot, aes(x = position / 1000, y = log2(score_ratio),
                           colour = strain)) +
    geom_hline(aes(yintercept = 0), linetype = 3) +
    geom_point(size = 0.75, alpha = 0.1) +
    geom_smooth(span=fit_span, se = FALSE, size = 0.75) +
    geom_point(data = cen_mid, aes(cen_mid/1000, 0),
               size = 2, colour = 'black') +
    annotate("text", x=cen_mid$cen_mid / 1000, y=-0.2, label = "CEN",
             size=2, colour='black') +
    scale_color_manual('', values=strain_colors) +
    ylim(y_lim) +
    labs(title = '', x = paste0("Position on ", chr_name, ' (Kb)'),
         y = "Red1 occupancy\nlog2(sample/wild type)") +
    theme(axis.title.y = element_text(size=11)) 
    
    if (!show_legend) p <- p + guides(colour=FALSE)
  
  p
}

plot_signal_in_peak_ratio(
  chrI_quants, chr='chrI_S288C', chr_name='chr I', trim_chr=c(0, 0),
  y_lim=c(-1.5, 1.5), fit_span=0.75)

plot_signal_in_peak_ratio(
  chrIII_quants, chr='chrIII_S288C', chr_name='chr III', trim_chr=c(0, 0),
  y_lim=c(-1.5, 1.5), fit_span=0.75)

plot_signal_in_peak_ratio(
  chrIX_quants, chr='chrIX_S288C', chr_name='chr IX', trim_chr=c(0, 0),
  y_lim=c(-1.5, 1.5), fit_span=0.5, show_legend=TRUE)

plot_signal_in_peak_ratio(
  chrIV_quants, chr='chrIV_S288C', chr_name='chr IV', trim_chr=c(0, 0),
  y_lim=c(-1.5, 1.5), fit_span=0.25)

