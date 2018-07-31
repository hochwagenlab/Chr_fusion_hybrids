# This file contains paths to all data files. It is used by the other files as #
# the centralized source of required file paths.                               #
# File locations need to be updated in order to point to the actual files.     #

library(here)
library(stringr)

#------------------------------------------------------------------------------#
#                            MACS2 fragment pileup                             #
#                               bedGraph files                                 #
#------------------------------------------------------------------------------#
base_dir <- '/Volumes/LabShare/HTGenomics/HiSeqOutputs'
replicate_dir <- file.path(base_dir, 'AveReps_SK1_S288c_Yue_hybrid_MACS2_FE')
new_strain_dir <- file.path(base_dir, '2018-03-06_Chr_fusion_hybrids_JefBoeke')

MACS2_pileups <- list(
  WT=file.path(
    # replicate_dir, 'Red1-AH9029_reps_SK1_S288C_Yue_PM_SPMR',
    # 'Red1-AH9029_reps_SK1_S288C_Yue_PM_SPMR_FE.bdg'),
    new_strain_dir, 'AH9029-030618_SK1_S288C_Yue_PM_SPMR',
    'AH9029-030618_SK1_S288C_Yue_PM_SPMR_FE.bdg'),
  'chrIV-I(cen1∆)'=file.path(
    # replicate_dir, 'Red1-AH9033_AH9763_reps_SK1_S288C_Yue_PM_SPMR',
    # 'Red1-AH9033_AH9763_reps_SK1_S288C_Yue_PM_SPMR_FE.bdg'),
    new_strain_dir, 'AH9763-030618_SK1_S288C_Yue_PM_SPMR',
    'AH9763-030618_SK1_S288C_Yue_PM_SPMR_FE.bdg'),
  'chrIV(cen4∆)-I'=file.path(
    new_strain_dir, 'AH9764-030618_SK1_S288C_Yue_PM_SPMR',
    'AH9764-030618_SK1_S288C_Yue_PM_SPMR_FE.bdg'),
  'chrIX(cen9∆)-III-I(cen1∆)'=file.path(
    # replicate_dir, 'Red1-AH9031_AH9765_reps_SK1_S288C_Yue_PM_SPMR',
    # 'Red1-AH9031_AH9765_reps_SK1_S288C_Yue_PM_SPMR_FE.bdg'),
    new_strain_dir, 'AH9765-030618_SK1_S288C_Yue_PM_SPMR',
    'AH9765-030618_SK1_S288C_Yue_PM_SPMR_FE.bdg'),
  'chrIX-III(cen3∆)-I(cen1∆)'=file.path(
    new_strain_dir, 'AH9766-030618_SK1_S288C_Yue_PM_SPMR',
    'AH9766-030618_SK1_S288C_Yue_PM_SPMR_FE.bdg')
)


#------------------------------------------------------------------------------#
#                                   MACS2                                      #
#                                 peak files                                   #
#------------------------------------------------------------------------------#
base_dir <- '/Volumes/LabShare/HTGenomics/HiSeqOutputs'
replicate_dir <- file.path(base_dir, 'AveReps_SK1_S288c_Yue_hybrid_MACS2_FE')

MACS2_peaks <- list(
  'WT'=file.path(
    # replicate_dir, 'Red1-AH9029_reps_SK1_S288C_Yue_PM_SPMR',
    # 'Red1-AH9029_reps_SK1_S288C_Yue_PM_SPMR_peaks.narrowPeak'),
    new_strain_dir, 'AH9029-030618_SK1_S288C_Yue_PM_SPMR',
    'AH9029-030618_SK1_S288C_Yue_PM_SPMR_peaks.narrowPeak'),
  'chrIV-I(cen1∆)'=file.path(
    # replicate_dir, 'Red1-AH9033_AH9763_reps_SK1_S288C_Yue_PM_SPMR',
    # 'Red1-AH9033_AH9763_reps_SK1_S288C_Yue_PM_SPMR_peaks.narrowPeak'),
    new_strain_dir, 'AH9763-030618_SK1_S288C_Yue_PM_SPMR',
    'AH9763-030618_SK1_S288C_Yue_PM_SPMR_peaks.narrowPeak'),
  'chrIV(cen4∆)-I'=file.path(
    new_strain_dir, 'AH9764-030618_SK1_S288C_Yue_PM_SPMR',
    'AH9764-030618_SK1_S288C_Yue_PM_SPMR_peaks.narrowPeak'),
  'chrIX(cen9∆)-III-I(cen1∆)'=file.path(
    # replicate_dir, 'Red1-AH9031_AH9765_reps_SK1_S288C_Yue_PM_SPMR',
    # 'Red1-AH9031_AH9765_reps_SK1_S288C_Yue_PM_SPMR_peaks.narrowPeak'),
    new_strain_dir, 'AH9765-030618_SK1_S288C_Yue_PM_SPMR',
    'AH9765-030618_SK1_S288C_Yue_PM_SPMR_peaks.narrowPeak'),
  'chrIX-III(cen3∆)-I(cen1∆)'=file.path(
    new_strain_dir, 'AH9766-030618_SK1_S288C_Yue_PM_SPMR',
    'AH9766-030618_SK1_S288C_Yue_PM_SPMR_peaks.narrowPeak')
)


#------------------------------------------------------------------------------#
#                                     GFF                                      #
#                                    files                                     #
#------------------------------------------------------------------------------#
base_dir <- file.path('/Volumes/LabShare/GenomeSequences/',
                      'S288C_SK1_Yue_hybrid_genome')
gff_files <- list(
  SK1=file.path(base_dir, 'SK1.all_feature.gff'),
  S288C=file.path(base_dir, 'S288c.all_feature.gff')
)