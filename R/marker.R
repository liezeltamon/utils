# From Braun et al. 2023 Science Comprehensive cell atlas of the first-trimester developing human brain
brain_class <- c("CDK1", "INA", "NHLH1", "SOX2", "HES1", "BCAN", "SOX10", "PDGFRA", "PTPRC", "CLDN5", "RGS5", "HBB")

radial_glia_like_ipc <- c("SOX2",
                          "SOX21", "FGFR2", "FOXP1" ## These genes, promoting cell cycle reentry, are among the ones upregulated in radial glia-like subtype of IPCs
                          )
neuronal_like_ipc <- "NEUROD6" ## Increasing expression as neuronal IPC goes through cell cycle till post mitotic cells exiting to become neuroblasts (See Braun et al. 2023 Fig. 3E)

glioblasts <- c("BCAN", "TNC") ## Together identified putative glioblasts in all brain regions

## About half of all glioblasts additionally expressed the astrocyte-specific water channel AQP4 and tight-junction Connexin-43 (encoded by GJA1), which we defined as pre-astrocytes (fig. S12D)
pre_astrocytes <- c("AQP4", "GJA1")

## Emerged later in small numbers compared to glioblasts and pre-astrocytes
pre_opcs <- c("EGFR", "DLL3")

## Patterning genes underlying the regionalisation of radial glia, OPCs, and neurons conserged between human and mouse
forebrain_patterning <- "FOXG1"
mid_hindbrain_boundary_patterning <- "EN1"
midbrain_patterning <- "GATA3"
hindbrain_patterning <- "HOXD3"

## OPC early prototype markers
early_opcs <- "PDGFRA"
early_cops <- c("ENPP6", "MAG") ### COPs - commited OPCs

## Ganglionic eminences markers based on Braun et al. 2023 Fig. S8; rhere are different ganglionic eminences within ventral telecenphalon
lge <- c("MSN",
         "MEIS2" # LGE/CGE marker from Fleck
         )
cge <- c("CALB2", "MEIS2")
mge <- c("SST", "CRABP1")
ge_unclassified <- c("OTP", "LHX1", "FOXD1")
ventral_midbrain <- "FOXA1"
ventral_telencephalon <- c("DLX2", "FOXG1")

# From Pollen et al. 2014

neuron_overRG <- c("STMN2")
RG_overNeuron <- c("SOX2", "HES1")

# From Fleck et al. 2022 Inferring and perturbing cell fate regulomes in human brain organoids

## Based on text and extended data Fig. 9

# POU5F1 - early
# ZIC2 - early before FOXG1 expression
# PAX6, VIM - neural progenitor cell (NPC)
# HES1, GLI3 - dorsoventral branchpoint
# FOXG1 - telencephalon, non-telencephalon (nt), ventral telencephalon, dorsal telencephalon
# TCF7L2, LHX9 - nt
# RSPO1 - earlier nt (based on Fleck Fig.1f) 
# LHX1 - later nt (based on Fleck Fig.1f) 
# SIX6 - telencephalon (before branching)
# DLX5, DLX2 - ventral telencenphalon / ganglionic eminences
# MEIS1 - LGE/CGE ganglionic eminences
# NEUROD6, EMX1 - dorsal
# STMN2, DCX - marker of more mature cell types regardless of branch

branch_markers <- c("POU5F1", "ZIC2", "PAX6", "VIM", "HES1", "GLI3", "STMN2", "DCX", "TCF7L2", "LHX9", "RSPO1", "LHX1", "FOXG1", "SIX6", "TTF1", "DLX2", "DLX5", "MEIS2", "NEUROD6", "EMX1") # NKX2-1 = TTF1

## gRNAs targreting thse TFs showed enrichment of ventral telencenphalon branch + depletion of other regions e.g. cortex
dorsal_over_others <- c("GLI3", "TBR1")

## Opposing effect as targeting GLI3 or TBR1
ventral_over_others <- c("HES1", "HOPX")

## Two genes have opposing effect on dorsal telencephalon commitment
## Genes activated by GLI3 correlate with positive transition to dorsal branch
## HES1 has repressive effect on those genes
dorsoventral_branchpoint <- c("GLI3", "HES1")

## This suggests that SHH promotes ventralization predominantly by preventing GLI3-induced dorsalization39,44. (Fleck et al. 2022 Inferring)
ventralisation <- "SHH"
dorsalisation <- "GLI3"

## TFs expressed across primary and organoid cell types and targeted in Fleck targeted in the single-cell genomic perturbation experiment

fleck_tftargets_across_stages <- c("PAX6", "SOX1", "GLI3", "SOX9", "HES1", "HOPX", "LHX2", "MEIS1", "E2F2", "FOXN4", "NEUROD1", "MYT1L", "BACH2", "ZFPM2", "ST18", "NEUROD6", "SOX5", "ZEB2", "TBR1", "BCL11B")

## Based on Braun 2023 Fig. 3N

radial_glia_and_glioblast_early <- c("LIX1", "HMGA2", "CCND1", "LEF1", "LRRN1", "DDIT4")
radial_glia_and_glioblast_late <- c("TCIM", "HOPX", "FAM107A", "PON2", "FGFBP3", "FOS", "CDO1", "ITGB8", "SLC1A3", "BCAN", "PTN", "FABP7")
neuronal_ipc_early <- c("DLL1", "RGS16", "BNIP3", "CFAP298", "DLL3")
neuronal_ipc_late <- c("WIPF3", "POU3F2", "SMOC1")
neuroblast_early <- c("PLPPR1", "PPP1R17", "NHLH2", "CRABP1", "NXPH4")
neuroblast_late <- c("ZNF704", "ROBO2",
                     "MN1", "DPY19L1", "EPHB6")
neuron_early <- c("CALB2", "CNTNAP2", "ARG2", "CDH7")
neuron_late <- c("SLA", "NTM",
                 "LIMCH1", "CHL1")