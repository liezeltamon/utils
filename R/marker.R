# From Braun et al. 2023 Science Comprehensive cell atlas of the first-trimester developing human brain
brain_class <- c("CDK1", "INA", "NHLH1", "SOX2", "HES1", "BCAN", "SOX10", "PDGFRA", "PTPRC", "CLDN5", "RGS5", "HBB")

radial_glia_like_IPC <- c("SOX2",
                          "SOX21", "FGFR2", "FOXP1" ## These genes, promoting cell cycle reentry, are among the ones upregulated in radial glia-like subtype of IPCs
                          )
neuronal_like_IPC <- "NEUROD6" ## Increasing expression as neuronal IPC goes through cell cycle till post mitotic cells exiting to become neuroblasts (See Braun et al. 2023 Fig. 3E)

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

# From Pollen et al. 2014

neuron_overRG <- c("STMN2")
RG_overNeuron <- c("SOX2", "HES1")

# From Fleck et al. 2022 Inferring and perturbing cell fate regulomes in human brain organoids

## gRNAs targreting thse TFs showed enrichment of ventral telencenphalon branch + depletion of other regions e.g. cortex
dorsal_over_others <- c("GLI3", "TBR1")

## Opposing effect as targeting GLI3 or TBR1
ventral_over_others <- c("HES1", "HOPX")

## Two genes have opposing effect on dorsal telencephalon commitment
## Genes activated by GLI3 correlate with positive transition to dorsal branch
## HES1 has repressive effect on those genes
dorsoventral_branchpoint <- c("GLI3", "HES1")

# This suggests that SHH promotes ventralization predominantly by preventing GLI3-induced dorsalization39,44. (Fleck et al. 2022 Inferring)
ventralisation <- "SHH"
dorsalisation <- "GLI3"
