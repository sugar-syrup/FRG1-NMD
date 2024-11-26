#!/usr/bin/env Rscript


# Title: ISAR_gencode.v33.SIRVomeERCCome
# Objective: Standard IsoformSwitchAnalyseR (ISAR) pipeline giving differential transcript usage (DTU) with gencode annotations and first output analyses
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

######
# Load libraries
######

library(optparse)
library(IsoformSwitchAnalyzeR)

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 2)

# mydir is the first argument
mydir=arguments$args[1]

# condition string is the second, gets converted to character
cond=arguments$args[2]
cond <- unlist(strsplit(cond,","))

# Get design table
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE)
myDesign <- data.frame(
    sampleID = samples$sample,
    condition = samples$condition
)

# Compile sampleVector from experiment file
myFiles <- file.path(mydir, "Salmon", samples$sample, "quant.sf") 
names(myFiles) <- samples$sample

print(myFiles)

# Import quantifications
salmonQuant <- importIsoformExpression(
    sampleVector = myFiles,
    addIsofomIdAsColumn = TRUE,
    readLength = 36,
)

# Get comparison data frame
myComparison <- data.frame(
    condition_1 = "control",
    condition_2 = cond[-1]
)

# Generate switchlist
SwitchList <- importRdata(
		isoformCountMatrix   = salmonQuant$counts,
		isoformRepExpression = salmonQuant$abundance,
		designMatrix         = myDesign,
		addAnnotatedORFs     = TRUE,
		onlyConsiderFullORF = FALSE,
		removeNonConvensionalChr = TRUE,
		isoformExonAnnoation = file.path("/home", "shithij", "projects", "FRG1-NMD","reference", "gencode.v47.annotation.gtf"),
		comparisonsToMake= myComparison,	
		showProgress = TRUE,
		ignoreAfterBar = TRUE,
    		ignoreAfterSpace = TRUE,
    		ignoreAfterPeriod = FALSE,	# Set to FALSE for Gencode
    		removeTECgenes = TRUE		# If set to TRUE, spike_ins need "gene_name" and "gene_type"
)

# SwitchList <- isoformSwitchAnalysisPart1(
#     switchAnalyzeRlist   = SwitchList,
#     # pathToOutput = 'path/to/where/output/should/be/'
#     outputSequences      = FALSE, # change to TRUE whan analyzing your own data 
#     prepareForWebServers = FALSE  # change to TRUE if you will use webservers for external sequence analysis
# )

# extractSwitchSummary(SwitchList)

# SwitchList <- isoformSwitchAnalysisPart2(
#   switchAnalyzeRlist        = SwitchList, 
#   dIFcutoff                 = 0.3,   # Cutoff for defining switch size - set high for short runtime in example data
#   n                         = 10,    # if plotting was enabled, it would only output the top 10 switches
#   removeNoncodinORFs        = TRUE,  # Because ORF was predicted de novo
#   pathToCPC2resultFile      = system.file("extdata/cpc2_result.txt"         , package = "IsoformSwitchAnalyzeR"),
#   pathToPFAMresultFile      = system.file("extdata/pfam_results.txt"        , package = "IsoformSwitchAnalyzeR"),
#   pathToIUPred2AresultFile  = system.file("extdata/iupred2a_result.txt.gz"  , package = "IsoformSwitchAnalyzeR"),
#   pathToSignalPresultFile   = system.file("extdata/signalP_results.txt"     , package = "IsoformSwitchAnalyzeR"),
#   outputPlots               = FALSE  # keeps the function from outputting the plots from this example
# )

# What is in the Switchlist
str(SwitchList)
write.csv(SwitchList$isoformFeatures, file = "isoformFeatures.csv")
write.csv(SwitchList$geneFeatures, file = "geneFeatures.csv")

# Filtering
SwitchList_filt <- preFilter(
  SwitchList,
  geneExpressionCutoff = 1, # FPMK threshold
  isoformExpressionCutoff = 0, # FPMK threshold
  IFcutoff=0.01,
  removeSingleIsoformGenes = TRUE,
  reduceToSwitchingGenes=FALSE,
  alpha=0.05,
  dIFcutoff = 0.1,
  quiet=FALSE
)

print("Filtering done")

# Testing for Isoform Switches via DEXSeq
SwitchList_filt_Analyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = SwitchList_filt,
    reduceToSwitchingGenes=TRUE
)

# Set WD
setwd((file.path(mydir, "ISAR")))

save.image(file='ISAR_session.RData')

# Summarize switching features
extractSwitchSummary(SwitchList_filt_Analyzed)

# Predicting Switch Consequences
SwitchList_filt_Analyzed <- analyzeSwitchConsequences(
    SwitchList_filt_Analyzed,
    consequencesToAnalyze = c('NMD_status'), 
    dIFcutoff = 0.1, 
    showProgress=TRUE
)


# Summarize switching features without consequences
extractSwitchSummary(SwitchList_filt_Analyzed ) # , dIFcutoff = 0.1, filterForConsequences = FALSE)

# Summarize switching features with consequences
extractSwitchSummary(SwitchList_filt_Analyzed ) #, dIFcutoff = 0.1, filterForConsequences = TRUE)

# Plot Top 10 Switches
switchPlotTopSwitches(
    switchAnalyzeRlist = SwitchList_filt_Analyzed, 
    n = 10,
    pathToOutput = file.path(mydir, "ISAR","Plots"),
    filterForConsequences = FALSE, 
    splitFunctionalConsequences = TRUE
)

# Summarize switching features with consequences
extractTopSwitches(
    SwitchList_filt_Analyzed, 
    filterForConsequences = TRUE, 
    n = 10, 
    sortByQvals = TRUE
)


extractConsequenceSummary(
    SwitchList_filt_Analyzed,
    consequencesToAnalyze='all',
    plotGenes = FALSE, 
    asFractionTotal = FALSE      # enables analysis of fraction of significant features
)


# Volcano plot
pdf("Volcano_plot.pdf")
ggplot(data=SwitchList_filt_Analyzed$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
        size=0.5
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    facet_wrap( ~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('Signficant\nIsoform Switch', values = c('gray','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
dev.off()

pdf("Volcano_plot_PTC.pdf")
ggplot(data=SwitchList_filt_Analyzed$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
     geom_point(
        aes( color=PTC & isoform_switch_q_value < 0.05 & abs(dIF) > 0.1 ), # default cutoff
        size=0.5
    ) +
    geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
    facet_wrap( ~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    scale_color_manual('PTC', values = c('gray','red')) +
    labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
dev.off()

write.csv(SwitchList_filt_Analyzed$isoformFeatures, file="SwitchList_filt_Analyzed.csv")
save.image(file='ISAR_session.RData')

writeLines(capture.output(sessionInfo()), paste0(mydir, "/ISAR/ISAR_session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
# Based on: http://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file

print("File 1 completed running without errors")
