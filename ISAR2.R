#!/usr/bin/env Rscript

######
# Load libraries
######

library(optparse)
library(dplyr)
library(limma)
library(ggplot2)
library(ggrepel)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))
library(extrafont)
library(cowplot)

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 2)

# Assign arguments to variables
mydir <- arguments$args[1]
myComp <- arguments$args[2]

# Create necessary folders if they don’t exist
dir.create(file.path(mydir, "Plots"), showWarnings = FALSE)
dir.create(file.path(mydir, "Plots", "ISAR"), showWarnings = FALSE)

# Read condition file
df <- read.table(file = myComp, sep = '\t', header = TRUE)

# Initialize empty data frames and vectors
combined <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(combined) <- c("isoform_id", "gene_name", "PTC")
dIF_df <- data.frame(matrix(ncol = 1, nrow = 0))
colnames(dIF_df) <- c("isoform_id")
AllCond <- c()

# Get IsoformSwitchAnalyzeR version
ISAR_version <- packageVersion("IsoformSwitchAnalyzeR")

# Loop through conditions to process data and generate plots
for (row in 1:nrow(df)) {
  cond <- unlist(strsplit(as.character(df[row, "cond"]), " "))
  path <- df[row, "path"]
  
  # Load and filter data
  d <- read.table(file = file.path(path), sep = ',', header = TRUE)
  d_sig <- d %>% filter(isoform_switch_q_value < 0.05, !grepl("^SIRV|^ERCC", isoform_id))
  d_all <- d %>% filter(!grepl("^SIRV|^ERCC", isoform_id))
  
  # Iterate over each condition
  for (i in cond) {
    AllCond <- c(AllCond, i)
    dIF <- d_sig %>% filter(condition_2 == i)
    dIF_all <- d_all %>% filter(condition_2 == i) %>%
      mutate(isoform_switch_q_value = replace(isoform_switch_q_value, isoform_switch_q_value == 0, 1e-320))
    dIF_sig <- dIF_all %>% filter(isoform_switch_q_value < 0.05 & abs(dIF) > 0.1)
    
    # Create a directory for each condition
    dir.create(file.path(mydir, "Plots", "ISAR", i), showWarnings = FALSE)

    # Volcano plot
    myplot <- ggplot(dIF_all, aes(x = dIF, y = -log10(isoform_switch_q_value))) + 
      theme_classic() + 
      theme(
        legend.position = "top",
        strip.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_rect(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(size = 10),
        text = element_text(family = "Arial")
      ) +
      coord_cartesian(xlim = c(-1, 1), ylim = c(0, 330), clip = 'off') +
      geom_point(data = filter(dIF_all, PTC == FALSE), color = "royalblue", alpha = 0.5) +
      geom_point(data = filter(dIF_all, PTC == TRUE), color = "red", alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_hline(yintercept = 320, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = 0.1, linetype = "dashed") +
      geom_vline(xintercept = -0.1, linetype = "dashed") +    
      labs(
        title = paste0("Control_vs_", i),
        subtitle = "ISAR DTU analysis",
        x = "delta isoform fraction (dIF)",
        y = "-log10 q-value"
      ) +
      scale_color_identity(
        name = "PTC",
        breaks = c("royalblue", "red"),
        labels = c("FALSE", "TRUE"),
        guide = "legend"
      )
    
    # Density plot
    myplot2 <- ggplot(dIF_sig, aes(dIF)) + 
      theme_classic() +
      theme(
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)
      ) +
      xlim(-1, 1) +
      geom_density(data = filter(dIF_sig, PTC == FALSE), color = "royalblue", fill = "royalblue", alpha = 0.2) +
      geom_density(data = filter(dIF_sig, PTC == TRUE), color = "red", fill = "red", alpha = 0.2) +
      labs(
        x = "delta isoform fraction (dIF)",
        y = "filtered density",
        caption = paste0("Filter: |dIF| > 0.1 & q-value < 0.05 \nISAR (", ISAR_version, ")")
      ) +
      scale_color_identity(
        name = "PTC",
        breaks = c("royalblue", "red"),
        labels = c("FALSE", "TRUE"),
        guide = "legend"
      )

    # Combine and save plots
    finalplot <- plot_grid(myplot, myplot2, ncol = 1, align = "v", rel_heights = c(2, 1))
    ggsave(
      file.path(mydir, "Plots", "ISAR", i, paste0(i, "_Volcano_Density.pdf")),
      plot = finalplot, width = 10, height = 15, units = "cm", device = cairo_pdf, bg = "transparent"
    )

    # Add data to combined and dIF_df data frames for further analysis
    dIF_Barcode <- dIF %>% select(isoform_id, dIF) %>%
      rename_at(vars(dIF), ~ paste0("dIF_", i))
    dIF_df <- merge(dIF_df, dIF_Barcode, by = "isoform_id", all = TRUE)
    
    dIF_filt <- dIF %>% filter(abs(dIF) > 0.1) %>%
      rename_at(vars(dIF), ~ paste0("dIF_", i)) %>%
      rename_at(vars(isoform_switch_q_value), ~ paste0("qvalue_", i))
    combined <- merge(combined, dIF_filt, by = c("isoform_id", "gene_name", "PTC"), all = TRUE)
  }
}

# Write combined results to file
write.csv(combined, file = file.path(mydir, "Plots", "ISAR", "ISAR_combined.csv"))

