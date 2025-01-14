# Load necessary libraries
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)

#Read files FRK1 KO and KD files seperately
a <- read.csv("kd - kd.csv", row.names = 1)
b<- read.csv("ko - ko.csv", row.names = 1)

# Remove rows where PTC is FALSE or NA
a <- a %>% filter(PTC != FALSE)
b <- b %>% filter(PTC != FALSE)
#remove PTC column
a <- a[, !colnames(a) %in% c("PTC")]
b <- b[, !colnames(b) %in% c("PTC")]

#order according to q value
a <- a[order(a$isoform_switch_q_value), ]
b <- b[order(b$isoform_switch_q_value), ]

# Filter rows where isoform_switch_q_value <= 0.05
a <- a %>% filter(isoform_switch_q_value <= 0.05)
b <- b %>% filter(isoform_switch_q_value <= 0.05)

# Renaming the column
colnames(a)[colnames(a) == "isoform_switch_q_value"] <- "isoform_switch_q_value_kd"
colnames(a)[colnames(a) == "dIF"] <- "dIF_kd"
colnames(b)[colnames(b) == "isoform_switch_q_value"] <- "isoform_switch_q_value_ko"
colnames(b)[colnames(b) == "dIF"] <- "dIF_ko"

# Preserve rownames as a column for merging
a$RowNames <- rownames(a)
b$RowNames <- rownames(b)
merged_df <- merge(a, b, by = "RowNames", all = TRUE)

#removing q value column
merged_df <- merged_df[, !colnames(merged_df) %in% c("isoform_switch_q_value_kd")]
merged_df <- merged_df[, !colnames(merged_df) %in% c("isoform_switch_q_value_ko")]

# Restore row names if needed
rownames(merged_df) <- merged_df$RowNames
merged_df$RowNames <- NULL
write.csv(merged_df, "merged.csv")


#read up and down-regulated files
up <- read.csv("merged_upset plot - up.csv", row.names = 1)
down <- read.csv("merged_upset plot - down.csv", row.names = 1)

#replace na with 0
down$dIF_kd[is.na(down$dIF_kd)] <- "0"
down$dIF_ko[is.na(down$dIF_ko)] <- "0"
up$dIF_kd[is.na(up$dIF_kd)] <- "0"
up$dIF_ko[is.na(up$dIF_ko)] <- "0"
# Replace all values except 0 with 1
up$dIF_kd[up$dIF_kd != 0] <- 1
up$dIF_ko[up$dIF_ko != 0] <- 1
down$dIF_kd[down$dIF_kd != 0] <- 1
down$dIF_ko[down$dIF_ko != 0] <- 1
#Renaming the columns
colnames(up)[colnames(up) == "dIF_kd"] <- "FRG1_KD"
colnames(down)[colnames(down) == "dIF_kd"] <- "FRG1_KD"
colnames(up)[colnames(up) == "dIF_ko"] <- "FRG1_KO"
colnames(down)[colnames(down) == "dIF_ko"] <- "FRG1_KO"


# set as variable
data <- up
data <- down
# Convert to a binary matrix
data_matrix <- as.matrix(data)
data_matrix <- data_matrix > 0  # Ensure logical matrix
# Generate a combination matrix (required for UpSet plot)
comb_mat <- make_comb_mat(data_matrix)

# Visualize the UpSet plot with numbers
ht3 <- UpSet(
  comb_mat,
  set_order = colnames(data),
  comb_col = c("red", "red"), # Customize colors for combinations
  column_title = "PTC-containing isoforms",
  row_title = "Set Size",
  top_annotation = upset_top_annotation(
    comb_mat, 
    gp = gpar(col = "black"),
    add_numbers = TRUE,          # Add numbers for intersection sizes
    numbers_gp = gpar(fontsize = 10, col = "black") # Customize appearance
  ),
  right_annotation = upset_right_annotation(
    comb_mat, 
    gp = gpar(col = "black"),
    add_numbers = TRUE,          # Add numbers for set sizes
    numbers_gp = gpar(fontsize = 10, col = "black") # Customize appearance
  )
)
# Draw the UpSet plot
draw(ht3)
