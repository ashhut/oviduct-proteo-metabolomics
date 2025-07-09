## Proteomics Analysis

library(DEP)
library(dplyr)
library(data.table)
setwd("proteomics/")
inputDir = "analyses/maxquant/combined/txt/"
outputDir = "analyses/dep/"
plotDir = "output/plots/"

data <- read.table(paste0(inputDir,"proteinGroups.txt"), header=TRUE, sep="\t") ## import data table
data <- filter(data, Reverse != "+", Potential.contaminant != "+")  ## # We filter for contaminant proteins and
## decoy database hits, which are indicated by "+" in the columns "Potential.contaminants" and "Reverse", respectively.

data$Protein.IDs_1 <- gsub('scaffold.*', '', data$Protein.IDs)
data$Protein.IDs_1 <- gsub('tr', '', data$Protein.IDs_1)
data$Protein.IDs_1[3] <- "|CCM|CCM"

data$Protein.IDs_2 <- unlist(lapply(strsplit(as.character(data$Protein.IDs_1), "[|]"), '[[', 2))
data$Protein.IDs_3 <- unlist(lapply(strsplit(as.character(data$Protein.IDs_1), "[|]"), '[[', 3))

data$Protein.IDs_3 %>% duplicated() %>% any()
data %>% group_by(Protein.IDs_3) %>% summarize(frequency = n()) %>%
  arrange(desc(frequency)) %>% filter(frequency > 1)
data_unique <- make_unique(data, "Protein.IDs_3", "Protein.IDs_3", delim = ";")

LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
experimental_design <- fread(paste0(outputDir,"experimental_design"))
data_se <- make_se(data_unique, LFQ_columns, experimental_design)

# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se_parsed <- make_se_parse(data_unique, LFQ_columns)

# Let's have a look at the SummarizedExperiment object
data_se
plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
## This is VERY stringent, also given the fact that your non-oestrus group isn't super tight
## It might be worth making it less stringent for filtering proteins
data_filt <- filter_missval(data_se, thr = 0)

# ASH Could try the Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt2 <- filter_missval(data_se, thr = 1)

# Filter for proteins that are quantified in at least 2/3 of the samples.
frac_filtered <- filter_proteins(data_se, "fraction", min = 0.66)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
# pdf("proteins_overlap_between_samples.pdf")
plot_coverage(data_filt)
# dev.off()

# Normalize the data
# The data is background corrected and normalized by variance stabilizing transformation (vsn).
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
pdf("normalised_data_plot.pdf")
plot_normalization(data_filt, data_norm)
dev.off()

# Plot a heatmap of proteins with missing values
pdf("missing_values.pdf")
plot_missval(data_filt)
dev.off()

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
pdf("intensity_dist_of_missing_values.pdf")
plot_detect(data_filt)
dev.off()

# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(data_norm, fun = "")

## ASH: You could try the different imputation methods
## Not sure what is best to use, I used the first imputation method below for the following analyses
# 1. Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# 2. Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# 3. Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)


# Plot intensity distributions before and after imputation
pdf("intensity_dist_before_and_after_imputation.pdf")
plot_imputation(data_norm, data_imp)
dev.off()

## Differential Analysis
# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "NOE")

# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")


# Test manually defined comparisons
## ASH: This is perhaps where you could try and include the comparisons with the blood and tissue
## To see where they fit with the samples
## I'm not sure how to "control" or "minus" these from the analysis but maybe you can see what
## overlaps with the fluid samples for both control and oestrus and then look at differential proteins between
## the fluid and tissues. Because I guess you would expect the proteins that are in both to be higher in the
## tissue?
## And then can look at what is not in both groups as well?

data_diff_manual <- test_diff(data_imp, type = "manual",
                              test = c("Blood", "OE")) ## put in the actual IDs you want to compare


# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

# Plot the first and second principal components
pdf('pca.pdf')
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)
dev.off()

## new code for PCA
plot_pca(dep, x = 1, y = 2, n = 500, label=TRUE, indicate ="condition")
## to change label make it TRUE or FALSE

# Plot the Pearson correlation matrix
pdf("plot_pearson_correlaton.pdf")
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
dev.off()

# Plot a heatmap of all significant proteins with the data centered per protein
pdf("heatmap_all_sig_proteins.pdf")
plot_heatmap(dep, type = "centered", kmeans = TRUE, k=6,
             col_limit = 2, show_row_names = FALSE, indicate = c("condition"))
dev.off()

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
pdf("heatmap_all_sig_contrast_proteins.pdf")
plot_heatmap(dep, type = "contrast", kmeans = TRUE,
             k = 6, col_limit = 10, show_row_names = FALSE, )
dev.off()

# Plot a volcano plot for the contrast "OE_vs_NOE""
pdf("volcano_plot.pdf")
plot_volcano(dep, contrast = "OE_vs_NOE", label_size = 2, add_names = TRUE)
dev.off()


# Plot a barplot for whatever proteins of interest you want to look at
## If you have any that you want to specifically see if they're differentially expressed.
pdf("test.pdf")
plot_single(dep, proteins = c("CCM", "AGPS"))
dev.off()

pdf("CCM_log2centred_intensity.pdf")
plot_single(dep, proteins = "CCM", type = "centered")
dev.off()

# Plot a frequency plot of significant proteins for the different conditions
pdf("freq_sig_proteins.pdf")
plot_cond(dep)
dev.off()

# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()
colnames(data_results)
write.table(data_results, file="data_differential_results.txt", quote=FALSE, sep="\t")

# Generate a wide data.frame
df_wide <- get_df_wide(dep)
write.table(df_wide, file="data_wide_results.txt", quote=FALSE, sep="\t")

# Generate a long data.frame
df_long <- get_df_long(dep)
write.table(df_long, file="data_long_results.txt", quote=FALSE, sep="\t")

# Save analyzed data
save(data_se, data_norm, data_imp, data_diff, dep, file = "data.RData")
# These data can be loaded in future R sessions using this command
load("data.RData")


### LFQ-based DEP analysis ###
# The wrapper function performs the full analysis
data$Gene.names <- unlist(lapply(strsplit(as.character(data$Protein.IDs_1), "[|]"), '[[', 3))

data_results <- LFQ(data, experimental_design, fun = "MinProb",
                    type = "control", control = "NOE", alpha = 0.05, lfc = 1)

# Make a markdown report and save the results
report(data_results)

# See all objects saved within the results object
names(data_results)

# Extract the results table
results_table <- data_results$results

# Number of significant proteins
results_table %>% filter(significant) %>% nrow()

# Extract the sign object
full_data <- data_results$dep

# Use the full data to generate a heatmap
pdf("heatmap_LFQ_based_diff.pdf")
plot_heatmap(full_data, type = "contrast", kmeans = TRUE,
             k = 6, col_limit = 2, show_row_names = FALSE)
dev.off()



