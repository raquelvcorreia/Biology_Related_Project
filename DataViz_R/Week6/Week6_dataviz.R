library(UpSetR)


# load in data
#mutations <- read.csv( system.file("extdata", "DataViz/Labs/Week6/mutations.csv", package = "UpSetR"), header=TRUE, sep = ",")
mutations <- read.csv("DataViz/Labs/Week6/mutations.csv", header=TRUE, sep = ",")
str(mutations)

upset(mutations, sets = c("PTEN", "TP53", "EGFR", "PIK3R1", "RB1"), sets.bar.color = "#56B4E9", order.by = "freq", empty.intersections = "on")


# Packages to help tidy our data
library(tidyverse)
library(readxl)
# Packages for the graphical analysis section
library(RColorBrewer)
# Data projection packages
library(umap)

# Read in our RNAseq data (.gz file read as text, no need to unzip)
tissue_data.df <- read.table(file = "DataViz/Labs/Week6/GSE150316_DeseqNormCounts_final.txt.gz",
                             header = TRUE,
                             row.names = 1)
# Take a quick look at it
head(tissue_data.df)
dim(tissue_data.df)
#Read in some additional patient data
patient_data.df <- read_excel("DataViz/Labs/Week6/2020.07.30.20165241-1_supp_table3.xlsx", sheet=1)
# Take a quick look at it
head(patient_data.df)
dim(patient_data.df)

# Reformat our patient data, store in patient_viral_load.df
patient_viral_load.df <-
  patient_data.df %>%
  rename_with(str_replace_all, pattern="\\r\\n|\\s", replacement = "_") %>%
  select(1:3)


patient_viral_load.df



# Genes where transcript levels are very low, or which don't show a lot of 
#variance, are not that informative for our UMAP analysis, and will just slow it down. 
# Trim the tissue data down by removing the data on these genes

tissue_data_filtered.df <-
  tissue_data.df %>%
  # Convert the row names to a column
  rownames_to_column(var="gene") %>%
  # Set up the table to perform row-wise operations
  rowwise() %>%
  # Calculate the mean expression of each gene across all 
  # tissue samples
  mutate(mean = mean(c_across(where(is.numeric)))) %>%
  # Filter for samples with low expression
  filter(mean > 0.5) %>%
  # Calculate overall variance in case we need to make our dataset smaller
  mutate(variance = var(c_across(where(is.numeric)))) %>%
  # Arrange samples by descending variance
  arrange(desc(variance)) %>%
  # Remove the grouping specification
  ungroup()

head(tissue_data_filtered.df)
dim(tissue_data_filtered.df)

# We need to transpose the data.
# We can do it with dplyr to keep it as a data frame and to add some info
tissue_RNAseq.df <-
  tissue_data_filtered.df %>%
  select(1:89) %>% # trim down the columns
  pivot_longer(cols=c(2:89), names_to = "sample", values_to = "norm_counts") %>%
  pivot_wider(names_from = gene, values_from = norm_counts)

tissue_RNAseq.df

# We want to add some additional sample information before
# assessing the data
tissue_RNAseq.df <-
  tissue_RNAseq.df %>%
  # Grab just the sample names
  select(sample) %>%
  # Grab information from it like case number, tissue, and
  # tissue number
  str_match_all(., pattern=c("case([\\w]+)\\.([a-z]+)([\\d|\\.NYC]*)|(NegControl\\d)")) %>%
  # Bind it all together
  do.call(rbind, .) %>%
  # Convert the results to a data frame and DO NOT convert strings as factors
  as.data.frame(stringsAsFactors = FALSE) %>%
  # Rename the columns based on the capture groups
  rename(., sample = V1, case_num = V2, tissue = V3, tissue_num = V4, neg_num = V5) %>%
  # Coalesce some of the info due to negative control samples
  # and clean up a column
  mutate(case_num = coalesce(case_num, neg_num),
         tissue_num = str_replace_all(.$tissue_num, pattern = "\\.", replace = ""), 
         tissue = replace_na(tissue, "None")  # Replace the negative control tissues to "None"
  ) %>%
  # Drop the neg_num column
  select(1:4) %>%
  # Join this result to the RNA-seq info
  full_join(., y=tissue_RNAseq.df, by=c("sample" = "sample")) %>%
  # Join that result to grab viral load information
  right_join(patient_viral_load.df, y=., by=c("Case_No" = "case_num")) %>%
  # Fix some column names
  rename(case_num = Case_No,
         viral_load = `Viral_high_vs._viral_low*`,
         viral_load_percent = `Viral_load`)
head(tissue_RNAseq.df)



# How many tissue types do we have? and how many samples from each tissue type
table(tissue_RNAseq.df$tissue)



# Generate a matrix version of our data but drop the sample
# information!
tissue_RNAseq.mx <- as.matrix(tissue_RNAseq.df[,c(-1:-6)])



# Set our seed
set.seed(1981)
# Generate our projection
tissue_umap <- umap(tissue_RNAseq.mx)


# load in data
str(tissue_umap)

tissue_umap$layout[,1]
tissue_umap$layout[,2]

# Re-map our projection points with our tissue data
tissue_umap.df <- data.frame(x.coord = tissue_umap$layout[,1],
                             y.coord = tissue_umap$layout[,2])
tissue_umap.df <- cbind(tissue_RNAseq.df[,1:6], tissue_umap.df)
tissue_umap.df <-
  tissue_umap.df %>%
  mutate(viral_load = replace_na(viral_load, replace = "DNW"))



# Adjust our plot window size according to the expected
# output
options(repr.plot.width=20, repr.plot.height=20)
combo.colours = c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"), brewer.pal(9, "Set1"))
# try combining some other palettes, see
# https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/#rcolorbrewer-palettes
# combo.colours = c(brewer.pal(8, "Dark2"),brewer.pal(8, "Paired"))
# 1. Data
ggplot(data = tissue_umap.df) +
  # 2. Aesthetics
  aes(x = x.coord, y = y.coord, colour = tissue, shape = viral_load, ) +
  # Themes
  theme_bw() +
  theme(text = element_text(size=20)) +
  # 3. Scaling
  scale_colour_manual(values = combo.colours) +
  # 4. Geoms
  geom_text(aes(label = case_num), size = 10)
# make the data point markers bold:
# geom_text(aes(label = case_num), size = 10, fontface="bold")

# Adjust our plot window size according to the expected
# output
options(repr.plot.width=20, repr.plot.height=20)
# combo.colours = c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"), brewer.pal(9, "Set1"))
# try combining some other palettes, see
# https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/#rcolorbrewer-palettes
combo.colours = c(brewer.pal(8, "Dark2"),brewer.pal(8, "Paired"))
# 1. Data
ggplot(data = tissue_umap.df) +
  # 2. Aesthetics
  aes(x = x.coord, y = y.coord, colour = tissue, shape = viral_load, ) +
  # Themes
  theme_bw() +
  theme(text = element_text(size=20)) +
  # 3. Scaling
  scale_colour_manual(values = combo.colours) +
  # 4. Geoms
  # geom_text(aes(label = case_num), size = 10)
  # make the data point markers bold:
  geom_text(aes(label = case_num), size = 10, fontface="bold")

