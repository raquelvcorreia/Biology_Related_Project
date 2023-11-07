# Packages to help tidy our data
library (tidyverse)
# Packages for the graphical analysis section
library (repr)
# packages used for loading data from Excel
library (readxl)



# read in FPKM-summarized expression data from an Excel file
expression_data.df <- read_excel("DataViz/Labs/Week3/rpkmDFeByg_new.xlsx", sheet=1)
str(expression_data.df)

# read in an output from DESeq2 for our experiment, with an added column of log10 of average expression value
DEGs.df <- read_excel("DataViz/Labs/Week3/DEGs.xlsx", sheet=1)
str(DEGs.df)



# merging significance values (FDRs) in DEG list into our
# expression_data data frame
# merge data sets using AGI ID as a key. These were in unlabeled columns 
# that were renamed to '0001' when importing the data from excel
expression_data_w_FDR.df <- merge(expression_data.df, DEGs.df, by.x = "...1", by.y = "...1", all.x = TRUE, all.y = TRUE)
str(expression_data_w_FDR.df)



# create a subset of genes that have an FDR rate of less 
# than 10% and 2-fold change or more
FDR10.df <- subset(expression_data_w_FDR.df, AP3_TRL_FDR < .1 & abs(AP3_TRL_logFC) > 1)
str(FDR10.df)


# create a subset of genes that have an FDR rate of less 
# than 1% and 2-fold change or more
FDR1.df <- subset(expression_data_w_FDR.df, AP3_TRL_FDR < .01 & abs(AP3_TRL_logFC) > 1)
str(FDR1.df)


# create new data frame with just expression values of FDR10 or FDR1 genes
top_FDR10.df <- FDR10.df[ , c(1,2,3,4,5)]
str(top_FDR10.df)
top_FDR1.df <- FDR1.df[ , c(1,2,3,4,5)]
str(top_FDR1.df)



# the heatmap function requires a data matrix
m <- (as.matrix(top_FDR10.df[, -1]))
top_FDR10.df$...1
rownames(m) <- top_FDR10.df$...1
heatmap(m, Colv=NA, Rowv=NA, col=rev(heat.colors(256)), mar = c(8, 6))


m <- (as.matrix(top_FDR1.df[, -1]))
rownames(m) <- top_FDR1.df$...1
heatmap(m, Colv=NA, Rowv=NA, col=rev(heat.colors(256)),  mar = c(8, 6))





# let's add a sidebar with colours denoting baseMean (~average) expression level, after log10-transformation
# We'll cluster by rows by excluding the Rowv = NA parameter
library(RColorBrewer) # we explored this library in Lab 1! 
my_group <- as.numeric(FDR1.df[FDR1.df$...1 %in% rownames(m), "Log_level"])

# You need to add 1 when selecting colour groups as our log_level values can range from 0 .. n
colSide <- brewer.pal(9, "Greys")[my_group + 1]

heatmap(m, Colv=NA, RowSideColors=colSide, col=rev(heat.colors(256)), scale = "row")

# Calculate the min, median, and max values
val_range <- c(min(m), median(m), max(m))

# Change the margins of the plot (the first is the bottom margin)
par(mar = c(6, 4.1, 4.1, 2.1), xpd=TRUE)

# Plot a corresponding legend for colour range
legend(x=-0.1, y = 0.2, legend=c("min", "med", "max"),fill=rev(heat.colors(3)), bty="n")

# Make the legend for the log_level of base mean expression
log_legend <- data.frame(log_level = unique(my_group), 
                         log_col = unique(colSide), 
                         stringsAsFactors = FALSE) %>% 
  # Sort ascending by log level
  arrange(log_level)
str(log_legend)
legend(x=-0.1, y=0.05, 
       legend=log_legend$log_level, fill=log_legend$log_col, 
       title = "Log10 level", bty="n", ncol = 2)


# let's add a sidebar with colours denoting  baseMean (~average) expression level, after log10-transformation
# We'll cluster by rows by excluding the Rowv = NA parameter
library(RColorBrewer)
# my_group <- as.numeric(as.factor(substr(rownames(m), 1 , 3)))
my_group <- as.numeric(FDR1.df[FDR1.df$...1 %in% rownames(m), "Log_level"])

# You need to add 1 when selecting colour groups as our log_level values can range from 0 .. n
colSide <- brewer.pal(9, "Greys")[my_group + 1]
heatmap(m, Colv=NA, RowSideColors=colSide, col=rev(heat.colors(256)), scale = "none")

# Calculate the min, median, and max values
val_range <- c(min(m), median(m), max(m))

# Change the margins of the plot (the first is the bottom margin)
par(mar = c(6, 4.1, 4.1, 2.1), xpd=TRUE)

# Plot a corresponding legend for colour range
# legend(x="bottomright", inset=c(-0.1,-0.15), legend=c("min", "med", "max"),fill=rev(heat.colors(3)), bty="n")
legend(x=-0.1, y = 0.2, legend=c("min", "med", "max"),fill=rev(heat.colors(3)), bty="n")

# Make the legend for the log_level of base mean expression
log_legend <- data.frame(log_level = unique(my_group), 
                         log_col = unique(colSide), 
                         stringsAsFactors = FALSE) %>% 
  # Sort ascending by log level
  arrange(log_level)

legend(x=-0.1, y=0.05, 
       legend=log_legend$log_level, fill=log_legend$log_col, 
       title = "Log10 level", bty="n", ncol = 2)

