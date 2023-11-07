library(tidyverse)
library(repr)
library(viridis)


# Create a vector with your file names
fileNames <- c("K60dvsw", "K40dvsw", "K20dvsw", "Krwvsw", "L60dvsw", "L40dvsw", "L20dvsw", "Lrwvsw")


# While we import the files, we'll add a column to denote the dataset that it comes from. 
# This will help us put the data into a long format
data.list <- lapply(fileNames, FUN = function(x) mutate(read_tsv(paste0("DataViz/Labs/Week4/Lab_4_R_files/", x, "_up_GO.txt")), dataSet = x))
names(data.list) <- fileNames

# What is the structure of the data set
str(data.list, give.attr = FALSE)
# Note that L60dvsw_up_GO.txt has no data in it.

# How many rows are there across all data sets?
sum(unlist(lapply(data.list, FUN = function(x) dim(x)[1])))


# Now that all the data is in a list, we can combine it into a single data frame
allData.df <- 
  do.call(rbind, data.list) %>% 
  mutate(enrich = queryitem/bgitem)

# What does the dataframe look like?
str(allData.df)




# Now that all the data is in a list, we can combine it into a single data frame
allData.df <- 
  do.call(rbind, data.list) %>% 
  mutate(enrich = queryitem/bgitem)

# What does the dataframe look like?
str(allData.df)


# Find the significant GO terms you're interested in
allSigTerms <-
  allData.df %>% 
  # Filter the data
  filter(FDR <= 0.01, term_type == "P") %>% 
  # What are the terms left after filtering?
  pull(Term)

allSigTerms
# How many unique terms are there?
str(unique(allSigTerms))




# Import the specific term list

goTermKeep <- read_file("DataViz/Labs/Week4/Lab_4_R_files/keep_formakinggofigures.txt") %>% str_split(pattern = "\r\n") %>% unlist()

# What does it look like
goTermKeep

# What are the attributes
str(goTermKeep)



# What happens when we filter by these terms?
sigData.df <-
  allData.df %>% 
  # This mutation step will return L60dvsw back as an x-axis value
  mutate(dataSet = factor(dataSet, levels = fileNames)) %>% 
  filter(Term %in% goTermKeep) %>% 
  mutate(Term = factor(Term, levels = goTermKeep)) %>% 
  # Filter the data again based on the lab notes
  filter(FDR <= 0.01, term_type == "P")

str(sigData.df)


options(repr.plot.width = 16, repr.plot.height = 10)
# Set a maximum dot size
max_dotsize <-6
enrichMax <- max(sigData.df$enrich)
min_FDR = min(sigData.df$FDR)

# The hardest part of printing this is you want to maintain (as much of) the original goTermKeep values
# On the other hand you want to drop any data where FDR > 0.01 or the term_type != "P"

# If you filter the dataset this way, you'll lose y-axis values, even if converted to a factor with levels!
# To retain those missing values/levels from the factor, set scale_y_discrete(drop = FALSE)

ggplot(sigData.df) +
  # 2. Aesthetics and Theme
  aes(x = dataSet, y = Term, colour = FDR, size = enrich) + #Set the y and x axes
  
  labs(x = NULL, y = NULL) +
  theme(axis.title.x = element_text(size = rel(0.7))) +
  
  # 3. Scaling
  scale_size_continuous(range = c(1, max_dotsize), 
                        limits = c(0, enrichMax), 
                        breaks = seq(0, enrichMax, 0.1),
                        guide = guide_legend(order=2)
  ) +
  scale_y_discrete(drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  #    scale_colour_gradient(limits = c(min_FDR, 0.01), 
  #                           low="gray36", 
  #                           high="gray87",
  #                           guide = guide_colorbar(order=1)
  #                          ) +
  scale_colour_viridis(guide = guide_colorbar(order=1)) +
  
  # 4. Geoms
  geom_point(aes(size = enrich, color = FDR))


# In the above graph the spacing between the samples is a bit wide (depending on your screen)
# When we save this plot as a pdf, try playing around with the width and height parameters
# The ggsave function will automatically adjust column spacing, fonts etc.
# The following width and height settings result in the plot you see in the lab manual
ggsave("bubble_plot.pdf", width = 25,  height = 30,  units = "cm")