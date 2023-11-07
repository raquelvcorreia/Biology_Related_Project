# Packages to help tidy our data
library (tidyverse)
# Packages for the graphical analysis section
library (repr)
library (viridis)
# packages used for working with/formatting dates in R
library (lubridate)
library (zoo)

# Initialize a plot with our data
plot_1.df <- read_tsv('DataViz/Labs/Week1/Lab1_scatterplot_set_1.txt')

# standard deviation of set 1 
sd(plot_1.df$x_values)  

# Take a quick look at the structure of the data
str(plot_1.df)

# Instantiate ggplot object
plot_1.plot <- ggplot(plot_1.df)

# Update the aesthetics with axis and colour information, then add a scatter graph!
plot_1.plot <- plot_1.plot +
  aes(x = x_values, y = y_values, colour = dataset) +
  geom_point() +
  theme(text = element_text(size = 20)) + # set text size
  guides(colour = guide_legend(title="Set 1")) + # Legend title
  xlab("x value") + # Set the x-axis label
  ylab("y value") + # Set the y-axis label
  xlim(0,100) +     # specify the minimum and maximum x axis values
  ylim(0,100)       # ditto for the y axis

#display our plot
plot_1.plot


# Save the plot we've generated to the root directory of the files.
ggsave(plot = plot_1.plot, filename = "Set_1.png", scale=2, device = "png", units = c("cm"))



# Initialize a plot with our data
plot_2.df <- read_tsv('DataViz/Labs/Week1/Lab1_scatterplot_set_2.txt')


# Instantiate ggplot object
plot_2.plot <- ggplot(plot_2.df)

# Update the aesthetics with axis and colour information, then add a scatter graph!
plot_2.plot <- plot_2.plot +
  aes(x = x_values, y = y_values, colour = dataset) +
  geom_point() +
  theme_bw() + #check out options https://ggplot2-book.org/polishing.html
  theme(text = element_text(size = 20)) + # set text size
  guides(colour = guide_legend(title="Set 2")) + # Legend title
  xlab("x value") + # Set the x-axis label
  ylab("y value") + # Set the y-axis label
  xlim(0,100) +     # specify the minimum and maximum x axis values
  ylim(0,100)       # ditto for the y axis

#display our plot
plot_2.plot

# Initialize a plot with our data
plot_3.df <- read_tsv('DataViz/Labs/Week1/Lab1_scatterplot_set_3.txt')

# Instantiate ggplot object
plot_3.plot <- ggplot(plot_3.df)


# Update the aesthetics with axis and colour information, then add a scatter graph!
plot_3.plot <- plot_3.plot +
  aes(x = x_values, y = y_values, colour = dataset) +
  geom_point() +
  theme_bw() + #check out options https://ggplot2-book.org/polishing.html
  theme(text = element_text(size = 20)) + # set text size
  guides(colour = guide_legend(title="Set 3")) + # Legend title
  xlab("x value") + # Set the x-axis label
  ylab("y value") + # Set the y-axis label
  xlim(0,100) +     # specify the minimum and maximum x axis values
  ylim(0,100)       # ditto for the y axis

# use a different colour palette, see options https://ggplot2-book.org/scale-colour.html
plot_3.plot <- plot_3.plot + 
  scale_colour_brewer(palette = "Pastel2")

#display our plot
plot_3.plot



# Initialize a plot with our data
plot_4.df <- read_tsv('DataViz/Labs/Week1/Lab1_scatterplot_set_4.txt')


# Instantiate ggplot object
plot_4.plot <- ggplot(plot_4.df)


# Update the aesthetics with axis and colour information, then add a scatter graph!
plot_4.plot <- plot_4.plot +
  aes(x = x_values, y = y_values, colour = dataset) +
  geom_point() +
  theme_bw() + #check out options https://ggplot2-book.org/polishing.html
  theme(text = element_text(size = 20)) + # set text size
  guides(colour = guide_legend(title="Set 4")) + # Legend title
  xlab("x value") + # Set the x-axis label
  ylab("y value") + # Set the y-axis label
  xlim(0,100) +     # specify the minimum and maximum x axis values
  ylim(0,100)       # ditto for the y axis

# use a different colour palette, see options https://ggplot2-book.org/scale-colour.html
plot_4.plot <- plot_4.plot + 
  scale_colour_manual(values= c("blue1", "blue2", "blue3", "blue4"))
# we only use the first colour ("blue1") here, as there is only one dataset

#display our plot
plot_4.plot

