# Packages to help tidy our data
library (tidyverse)
# Packages for the graphical analysis section
library (repr)
# packages used for working with/formatting dates in R
library (lubridate)
library (zoo)


covid_demographics.df <- 
  read_csv("DataViz/Labs/Week2/Ontario_age_and_sex_COVID-19_data_merged_filtered.csv") 
str(covid_demographics.df)
tail(covid_demographics.df)
head(covid_demographics.df, n=50)


# Format our dataset to focus on total_cases and total_deaths 
#(chain operator (%>%) allows the output of one operation to feed into the next.)
covid_demographics_total.df <-
  covid_demographics.df %>%
  # Pare down the dataset to just total_cases, total_deaths, and
  # total_hospitalizations
  select(from_date, to_date, geographic_area, age_group,
         total_cases, total_deaths, total_hospitalizations_count) %>%
  # Convert the age_group into a "factor"
  mutate(age_group = factor(age_group)) %>%
  rename(public_health_unit = geographic_area) %>%
  # Group the data so you can "summarize" in the next steps
  group_by(public_health_unit) %>%
  # Generate percent cases for each age group within a PHU
  mutate(percent_cases = total_cases/sum(total_cases),
         # Generate percent deaths for each age group within a PHU
         percent_deaths = total_deaths/sum(total_deaths),
         # Generate % hospitalizations for each age group within a PHU
         percent_hospitalizations =
           total_hospitalizations_count/sum(total_hospitalizations_count))
# Take a look at the different age demographics
levels(covid_demographics_total.df$age_group)

## Plotting the distribution of cases by age group
# Adjust our plot window size according to the expected output
options(repr.plot.width=20, repr.plot.height=10)
covid_demographics_total.df %>%
  # Filter for some of our age groups
  filter(age_group %in% c("20 to 39","40 to 59", "60 to 79", "80+")) %>%
  # 1. Data
  ggplot(.) +
  # 2. Aesthetics
  aes(x = percent_cases, fill = age_group) +
  theme(text = element_text(size = 20)) + # set text size
  # 4. Geoms
  geom_density(alpha = 0.5) + # generate kernel density estimate
  geom_rug(aes(colour = age_group)) 
# confirm our data values with a geom_rug

##Facet Plots
# Instantiate ggplot object
# Adjust our plot window size according to the expected output
options(repr.plot.width=20, repr.plot.height=20)
covid_demographics_total.df %>%
  # Select for just the important columns
  select(public_health_unit, age_group, percent_cases, percent_deaths) %>%
  # 1. Data
  ggplot(.) +
  # 2. Aesthetics
  aes(x = percent_cases, fill = age_group) +
  theme (text = element_text(size = 20)) + # set text size
  # 4. Geoms
  geom_density(alpha = 0.5) +
  # 6. Facet
  facet_wrap(~age_group, ncol = 1, scales = "free_y" ) 
# try with scales = "free_y" as a parameter


## Box Plots
# Adjust our plot window size according to the expected output
options(repr.plot.width=20, repr.plot.height=10)
covid_demographics_total.df %>%
  # Select for just the important columns
  select(public_health_unit, age_group, percent_cases, percent_deaths) %>%
  # Pivot the modified table to capture the "stat_group" of 
  # percent_cases vs percent_deaths
  pivot_longer(cols=c(3:4), names_to = "stat_group", values_to = "percent_PHU_total") %>%
  # Plot the data as a grouped boxplot
  # 1. Data
  ggplot(.) +
  # 2. Aesthetics
  aes(x=age_group, y = percent_PHU_total, fill = stat_group) +
  theme(text = element_text(size = 20)) + # set text size
  theme(legend.position = "right") +
  scale_y_continuous(limits = c(0, 0.6)) + # Set the y-axis limit
  # 4. Geoms
  geom_boxplot(outlier.shape=NA, notch = FALSE) + 
  # Add the boxplot geom
  geom_point(position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.75, seed = NA), alpha = 0.25) # Add data points for comparison



## Violin plots. It is recommended to use them when there is at least 30 datapoints
# Adjust our plot window size according to the expected output
options(repr.plot.width=20, repr.plot.height=10)
# Generate a basic box plot with outliers present
covid_demographics_total.df %>%
  # 1. Data
  ggplot(.) +
  # 2. Aesthetics
  aes(x=age_group, y = percent_cases) +
  theme(text = element_text(size = 20)) + # set text size
  theme(legend.position = "none") +
  # 4. Geoms
  geom_violin() + # Add the boxplot geom
  geom_jitter(aes(colour = public_health_unit), width = 0.1, height = 0)

