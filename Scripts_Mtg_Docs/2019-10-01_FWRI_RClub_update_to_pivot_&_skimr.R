#### October 1, 2019 ####
#### 1. Announcements ####

# We're scheduled to meet the 2nd Tuesday of each month, 2 pm in the FWRI 4th Floor Conference Room
  # November 12, 2019 at 2 pm
  # December 10, 2019 at 2 pm

# Upcoming RStudio Webinars
  # https://resources.rstudio.com/upcoming-webinars/webinar-registration
  # Check out "Interactivity in Production", October 2, 2019 1 - 2 pm
  # Register on the website if you'd like to attend

# I'll try to keep the SharePoint page as updated as possible, including scripts from each R Club meeting
  # https://fwcc.sharepoint.com/sites/teams/GISUF/SitePages/R%20Users%27%20Forum.aspx

# Send topics of or post them on the SharePoint page

#### 2. A Couple R Studio Tips ####

# In your files pane (by your Plots, Packages, Help, etc.), you can:
  # View all files in your current working directory
  # Navigate folders, create folders
  # Set & Go To working directories with the "More" tab

# Filter/Explore your data when viewing
  # View your dataset by selecting from the environment pane
  # Click the "Filter" button
  # Filter by categorical or numerical variables
  # Get histograms of numerical variables
  # Test it out below:
chicks <- data.frame(ChickWeight) # Built-in R data set, Weight versus age of chicks on different diets

#### 3. Tidy Tuesday ####
# https://thomasmock.netlify.com/post/tidytuesday-a-weekly-social-data-project-in-r/ 
# A weekly social data project in R
# Could be good practice if you're looking for practicing in the tidyverse

#### 4. Update to pivot() in tidyr package ####
  # The vignette: https://tidyr.tidyverse.org/dev/articles/pivot.html 
  # 2 New functions provide more flexibility for gather() and spread()
  # Inspired by melt() and cast() from data.table package and attributes of the cdata package
    # pivot_wider()
    # pivor_longer()
library(tidyr)
library(dplyr)
library(readr)
# Load an R dataset
relig_income

test <- relig_income %>% 
  # reshape ever column apart from religion
  pivot_longer(-religion, names_to = "income", values_to = "count")
  #names_to gives name of variable that will be created and stored in column names
  #values_to gives name of variable created from the data stored in the cell value
str(test$income)

#Try another R dataset
billboard

# the data encoded in the column names is really a number, not a string
billboard %>% 
  pivot_longer(
    cols = starts_with("wk"), 
    names_to = "week", 
    values_to = "rank",
    values_drop_na = TRUE
  )

#determine how long each song stayed in the charts
#need to convert week to integer
billboard %>% 
  pivot_longer(
    cols = starts_with("wk"), 
    names_to = "week", 
    names_prefix = "wk",
    #provide name-prototype pair that defines the type, class, and attributes of a vector, in this case integer()
    names_ptypes = list(week = integer()),
    values_to = "rank",
    values_drop_na = TRUE,
  )

# What if you have multiple values in the column names?
# Yet another R dataset
who
who %>% pivot_longer(
  cols = new_sp_m014:newrel_f65,
  #break the column name values apart with names to and
  names_to = c("diagnosis", "gender", "age"), 
  #then provide pattern for splitting; it's similar to extract
  #give it a regular expression containing groups (defined by () ) and it puts each group in a column
  names_pattern = "new_?(.*)_(.)(.*)",
  values_to = "count"
)

#### 5. skimr package ####
library(skimr)
# "Compact and Flexible Summaries of Data" . . . in the console
# alternative to default summary functions within R

# Summarize a whole dataset
skim(chicks)

# Select specific columns to summarize
skim(chicks, weight, Time)

# Also works with piping and grouped data!
chicks %>%
  dplyr::group_by(Diet) %>%
  skim()

# You can also specify your own stats using a list combined with the skim_with() function
# Supports any names class in your data

# Here let's try the interquartile range and the median absolute deviation
funs <- list(
  iqr = IQR,
  mad = mad
)

skim_with(numeric = funs, append = FALSE)
skim(chicks, weight)