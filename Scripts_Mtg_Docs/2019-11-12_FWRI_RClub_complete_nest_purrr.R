#### FWRI R Club: November 12, 2019 ####

#### 1. Announcements ####
 # Tampa Bay R Users Group Meetup = Nov. 19, 7 - 9 pm at Southern Brewing and Winemaking
 #    Topic = creating R packages

 # First meeting of 2020! == January 14, 2020

#### 2. tidyr::complete() ####
 # The default behavior in R dataframes is that if no data exists
 #  for a particular observation, then the row for that observation
 #  does not appear in the dataframe. This can cause problems when
 #  you need to use this dataframe as an input for something which
 #  expects to see values for all possible observations.
 # For me, I think abundance data for a species and we want both positive
 #  and zero catches to be retained for each sample . . . 
 #  Or if you are making future projections and the starting point has missing
 #  rows
 # Check out this blog post: http://www.imachordata.com/you-complete-me/
 # The complete() function within tidyr allows you to fill in the gaps
 #  for all observations that had no data. It allows you to define the
 #  observations that you want to complete and then declare what value to
 #  use to plug the gaps.
 # For example, if you were taking counts of different species of kelp for
 #  different years, and you had some combinations for which there were
 #  no individuals of a species, you could use the following to deal with it:

library(tidyverse)

 # Example from the blog:
kelpdf <- data.frame(
  Year = c(1999, 2000, 2004, 1999, 2004),
  Taxon = c("Saccharina", "Saccharina", "Saccharina", "Agarum", "Agarum"),
  Abundance = c(4,5,2,1,8)
)
kelpdf
# You will notice Agarum wasn't recorded in 2000, even though other
# species were recorded in that year

# So let's fix that with tidyr::complete()
library(tidyr)

kelpdf %>% complete(Year, Taxon)
# Complete looked at all possible combinations of Year and Taxon
# when there was a missing combination, it inserted NA

# We can specify what we want those missing values to be if not NA
# In this case, we want to change it to 0
kelpdf %>% complete(Year, Taxon, fill = list(Abundance = 0))

# But, what if we also know that other years were sampled and we want
# the zero data from those years? (The original data was a subset)
# Combine with full_seq()

kelpdf %>% complete(Year = full_seq(Year, period = 1),
                    Taxon,
                    fill = list(Abundance = 0))
# Note we had to define (or overwrite, really) the Year column

#### 3. tidyr::nest() ####

 # The nest option in tidyr can expand on the group_by functions
 # group_by() essentially creates a dataframe (actually, tibble) for each group value
 # nest(), in combination with group_by(), makes the smaller tibbles available to us as a list
 #  this can be extremely handy for downstream analyses

 # Example with R's gapminder package and associated dataset:
 # install.packages("gapminder")
library(gapminder)

# Let's take a brief look at the dataset using dplyr::glimpse() function
library(tidyverse)
dplyr::glimpse(gapminder)

# If we group by continent and then nest, we get a tibble with each row as a continent
# and a smaller tibble with other variables corresponding to the continent in a list
nested_by_continent <- gapminder %>%
  group_by(continent) %>%
  nest()
nested_by_continent

# If you want to access the first column, use double brackets
nested_by_continent[['continent']]

# To access the first continent's dataframe:
nested_by_continent[['data']][[1]]

# Let's say we want to run a linear model:
fit <- lm(lifeExp ~ gdpPercap, 
          data=nested_by_continent[['data']][[1]])

# Let's make it a function to take in any dataset that we feed in, rather than specifying, as above
le_vs_gdpPercap <- function(df) {
  lm(lifeExp ~ gdpPercap, data = df)
}

# The function above let's us generalize building the lm, and we can use the function on
# each one of our data frames in the list in the nested object

# Instead of using a for loop to go through each smaller dataframe, we can use purrr::map()
# to do the work for us

#install.packages("purrr")
library(purrr)

# The following will map through the first two continent's dataframes
map(nested_by_continent$data[1:2], le_vs_gdpPercap)

# Now we can combine map() with nest() and mutate() to add the results of the linear model
# for each continent to a tibble, just like we added smaller tibbles previoulsy

# We already have the tibble from nest() and we can feed that to mutate(), where we map data
# to the linear model function, written above

update <- nested_by_continent %>%
  mutate(model = purrr::map(data, le_vs_gdpPercap))
update
# This results in a tibble with the extra column containing the linear model objects, which
# are available as a list
# And just like that, we ran 5 linear regressions and saved the model objects

# Let's say we want to add a column of the slopes
# Write function to get slope from the fit
get.slope <- function(model){
  model$coefficients[2]
}

# Use the function to create new column that contains slope in the nested_by_continent dataframe
update.slopes <- update %>%
  mutate(slope = purrr::map(model, get.slope))  #fit here refers to the column "fit" in the nested df
update.slopes
#### 3b. Alternative to nest() ####
# Another option is to split() the dataset into pieces, fit a model to each piece, summarise
# and extract

gapminder %>%
  split(.$continent) %>% #note the .$ here is needed because split is not a tidyr function
  map(~ lm(lifeExp ~ gdpPercap, data = .x)) %>%
  map(summary) %>%
  map_dbl("r.squared")
