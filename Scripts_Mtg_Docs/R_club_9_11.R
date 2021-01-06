# TIDY DATA- from R for Data Science by Hadley Wickham####
# with tidyr()
# tidy data is required for many of tidyverse functionality including ggplot, dplyr, etc. 
#http://r4ds.had.co.nz/tidy-data.html

#tidy data is not end all be all- 
# Jeef Leeks blog post about the merits of 'untidy' data
# https://simplystatistics.org/2016/02/17/non-tidy-data/
#  "the goal is to solve a particular problem and the format I chose is the one that makes it most direct/easy to solve that problem, 
# rather than one that is theoretically optimal. "
  
library(tidyverse)

#the same dataset four different ways showing four variables: country, year, population, cases

tidyr::table1  #this is the only tidy dataset here

tidyr::table2

tidyr::table3

tidyr::table4a

tidyr::table4b

#Rules for a tidy dataset####
#1. Each variable has its own column
#2. Each observation has its own row



#How to make a tidy dataset?####
#Gathering and Spreading
# equivalent to melt and cast in baseR
# melt is to make something long as is gathering- makes wide tables narrow and long
# cast will make something wide, so will spreading- makes long tables shorter and wider

#GATHERING ####

tidyr::table4a # to tidy this dataset we need to GATHER columns into a new pair of variables

# three parameters to gather
# 1. the set of columns that represent VALUES  in the CURRENT layout -> here the columns are 1999 and 2000
#       -> columns to gather are specified with dplyr:select() style
#       -> c(Base, M_precip, M_temp)
# 2. the name of the variable whose values will form the column names in the NEW layout, called key
# 3. the name of the variable whose values are spread over the cells in the CURRENT layout, called value

table4a %>% gather('1999', '2000', key='year', value='cases')

# final result is that the gathered columns get dropped and we get the new key and value columns

#SPREADING ####
#opposite of gathering- used when observations are scattered across multiple rows

tidyr::table2

# two parameters to spread
#1. column that contains VARIABLE names in the CURRENT layout, called key
#2. column that cotntains values from multiple variables in the CURRENT layout, called value

spread(table2, key=type, value=count)

#options
?spread

# SEPARATING/PULLING####
#separate() and unite()

tidyr::table3
#one column rate has two varaibles (cases and population)

#seperate() ####
#similar to substr() in baseR
#pulls apart one column into multiple columns by splitting wherever a separator character appears
# dy default it splits values wherever it sees a nonalpha numeric character

#three variations of separate
#1. by separator with separate()
table3 %>% separate(rate, into=c('cases', 'population'))

#can use a specific character to split- sep
table3 %>% separate(rate, into=c("cases", "population"), sep="/")


#2. by position with separate()
#can also pass a vector of integers to separate.sep- interprets it as a positions to split at
# positive values start at 1 on the far left of the strings

table3 %>% separate(year, into=c("century", "year"), sep=2, convert=TRUE)

#3. by groups with extract()
#extract()- extracts one column into multiple columns
?extract
df <- data.frame(x = c(NA, "a-b", "a-d", "b-c", "d-e"))
df
df %>% extract(x, "A") #extract(df, x, "A")

df %>% extract(x, c("A", "B"), "([[:alnum:]]+)-([[:alnum:]]+)")

#other tidbits about separate
#use convert to convert to better types so cases and population are integers
table3 %>% separate(rate, into=c("cases", "population"), convert=TRUE)

#unite() ####
#inverse of separate()- combines multiple columns into a single column
#takes a dataframe, name of the new column and the columns to join

tidyr::table5


table5 %>% unite(new, century, year) 

# all of these can take a dataframe in the traditional sense w/o piping operator
unite(table5, new, century, year)

# default will place an underscore (_)
#if you dont want a separator

table5 %>% unite(new, century, year, sep="")

?unite

table5 %>% unite(new, century, year, sep="", remove=FALSE)

# Missing Values ####
#1. explicitly missing - flagged with NA
#2. implicitly missing- just not present in the data

stocks <- tibble(
  year= c(2015, 2015, 2015, 2015, 2016, 2016, 2016),
  qtr=c(1,2,3,4,2,3,4),
  return=c(1.88, 0.59, 0.35, NA, 0.92, 0.17, 2.66)
)

stocks

# explicit- return for 4th qtr in 2015 indicated by NA
# implicit- return for 1st qtr in 2016- just not there

#can make IMPLICT to EXPLICIT in a tidy format
#with complete() - takes a set of columns and finds all unique combinations and ensures the original dataset
# contains all values and fills in explicit NAs where necessary
stocks %>% complete(year, qtr)

#make EXPLICIT TO IMPLICIT in a tidy format
# two step process
# first untidy
# then retidy
stocks %>% spread(year, return) %>%
  gather('2015', '2016', key= "year", value="return", na.rm=TRUE)




