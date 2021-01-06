# A mini lesson on dplyr from the R For Data Science book by Hadley Wickham 

# Prepped by L. Herdter

trace(utils:::unpackPkgZip, edit=TRUE)
library(nycflights13)
library(tidyverse)

#flights dataset within the nycflights13 package
flights 

#dplyr Basics

# - pick observations by values using filter()
# reorder rows with arrange()
# pick variables by their names with select()
# create new variables that are functions of existing variables with mutate()
# collapse many values in to a single summary with summarize()

# all of these can be used in conjunction with group_by() by changing the scope of the function so they operate on a group
# by group basis

#combine multiple operations with the Pipe 

# FILTER ####
#flights on jan1
filter(flights, month==1, day==1)

#flights in nov and december
filter(flights, month %in% c(11,12))

#flights that were not delayed more than 120 minutes
filter(flights, arr_delay <=120, dep_delay <=120)

#ARRANGE ####
arrange(flights, arr_delay)

#arrange by descending
arrange(flights, desc(arr_delay))

#arrange to sort NAs
names(flights)[colSums(is.na(flights)) >0]

arrange(flights, desc(is.na(dep_time)), desc(is.na(dep_delay)))

#SELECT ####
select(flights, year, month, dep_delay, everything())

select(flights, year:day)

select(flights, -(year:day))

#additional helper functions with select
#starts_with('abc')
#ends_with()
#contains()
#matches(regex)
#num_range("x", 1:3)
#everything()

#select can be used to rename varaibles but its not terribly useful to do so because it will drop all of the variables 
# not explicitly mentioned
# use rename()
rename(flights, tail_num=tailnum)

#everything is helpful if you have some varaibles youd like to move to the start of thedataframe
select(flights, time_hour, air_time, everything())

#MUTATE ####
#adds new columns to dataset

flights_sml <- select(flights, year:day, ends_with("delay"), distance, air_time)

#can refer to ones you just created within the mutate statement
mutate(flights_sml, gain=arr_delay-dep_delay, speed=distance/air_time*60, hours= air_time/60, gain_per_hour=gain/hours )

#if you only want to keep new values use transmute()

flights_trans <- transmute(flights, gain=arr_delay - dep_delay, hours=air_time/60, gain_per_hour=gain/hours)

#all sorts of arithmetic operations including modular arithmetic (to calculate remainders)

#SUMMARIZE ####
#grouped summaries
by_day <- group_by(flights, year, month, day)
summarize(by_day, delay=mean(dep_delay, na.rm=TRUE))
# many summary functions
# mean, median, sd, IDR, mad, min, quantile, max, first, nth, last, counts with n_distinct, count with length() or just with count()

#PIPE ####
by_day <- group_by(flights, year, month, day) %>% summarize(delay=mean(dep_delay, na.rm=TRUE)) %>% 
  mutate() %>% select
  


# extra tidbits-- always useful to ungroup after grouping 
# can pipe a ggplot right onto a dplyr statement

not_cancelled <- flights %>% filter(!is.na(dep_delay), !is.na(arr_delay))

#measure delays
delays <- not_cancelled %>% group_by(tailnum) %>% summarize(delay=mean(arr_delay, na.rm=TRUE), n=n()) 

delays %>% ggplot(mapping=aes(x=n, y=delay)) + geom_point(alpha = 1/10)
