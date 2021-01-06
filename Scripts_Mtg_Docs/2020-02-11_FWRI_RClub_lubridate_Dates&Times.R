#### FWRI R Club: February 11, 2020 ####

#### 1. Announcements ####
  # Tampa Bay R Users Group Meetup: Feb. 25, 7 - 9 pm at Southern Brewing and Winemaking
  # Topic = creating presentations in R that look good with xaringanExtra

# R-Ladies Tampa: Sat., Mar. 13, 10a - 12p at the Entrepreneur Collaborative Center, Tampa
  # Topic = Tidyverse Skills: Data Transformation

# Our Next Meeting: March 10, 2020 at 2 pm 

#### 2. Dates and Times ####
 #Information comes from the book "R for Data Science", Chapter 16 and focuses on the "lubridate" package
 #Base-R has multiple ways to deal with dates and times too, but we will focus on "lubridate" today

 # Dates are just...complicated. They deal with the rotation of the Earth and it's orbit around the sun
 # Then we need to throw in months, time zones, daylight savings time, etc.
 # This is not all-encompassing but should help you get started

# Let's load the libraries we need

library(tidyverse)
library(lubridate)
library(nycflights13) #This is an R dataset, accessible to all R-users after you install the package

 # We'll focus on dates and date-times. If you need to storing times, try the "hms" package

# To get the current date, you can use today() or now()
today()
now()

#### Creating Dates and Times ####
# There are generally 3 ways you'll create a date/time
  # from a string
  # from individual date-time columns or components
  # from an existing date/time object from your database or data file

# From a string
  # Use the helpers provided by lubridate and 
  # identify the order of year, month, and day in your dates;
  # then arrange "y", "m" and "d" in same order
ymd("2017-01-31")

mdy("January 31st, 2017")

dmy("31-Jan-2017")

  # These functions also take unquoted numbers!
ymd(20170131)


  # ymd() and others create dates with times. Add an underscore
  # to one or more of "h", "m" and "s' to the name of the function
ymd_hms("2017-01-31 20:11:59")

mdy_hm("01/31/2017 08:01")

  # if needed, you can force a date-time by supplying a time zone
ymd(20170131, tz = "UTC")


# From individual date-time columns or components
  # Useful when your date and time components are spread across multiple columns
  # To create a date-time, use make_date() for dates or make_datetime() for date-times
flights %>% 
  select(year, month, day, hour, minute) %>% 
  mutate(departure = make_datetime(year, month, day, hour, minute))

# From existing date/time object
  # If you want to switch between date-time and date use
  # as_datetime() and as_date()
as_datetime(today())

as_date(now())

#### Extracting Components of Dates and Times ####
  # Pull out individual parts of the date with accessor functions
    # year()
    # month()
    # mday() # day of the month
    # yday() # day of the year - Julian day
    # wday() # day of the week
    # hour()
    # minute()
    # second()

datetime <- ymd_hms("2020-02-11 08:47:30", tz = "EST")
datetime

yr <- year(datetime)

month(datetime)

mday(datetime)

yday(datetime)

wday(datetime)

  # For month() and wday() you can set label = TRUE to return the month name abbreviation or day of week
    # set abbr = FALSE to return the full name

month(datetime, label = TRUE)

wday(datetime, label = TRUE, abbr = FALSE)

#### Rounding Dates ####
  # This can be useful for plotting, or calculating difference between a date and a rounded date/group
  # To round the date to a nearby unit of time, use floor_date(), round_date(), and ceiling_date()
  # Each function needs a vector of dates and the name of the unit to round down, up, or to

# Example: If you wanted to plot the number of flights per week

#First we need to make the time a little different to include the delay and pull out the hour and minute components
#It's just stored in a weird way in the flights dataset - the original flights database does not have an actual departure
 #time as we think of time
#This uses "modular math", but basically we are interested in using the flights_dt dataset below
make_datetime_100 <- function(year, month, day, time) {
  make_datetime(year, month, day, time %/% 100, time %% 100)
}

flights_dt <- flights %>% 
  filter(!is.na(dep_time), !is.na(arr_time)) %>% 
  mutate(dep_time = make_datetime_100(year, month, day, dep_time)) %>%
  select(origin, dest, ends_with("delay"), ends_with("time"))

flights_dt2 <- flights_dt %>%
  count(week = floor_date(dep_time, "week")) 

# If we want to plot the number of flights by week
flights_dt2 %>% 
  ggplot(aes(week, n)) +
  geom_line()


#### Setting Components ####
 # In case you need to modify parts of your date...
 # You can use the same accessory functions and assign the components

# Remind ourselves of what we set as our datetime for today's exercise
datetime

year(datetime) <- 2019
datetime

month(datetime) <- 01
datetime

hour(datetime) <- hour(datetime) + 1
datetime

# Rather than modifying in place, you can create a new date-time
  # with update()
update(datetime, year = 2020, month = 2, mday = 11, hour = 3)


#### Time Spans ####
     # durations (exact number of seconds)
     # periods (human units like weeks, months)
     # intervals (start and end point)

# When you subtract two dates you get a difftime object

work_years <- today() - ymd(20161026)
work_years

#This difftime ojects records a time span of seconds, minutes, hours, days, or weeks
# This can be difficult to work with so you can use as.duration() to always get seconds

as.duration(work_years) # notice  how the output here also gives you an approximation in years in this case

#Durations come with a bunch of constructors:
dseconds(15)
#> [1] "15s"
dminutes(10)
#> [1] "600s (~10 minutes)"
dhours(c(12, 24))
#> [1] "43200s (~12 hours)" "86400s (~1 days)"
ddays(0:5)
#> [1] "0s"                "86400s (~1 days)"  "172800s (~2 days)"
#> [4] "259200s (~3 days)" "345600s (~4 days)" "432000s (~5 days)"
dweeks(3)
#> [1] "1814400s (~3 weeks)"
dyears(1)
#> [1] "31536000s (~52.14 weeks)"

#You might get some unexpected results with durations when it comes to adding/subtracting dates and times,
# so use of periods can be helpful because they don't have a fixed timespan in seconds

one_pm <- ymd_hms("2016-03-12 13:00:00", tz = "America/New_York")

one_pm
#> [1] "2016-03-12 13:00:00 EST"
one_pm + days(1)
#> [1] "2016-03-13 13:00:00 EDT"

#Periods also have constructors:
seconds(15)
#> [1] "15S"
minutes(10)
#> [1] "10M 0S"
hours(c(12, 24))
#> [1] "12H 0M 0S" "24H 0M 0S"
days(7)
#> [1] "7d 0H 0M 0S"
months(1:6)
#> [1] "1m 0d 0H 0M 0S" "2m 0d 0H 0M 0S" "3m 0d 0H 0M 0S" "4m 0d 0H 0M 0S"
#> [5] "5m 0d 0H 0M 0S" "6m 0d 0H 0M 0S"
weeks(3)
#> [1] "21d 0H 0M 0S"
years(1)
#> [1] "1y 0m 0d 0H 0M 0S"

# Periods can be very helpful for things like Daylight Savings Time
# Daylight Savings Time

one_pm
one_pm + ddays(1)
#> [1] "2016-03-13 14:00:00 EDT"
one_pm + days(1)
#> [1] "2016-03-13 13:00:00 EDT"


# Intervals
  # An interval is a duration with a starting point
  # This makes it precise so you can determine exactly
   #how long it is

years(1) / days(1)
next_year <- today() + years(1)
(today() %--% next_year) / ddays(1)
#> [1] 366

#### Time Zones ####
# In R, time zones are an attribute of the date-time that
  # only controls printing
# The three objects below represent the same instant in time

(x1 <- ymd_hms("2015-06-01 12:00:00", tz = "America/New_York"))
#> [1] "2015-06-01 12:00:00 EDT"
(x2 <- ymd_hms("2015-06-01 18:00:00", tz = "Europe/Copenhagen"))
#> [1] "2015-06-01 18:00:00 CEST"
(x3 <- ymd_hms("2015-06-02 04:00:00", tz = "Pacific/Auckland"))
#> [1] "2015-06-02 04:00:00 NZST"

# Verify they are the same with subtraction:
x1 - x2
#> Time difference of 0 secs
x1 - x3
#> Time difference of 0 secs

# Unless otherwise specified, "lubridate" always uses UTC, which is the standard time used by the scientific community
  # UTC does not have DST, which  makes it convenient for computation
  # Operations like c() will often drop the time zone; date-times will display in your local time zone
x4 <- c(x1, x2, x3)
x4
#> [1] "2015-06-01 12:00:00 EDT" "2015-06-01 12:00:00 EDT"
#> [3] "2015-06-01 12:00:00 EDT"

# If you need to change the time zone
  # Option 1: keep the instant in time the same but change how displayed
x4a <- with_tz(x4, tzone = "Australia/Lord_Howe")
x4a
#> [1] "2015-06-02 02:30:00 +1030" "2015-06-02 02:30:00 +1030"
#> [3] "2015-06-02 02:30:00 +1030"
x4a - x4
#> Time differences in secs
#> [1] 0 0 0

  # Option 2: change the underlying instant in time (when an instant has been labelled with incorrect time zone)
x4b <- force_tz(x4, tzone = "Australia/Lord_Howe")
x4b
#> [1] "2015-06-01 12:00:00 +1030" "2015-06-01 12:00:00 +1030"
#> [3] "2015-06-01 12:00:00 +1030"
x4b - x4
#> Time differences in hours
#> [1] -14.5 -14.5 -14.5
 