#### FWRI R Club: January 14, 2020 ####

#### 1. Announcements ####

# Short meeting today because seminar at 2:30 pm in KAS (FWRI)

# Tampa Bay R Users Group Meetup: Jan. 21, 7 - 9 pm at Southern Brewing and Winemaking
  # Topic = Travis Gerke presents: Job market analysis for the R data scientist/unicoRn!

# R-Ladies Tampa: Sat., Feb. 8, 10a - 12p at the Entrepreneur Collaborative Center, Tampa
  # This will be the start of a set of R-Ladies meet-ups that go through the Tidyverse
  # They'll start with Data Visualization this month
  # Next month will be data transformation/wrangling

# Our Next Meeting: February 11, 2020 

#### 2. Code Snippets ####
# https://appsilon.com/r-studio-shortcuts-and-tips-part-2/

     # good for when writing same chunks over and over,
     # or remembering all the brackets or
     # required parameters for functions
     
     # How to use:
     # Recognized by auto-completion by a {snippet} tag
     # just type the snippet name and hit "Tab" twice to use it
     # Some of the snippets which are available by default include:
        # Declarations: lib, req, fun, ret, mat
        # Loops: for, while, switch
        # Conditionals: if, el, and ei for conditionals
        # Apply family functions: apply, lapply, sapply, etc.
        # S4 classes/methods definitions: sc, sm, and sg.
        # Shiny App template: shinyapp

     # Example of "if": type "if" and hit tab


     # Example of for loop: type "for" and hit tab


     # Example of inserting a comment code line with a time stamp: type "ts"


     # For customizing or creating your own snippets use Edit Snippets button under Snippets section in
     # Tools -> Global Options -> Code (see link at beginning of Section 2 for details)


#### 3. Editing with Multiple Cursors ####
# https://www.r-bloggers.com/rstudios-multiple-cursors-rule/
     # In RStudio, write and edit in more than one place at once
     # Can save you time with repetitive tasks that span multiple lines

     # Example: I'm often asked to check out some information for multiple "samples" or "references"
     #          and they are character strings. I'd like to be able to take a list of them from Excel,
     #          or an email or whatever I've been sent so I can make it into a character vector

     # There's a list of "references" below to try it out:
       # hold down Alt and click and drag the cursor down the left side of each row (or Ctrl + Alt + Down)
       # to highlight: hold down Alt+Shift, then press the right arrow
       # add quotes: use Shift+' (the quotation mark key)
       # while cursor still across all lines at the right side, press the comma key to add a comma
       # clean it up by removing the last column and passing to the c() function and continue on

#Try it out for yourself here
TBM2014111006
TBM2015010604
TBM2015011101
TBM2015040406
TBM2015050201
TBM2015051004
TBM2015061601

#Here's an example of what you can do with the quoted string:
RefCheck  <- c("TBM2014111006",
"TBM2015010604",
"TBM2015011101",
"TBM2015040406",
"TBM2015050201",
"TBM2015051004",
"TBM2015061601")
RefCheck