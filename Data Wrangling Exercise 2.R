#Data Wrangling Exercise 2
#Zach Young

install.packages("tidyr")
install.packages("dplyr")
library(tidyr)
library(dplyr)

# First, I changed all blank values to NA
dat <- read.csv("titanic_original.csv", na.strings = " ", header = TRUE)

# Port of Embarkation
titanic_original$embarked[is.na(titanic_original$embarked)] <- "S"

# Filling in missing age
titanic_original$age[is.na(titanic_original$age)] <- titanic_original %>% 
                                                     select(age) %>%
                                                     summarise(avg = mean(age, na.rm = TRUE))

# Filling in Lifeboat with "None"
titanic_original$boat[is.na(titanic_original$boat)] <- "None"

# Create new column has_cabin_number
titanic_original1 <- mutate(titanic_original, has_cabin_number = ifelse(cabin %in% NA, yes = 0, no = 1))

titanic_clean <- titanic_original1




