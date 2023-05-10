library(rhdf5)
library(tidyverse)

datafile <- '/home/gummi/banditExperiment/dataset.data'

# Create an empty dataframe
df <- tibble(Mouse=character(), Sex=character(), Cohort=character(), State=character(), Day=numeric(), Pump=character(), Prob=numeric())

# Load the hdf5 file
file <- h5readFile(datafile)

# Create a dataframe 
for(mouse in names(file[[cohort]])){  
  # Determine sex of mouse
  sex <- h5read(file, paste0(cohort, '/', mouse, '/attrs/Sex'))
  
  # Sort experiment days
  exp_days <- names(file[[cohort]][[mouse]][[state]])
  exp_days <- sort(exp_days)

  for(ei in seq(length(exp_days))){
    exp <- exp_days[ei]
    expgrp <- file[[cohort]][[mouse]][[state]][[exp]]
    
    choice <- as.character(expgrp$Action$Choice)
    
    pump_counts <- table(choice)
    probs <- as.numeric(pump_counts) / sum(pump_counts)
    pumps <- names(pump_counts)
    
    for(i in seq(along=pumps)){
      df <- rbind(df, tibble(Mouse=mouse, Sex=sex, Cohort=cohort, State=state, Day=ei, Pump=pumps[i],Prob=probs[i]))
    }
  }
}

# Plotting
ggplot(df, aes(x=Pump, y=Prob)) +
  geom_bar(stat="identity", fill="grey", width = 0.5) +
  geom_line(aes(group=Mouse), color='black', size=0.25, alpha=0.25) +
  geom_hline(aes(yintercept=1/3), color='crimson', linetype='dashed', size=0.5) +
  theme_minimal() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylim(0, 0.5) +
  ggtitle("Deterministic") +
  theme(plot.title = element_text(hjust = 0.5, size=8), axis.text.y = element_text(size=8))

