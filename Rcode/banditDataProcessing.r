# This file contains functions that are used to import and process the mouse bandit expiriment data
library(tidyverse)


# Onetime function for loading all the excel files into a list and creating a csv with the stats

loadExcelLocFileGenerateCSV <- function(folder_path) {
    files <- list.files(path = folder_path, pattern = "*.xlsx", full.names = TRUE)
    summaryData <- list()
    count <- 1
    for(file_path in files){
            sheet <- read_excel(file_path, col_names = FALSE,progress = FALSE)
            splitt <- str_split(as.character(sheet[19, 2]), "\\\\", simplify = TRUE)
            file_name <- splitt[length(splitt)]
            file_name <- str_remove(file_name, ".avi")
            parts <- str_split(file_name, "-", simplify = TRUE)
            date <- paste0(substr(parts[[2]], 3,4),substr(parts[[2]], 1,2),substr(parts[[2]], 5,8))
            string <- paste0(parts[[1]],"_",date)
            mazeType <- as.character(sheet[10, 2])
            start <- as.numeric(sheet[1, 2])
            # Process data
            data <- as.data.frame(sheet)
            colnames(data) <- data[start-1,]
            data <- data[-c(1:(start)),]
            data$MazeType <- mazeType
            
            # Save the data as a csv
            write.csv(data,file = paste0("RawTelemteryData/",string,".csv"))

            # Filter out the the first 300 seconds
            data <- data[which(as.numeric(data[["Recording time"]]) > 300),]

            # Compute total distance traveled
            totalDistance <- sum(as.numeric(data[["Distance moved"]]),na.rm = TRUE)

            # Compute the mean speed
            meanSpeed <- mean(as.numeric(data[["Velocity"]]),na.rm=TRUE)

            # Write data to list
            summaryData[[count]] <- c(string,totalDistance,meanSpeed,mazeType)
            count <- count + 1

    }
    # Convert the list to a dataframe
    summaryData <- do.call(rbind,summaryData)

    # Rename the columns
    colnames(summaryData) <- c("FileName","TotalDistance","MeanSpeed","Maze")

    # Write the data to a csv
    write.csv(summaryData,file = "TelemetryStatsProbPhase.csv")

  return(TRUE)
}


processLocationData <- function(experiment,attribute) { # nolint

        telemetryDataframe <- experiment[["LocationData"]]

        # Compute total distance traveled
        totalDistance <- sum(as.numeric(telemetryDataframe[["Distance moved"]]),na.rm = TRUE) %>% 
                                                                        addPrefix("totalDistance")
        # Compute the mean speed
        meanSpeed <- mean(as.numeric(telemetryDataframe[["Velocity"]]),na.rm=TRUE) %>% 
                                                                        addPrefix("meanSpeed")

        # Compute the distance reward gotten per distance traveled 
        rewardPerDistance <- sum(as.logical(experiment[["Reward"]]))/unname(totalDistance) %>% 
                                                                        addPrefix("rewardPerDistance")

        # Compute the distance from each station
        distsDF <- computStationDistance(telemetryDataframe,attribute)

        # Compute the mean distance from each station
        meanDists <- colMeans(distsDF[,c(1:3)],na.rm = TRUE) %>%
                                       addPrefix("meanDistanceFrom_")

        # Plot a line plot of the distances through time 
        # ggplot(pivot_longer(distsDFlowRes,cols = c("0.25","0.5","1"),names_to = "variable"), aes(x = Time, y = value)) + 
        #    geom_line(aes(color = variable), size = 1) +
        #    scale_color_manual(values = c("#00AFBB", "#E7B800","red")) +
        #    theme_minimal()

        # Compute the binned per minute avarage and plot it as line plot
        #distsDFlowRes %>%
        #   mutate(interval = Time%/%60) %>%
        #    group_by(interval) %>%
        #    summarise(mean = mean()) %>% 

        #print(paste0("Total Distance: ",totalDistance))	
        #print(paste0("Mean Speed: ",meanSpeed))
        #print(paste0("Mean Distance from 0.25: ",meanDists[["meanDistanceFrom_0.25"]]))
        #print(paste0("Mean Distance from 0.5: ",meanDists[["meanDistanceFrom_0.5"]]))
        #print(paste0("Mean Distance from 1: ",meanDists[["meanDistanceFrom_1"]]))

        return(c(totalDistance,meanSpeed,rewardPerDistance,meanDists))

}



# Function that converts the location data from x and y coords to distance from the different stations

computStationDistance <- function(telemetryDataframe,attribute) {

    # Map the maze type to the station locations    A     C      B                     A    C       B
    StationLocations <- list("Upright Y maze" = c(0,-34, 28,15.5, -32,18.6), "CCW Y maze" = c(41,3, -21,34, -17,-34)) 

    # Get the location vector for the current maze
    LocVec <- StationLocations[[telemetryDataframe[["MazeType"]][[1]]]]

    # Compute the distance from each station
    statA <- sqrt((as.numeric(telemetryDataframe[["X center"]]) - LocVec[1])^2 + (as.numeric(telemetryDataframe[["Y center"]]) - LocVec[2])^2)
    statC <- sqrt((as.numeric(telemetryDataframe[["X center"]]) - LocVec[3])^2 + (as.numeric(telemetryDataframe[["Y center"]]) - LocVec[4])^2)
    statB <- sqrt((as.numeric(telemetryDataframe[["X center"]]) - LocVec[5])^2 + (as.numeric(telemetryDataframe[["Y center"]]) - LocVec[6])^2)

    # Create a dataframe with the distances
    distsDF <- data.frame(statA,statB,statC) 

    if("0.25" %in% as.character(c(attribute$A_Prob, attribute$B_Prob, attribute$C_Prob)) ) {
            colnames(distsDF) <- as.character(c(attribute$A_Prob, attribute$B_Prob, attribute$C_Prob)) 
            # Reorder columns in the dataframe
            distsDF <- distsDF[,c("0.25","0.5","1")]
    } else {
        colnames(distsDF) <- c("A","B","C")
    }
    


    distsDF$Time <- as.numeric(telemetryDataframe[["Recording time"]])

    return(distsDF)
}

# Create synthetic data
syntheticData <- function(nTrials, LeftBias, ProbBandicVec, modeProbVec, timePerBandit, RewardProbPerStation) {
    
    # ProbBandicVec <- setNames(c(0.5,0.5,0.5),c("AC>A","AB>A","BC>B"))
    # LeftBias <- 0.5
    # LeftBiasVsBanditProb <- 0.5
    # RewardPerStation <- setNames(c(1,0.5,0.25),c("A","B","C"))

    experimentObj <- list()
    # Create the bandit sequence
    Bandit <-  c("ABC") # First bandit is always ABC
    Choice <- c(sample(c("A","B","C"),1)) # Randomly sample the first choice with equal probability	
    Time <- c(ceiling(runif(1)*10)) # Randomly sample the first time (0-10 seconds)
    Reward <- c(sample(c(TRUE,FALSE),1,prob = c(RewardProbPerStation[Choice[1]],1-RewardProbPerStation[Choice[1]]))) # Sample the first reward according to the first choice
    for(i in 2:nTrials){
        # Figure out which bandit we are in
        if(Choice[i-1] == "A"){
            Bandit[i] <- "BC"	
        }else if(Choice[i-1] == "B"){
            Bandit[i] <- "AC"
        }else if(Choice[i-1] == "C"){
            Bandit[i] <- "AB"
        }
        mode <- sample(c("LeftBias","Optimal","Explorer"),1,prob = modeProbVec)

        # Sample the pump choice either according to left bias or the bandit with some probability
        if(mode == "Optimal") {
            # Sample the pump choice according to the probility of the bandits
            if(Bandit[i] == "AB"){
                Choice[i] <- sample(c("A","B"),1,prob = c(ProbBandicVec["AB>A"],(1-ProbBandicVec["AB>A"])))
            }else if(Bandit[i] == "AC"){
                Choice[i] <- sample(c("A","C"),1,prob = c(ProbBandicVec["AC>A"],(1-ProbBandicVec["AC>A"])))
            }else if(Bandit[i] == "BC"){
                Choice[i] <- sample(c("B","C"),1,prob = c(ProbBandicVec["BC>B"],(1-ProbBandicVec["BC>B"])))
            }
        } else if (mode == "LeftBias") {
            # Sample the pump choice according to the left bias
            if(Bandit[i] == "AB"){
                Choice[i] <- sample(c("A","B"),1,prob = c(LeftBias,(1-LeftBias)))
            } else if (Bandit[i] == "AC"){
                Choice[i] <- sample(c("A","C"),1,prob = c((1-LeftBias),LeftBias))
            } else if (Bandit[i] == "BC"){
                Choice[i] <- sample(c("B","C"),1,prob = c(LeftBias,(1-LeftBias)))
            } 
         } else {
          # Find which choice has been made the least
          #choiceCount <- table(factor(Choice) , levels = c("A","B","C"))
          ## Depending on what the current bandit is, sample the choice with the lowest count
          #if(Bandit[i] == "AB"){
          #  Choice[i] <- names(which.min(choiceCount[1:2]))
          #} else if(Bandit == "AC"){
          #  Choice[i] <- names(which.min(choiceCount[c(1,3)]))
          #} else if(Bandit == "BC"){
          #  Choice[i] <- names(which.min(choiceCount[2:3]))
          #}
            # Sample the pump choice according inverse of GO
            if(Bandit[i] == "AB"){
                Choice[i] <- sample(c("A","B"),1,prob = c(1-ProbBandicVec["AB>A"],ProbBandicVec["AB>A"]))
            } else if(Bandit[i] == "AC"){
                Choice[i] <- sample(c("A","C"),1,prob = c(1-ProbBandicVec["AC>A"],ProbBandicVec["AC>A"]))
            } else if(Bandit[i] == "BC"){
                Choice[i] <- sample(c("B","C"),1,prob = c(1-ProbBandicVec["BC>B"],ProbBandicVec["BC>B"]))
            }

         }

        # Sample the time according to the bandit/choice combo and the banditTimeVec with std dev of 2
        Time[i] <- Time[i-1] + rnorm(1,timePerBandit[paste0(Bandit[i],"->",Choice[i])],2)

        # Sample the Reward probiltiy according to the choice made
        Reward[i] <- sample(c(TRUE,FALSE),1,prob = c(RewardProbPerStation[Choice[i]],1-RewardProbPerStation[Choice[i]])) 
    }
    # Store the data in the experiment object
    experimentObj$Bandit <- Bandit
    experimentObj$Choice <- Choice
    experimentObj$Reward <- Reward
    experimentObj$Time <- Time

    # Create the list of attributes
    attribute <- list()
    attribute$A_Prob <- RewardProbPerStation["A"]
    attribute$B_Prob <- RewardProbPerStation["B"]
    attribute$C_Prob <- RewardProbPerStation["C"]
    attribute$Duration <- Time[nTrials]
    attribute$Experiment <- "Synthetic_Data"
    attribute$Mode <- 1
    attribute$Mouse <- 9999
    attribute$Probabilistic_Switch <- 300 # 5 minutes

    return(list(experimentObj,attribute))
    
}

# Add prefix to named vecs
addPrefix <- function(x,prefix,suffix = ""){
    names(x) <- paste0(prefix,names(x),suffix)
    return(x)
}

# Function for computing cramersV for vector with a certain lagg
cramersVautoCorr <- function(x,lagg) {
    if(lagg >= length(x)){
        #print("Cramers V comparison has a lagg that is longer than the input vector, returning NA")
        return(NA)
    }

    if(length(unique(x[-c(1:lagg,length(x):(length(x)-lagg+1))])) <= 1){
        #print("Cramers V comparison only has one unique value, returning 1")
        return(1)
    }

    cramerV <- rcompanion::cramerV(x[-c(1:lagg)],lag(x,lagg)[-c(1:lagg)]) 
    return(cramerV)   
}

# Function to check if the mouse made the optimal choice given a bandit
.optimalChoice <- function(bandit) {
    # Seprate the string by the > symbol
    banditSplit <- str_split(bandit,"->")
    if(banditSplit[[1]][2] == "1"){
        return("TRUE")
    }else if(banditSplit[[1]][1] == "1" & banditSplit[[1]][2] == "0.5"){
        return("TRUE")
    }else {
        return("FALSE")}

}
# Vectorize the .optimalChoice function
optimalChoice <- Vectorize(.optimalChoice)

# Function that takes the pump choices and translates them into which bandit was chosen given the probabilities at each pump
.banditChoice <- function(pumpChoice,bandit,pumbProb_A,pumpProb_B,pumpProb_C){
    if (bandit == "AB"){	
        isAt <- pumpProb_C}
    else if (bandit == "AC"){
        isAt <- pumpProb_B}
    else if (bandit == "BC"){
        isAt <- pumbProb_A}
    else if (bandit == "ABC"){
        isAt <- NA}

    if (pumpChoice == "A"){
        return(paste0(isAt,"->",pumbProb_A))
    }else if (pumpChoice == "B"){
        return(paste0(isAt,"->",pumpProb_B))
    }else if (pumpChoice == "C"){
        return(paste0(isAt,"->",pumpProb_C))
    }else{
        return(NA)
        }
    }	

.banditChoice <- Vectorize(.banditChoice)

# Wrapper function for determining which bandit was chosen
banditChoice <- function(experiment,attribute,returnStats) {
    A_Prob <- attribute$A_Prob
    B_Prob <- attribute$B_Prob
    C_Prob <- attribute$C_Prob
    banditChoices <- .banditChoice(experiment[["Choice"]], experiment[["Bandit"]],A_Prob,B_Prob,C_Prob)

    # If first string starts with NA, remove it
    if(str_starts(banditChoices[1],"NA")){
        banditChoices <- banditChoices[-1]
        } 

    if(!returnStats){
        return(banditChoices)
    }
    #banditChoices <- factor(banditChoices,levels = c("0.25->0.5","0.25->1", "0.5->0.25","0.5->1","1->0.25","1->0.5"))
    banditChoiceCount <- table(banditChoices)

    # Reorder the vectors and make sure they have the all all the bandits
    banditChoiceCount <- banditChoiceCount[c("0.25->0.5","0.25->1", "0.5->0.25","0.5->1","1->0.25","1->0.5")] 

    banditChoiceCount[is.na(banditChoiceCount)] <- 0

    # Compute the probability of choosing each bandit
    banditChoiceProb <- as_tibble(banditChoiceCount) %>%
                            separate("banditChoices",c("from","to"),sep = "->",remove = FALSE) %>%
                            group_by(from) %>%
                            mutate(Prob = n/sum(n))

    # Convert tibble to a named vector
    banditChoiceProb <- setNames(banditChoiceProb$Prob,banditChoiceProb$banditChoices)

    banditChoiceProb <- banditChoiceProb[c("0.25->0.5","0.25->1", "0.5->0.25","0.5->1","1->0.25","1->0.5")] %>%
                            addPrefix("P(",")")

    banditChoiceProb[is.na(banditChoiceProb)] <- 0

    banditChoiceCount <- banditChoiceCount %>% addPrefix("banditChoiceCount_")


    # Compute stats around making the optimal choice
    # Check if the mouse made the optimal choice for each bandit in the sequence
    banditOptimalChoiceSeq <- optimalChoice(banditChoices) %>% as.logical()
    banditChoicesOptimalProb <- sum(banditOptimalChoiceSeq)/length(banditOptimalChoiceSeq) %>% 
                                addPrefix("P(optimalChoice)")
    
    # Compute the cramerV for optimal choice
    cramerOptimalChoice <- cramersVautoCorr(as.character(banditOptimalChoiceSeq),1) %>% 
                            addPrefix("cramerVOptimalChoice")
    
    # Compute ngram entropy
    entropyVec <- c()
    for(n in 2:7){
        if(length(banditChoices) < n){
            ent <- c(NA) %>% addPrefix(paste0("ngramEntropy_",n))
        } else {
            temp <- ngram::ngram(paste0(banditChoices,collapse = " "),n = n) %>% 
                ngram::get.phrasetable()
            
            ent <- (sum(temp$prop * log2(temp$prop)) * -1) / log2(length(temp$prop)) %>% addPrefix(paste0("ngramEntropy_",n))
            if(is.nan(ent)){
                ent <- c(NA) %>% addPrefix(paste0("ngramEntropy_",n))
            }
        }

        entropyVec <- c(entropyVec,ent)

    }

    return(c(banditChoiceCount,
            banditChoiceProb,
            entropyVec,
            banditChoicesOptimalProb,
            cramerOptimalChoice))

}


# Function for determining which direction the mouse chose
.directionChoice <- function(pumpChoice,bandit) {

    if (bandit == "AB"){
        if (pumpChoice == "A"){
            return("L") 
        }else {
            return("R")
            }
    }else if (bandit == "AC"){
        if (pumpChoice == "C"){
            return("L")
        }else {
            return("R")
            }
    }else if (bandit == "BC"){
        if (pumpChoice == "B"){
            return("L")
        }else {
            return("R")
        }

    } else {
       return(NA)
    }

}

# Wrapper function for determining which direction the mouse chose
.directionChoice <- Vectorize(.directionChoice)

biasStats <- function(experiment,attribute,returnDirectionChoiceOnly = FALSE) {
    directionChoices <- .directionChoice(experiment[["Choice"]],experiment[["Bandit"]])

    if(is.na(directionChoices[1])){
        directionChoices <- directionChoices[-1]
    }
    if(returnDirectionChoiceOnly){
        return(directionChoices)
    }

    directionChoicesCounts <- table(directionChoices)

    if(length(directionChoicesCounts) < 2){
        if(names(directionChoicesCounts) == "L") {
            leftBias <- 1
         } else {
            leftBias <- 0 }
    } else {
        leftBias <- directionChoicesCounts[["L"]]/sum(directionChoicesCounts)
    }
    leftBias <- addPrefix(leftBias,"leftBias")

    # Compute the cramersV for direction choice
    cramerDirectionChoice <- cramersVautoCorr(directionChoices,1) %>% addPrefix("cramerDirectionChoice")

    # Compute the pump bias
        # Compute the pump count choices
    pumpProbAbs <- table(experiment[["Choice"]])/sum(table(experiment[["Choice"]]))
    pumpProbAbs <- pumpProbAbs[c("A","B","C")] 

    pumpToProb <- setNames(c("A","B","C"),as.character(c(attribute$A_Prob,
                                        attribute$B_Prob,
                                        attribute$C_Prob)))

    pumpCountRel <- pumpProbAbs[pumpToProb[c("0.25","0.5","1")]] %>% 
                    setNames(c("0.25","0.5","1")) %>%
                    addPrefix("P(",")")

    pumpCountRel[is.na(pumpCountRel)] <- 0
    pumpProbAbs[is.na(pumpProbAbs)] <- 0


    pumpProbAbs <- pumpProbAbs %>% addPrefix("P(",")")
    
    return(c(leftBias,cramerDirectionChoice,pumpCountRel,pumpProbAbs))

}


# Create a function that takes in experiment and filters it for some a time window in seconds

filterSeqByTime <- function(experiment,windowStart,windowEnd) {
    if (length(experiment[["Choice"]]) == 0){
        errorCondition("Choice is empty")
    }
    if(length(experiment[["Choice"]]) != length(experiment[["Time"]])){
        errorCondition("Choice and Time must be the same length")
    }
    if(windowStart > windowEnd){
        errorCondition("windowStart must be less than windowEnd")
    }

    inds <- (as.numeric(experiment[["Time"]]) > windowStart) & (as.numeric(experiment[["Time"]]) < windowEnd)
    newExp <- experiment
    newExp[["Choice"]] <- experiment[["Choice"]][inds]
    newExp[["Bandit"]] <- experiment[["Bandit"]][inds]
    newExp[["Time"]] <- experiment[["Time"]][inds]
    newExp[["Reward"]] <- experiment[["Reward"]][inds]

    if(!is.data.frame(experiment[["LocationData"]])){
        newExp[["LocationData"]] <- NA
    } else {
        locInds <- (as.numeric(experiment[["LocationData"]][["Recording time"]]) > windowStart) & (as.numeric(experiment[["LocationData"]][["Recording time"]]) < windowEnd)
        newExp[["LocationData"]] <- experiment[["LocationData"]][inds,]
    }

    return(newExp)

}

# Function for splitting up the probabilistic part into into sections and computing the stats for each section
splitExperiment <- function(experiment,attribute,n_splits,Prob) {
    if(Prob){
        start <- attribute$Probabilistic_Switch
    } else {
       start <- 0
    }
    res <- list()
    starts <- seq(start,
                  attribute$Duration,length.out = n_splits+1)

    # Print split times
    print(paste0("Splitting experiment. Breaks at: ",paste0(starts,collapse = ",")))

    for(i in 1:n_splits){
        expTemp <- filterSeqByTime(experiment,starts[i],starts[i+1])
        res[[paste0("Split_",starts[i],"_",starts[i+1])]] <- mouseExpStats(expTemp,
                                                                  attribute,
                                                                  paste0("Split Probabalistic Experiment, Part:",starts[i],"-",starts[i+1]))

    }
    return(res)
}

# Create a function that computes the number of time a mouse goes left or right three times in row
circleStat <- function(experiment,attribute) {
    directionChoices <- biasStats(experiment,attribute,returnDirectionChoiceOnly = TRUE)

    # Check if the sequence is atleast 3 long else just return na
    if(length(directionChoices) > 2) {
        countCircles <- 0
        for (i in 1:(length(directionChoices)-2)){
            if (directionChoices[i] == directionChoices[i+1] & directionChoices[i+1] == directionChoices[i+2]){
                countCircles <- countCircles + 1
            }
        }
        ProbOfCircle <- countCircles/(length(directionChoices)-2) %>% addPrefix("P(Circle)")


    
        #Count number of times the mouse changes direction
        countDirectionChanges <- 0
        for(i in 1:(length(directionChoices)-1)){
                if (directionChoices[i] != directionChoices[i+1]){
                    countDirectionChanges <- countDirectionChanges + 1
                }
            }
    
        probDirectionChanges <- countDirectionChanges/length(directionChoices) %>% addPrefix("P(DirectionChange)")
        countCircles <- addPrefix(countCircles,"circleCount")
        countDirectionChanges <- addPrefix(countDirectionChanges,"countDirectionChanges")
    } else {
        countCircles <- c(NA) %>% addPrefix("countCircles")
        ProbOfCircle <- c(NA) %>% addPrefix("P(Circle))")
        probDirectionChanges <- c(NA) %>% addPrefix("P(DirectionChange)")
    }

    return(c(countCircles,ProbOfCircle,probDirectionChanges))

}

# Function for determining how fast the mouse is making choices
mouseSpeed <- function(experiment,attribute) {
    banditChoices <- banditChoice(experiment,attribute,FALSE)
    time <- as.numeric(experiment[["Time"]])
    timeDiff <- diff(time)
    meanTimeToChoose <- mean(timeDiff) %>% addPrefix("meanTimeToChoice")
    if(length(timeDiff) != length(banditChoices)){
         banditChoices <- banditChoices[-1]
    }

    timePerBandit <- data.frame("banditChoices" = banditChoices,
                                 "timeDiff" = timeDiff) %>%
                                separate("banditChoices",c("from","to"),sep = "->",remove = FALSE) %>%
                                group_by(from) %>%
                                summarise(meanTimeToChoose = mean(timeDiff))

    # Convert tibble to a named vector
    timePerBandit <- setNames(timePerBandit$meanTimeToChoose,timePerBandit$from)
    timePerBandit <- timePerBandit[c("0.25","0.5","1")] %>% addPrefix("MeanTimeToChoicePerBandit_")

    return(c(meanTimeToChoose,timePerBandit))

}


rewardAttainment <- function(experiment,attribute) {
    banditChoices <- banditChoice(experiment,attribute,FALSE)

    reward <- as.logical(experiment[["Reward"]])

    totalReward <- sum(reward) %>% addPrefix("SumRewardGotten")

    rewardProb <- sum(reward)/length(reward) %>% addPrefix("P(Reward)")


    if(length(reward) != length(banditChoices)){
         reward <- reward[-1]
    }
    # Get the reward per bandit
    rewardPerBanditProb <- data.frame("choice" = banditChoices,
                                    "reward" = reward) %>%
                                    group_by(choice) %>%
                                    summarise(rewardPerBandit = mean(reward))

    # Convert tibble to a named vector and order
    ObsRewardProbPerBandit <- setNames(rewardPerBanditProb$rewardPerBandit,rewardPerBanditProb$choice)
    ObsRewardProbPerBandit <- ObsRewardProbPerBandit[c("0.25->0.5","0.25->1", "0.5->0.25","0.5->1","1->0.25","1->0.5")] %>% 
                                addPrefix("P_obs(Reward|",")")

    # Compute the observed reward per station
    ObsProbRewardPerStation <- data.frame("choice" = banditChoices,
                                    "reward" = reward) %>%
                        separate("choice",c("from","to"),sep = "->",remove = FALSE) %>%
                        group_by(to) %>%
                        summarise(rewardPerStation = mean(reward))


    # Convert tibble to a named vector and order
    ObsProbRewardPerStation <- setNames(ObsProbRewardPerStation$rewardPerStation,ObsProbRewardPerStation$to)
    ObsProbRewardPerStation <- ObsProbRewardPerStation[c("0.25","0.5","1")] %>% addPrefix("P_obs(Reward|",")")

    # Compute the probability of getting a reward after making the optimal choice
    optimalOrNot <- optimalChoice(banditChoices) %>% as.logical()
    rewardProbAfterOptimalChoic <- sum(reward[optimalOrNot])/sum(optimalOrNot) %>% addPrefix("P(Reward|OptimalChoice)")
    if(is.nan(rewardProbAfterOptimalChoic)){
        rewardProbAfterOptimalChoic <- c(NA) %>% addPrefix("P(Reward|OptimalChoice)")
    }

    # Comput how evenly the reward is distributed

    ObsProbRewardSumPerStation <- data.frame("choice" = banditChoices,
                                "reward" = reward) %>%
                    separate("choice",c("from","to"),sep = "->",remove = FALSE) %>%
                    group_by(to) %>%
                    summarise(rewardPerStation = sum(reward))
    
    rewardEveness <- sd(ObsProbRewardSumPerStation$rewardPerStation)/mean(ObsProbRewardSumPerStation$rewardPerStation) %>% addPrefix("RewardEvennes")
    if(is.nan(rewardEveness)){
        rewardEveness <- c(NA) %>% addPrefix("RewardEvennes")
    }

    return(c(totalReward,
            rewardProb,
            ObsRewardProbPerBandit,
            ObsProbRewardPerStation,
            rewardProbAfterOptimalChoic,
            rewardEveness))

}




# Function for determining some stats given a experiment object and its attribute
mouseExpStats <- function(experiment,attribute,inputDataType,verbosey = TRUE){

    if(verbosey){print(paste0("Computing experiment stats on :" ,inputDataType,". Number of Trials: ",length(experiment[["Choice"]])))}
    # Get number of trials
    numberOfTrails <- length(experiment[["Choice"]]) %>% addPrefix("NTrials")

    if(length(experiment[["Choice"]]) < 3) { return(c(numberOfTrails,rep(NA,54)))}

    # Get bandit choices stats
    banditChoices <- banditChoice(experiment, attribute, TRUE)

    # Get bias stats
    biasStats <- biasStats(experiment, attribute)

    # Retrieve Circle stats
    circleStats <- circleStat(experiment,attribute)	

    # Retrieve Mouse speed stats
    mouseSpeeds <- mouseSpeed(experiment,attribute)

    # Retrieve Reward stats
    rewardAttainment <- rewardAttainment(experiment,attribute)

    # Compute the movement stats
    if(!is.data.frame(experiment[["LocationData"]])){
        movementStats <- c(NA,NA,NA,NA,NA,NA) %>% addPrefix("movementStats_")
    } else {
        movementStats <- processLocationData(experiment,attribute)
    }

    res <- c(numberOfTrails,banditChoices,biasStats,circleStats,mouseSpeeds,rewardAttainment,movementStats)
    #print((res)) 

    return(res)
}
computeRollingPumpPreference <- function(experiment,attribute,inputDataType){
  # Get the estimate 
    banditChoices <- banditChoice(experiment, attribute, FALSE)

}



tabulateExperiments <- function(data,datafile,numberOfSplits,telemtryFolder) {
    dataType = "Action"
    # Create lists for storing the results
    extractedFeaturesDet <- list()
    extractedFeaturesProbDetPart <- list()
    extractedFeaturesProbProbPart <- list()
    extractedFeaturesRollingProb <- list()
    extractedFeaturesRollingDet <- list()
    ## Create a nested list for the split probabilistic experiments
    splitExperiments <- list()
    for(i in 1:numberOfSplits) {
        splitExperiments[[i]] <- list()
    }

    print("Starting to tabulate experiments")
    detCounter <- 0
    probDetCounter <- 0
    missExcelCounter <- 0
    missingsTelemetryData <- c()
    colnamesList <- list()
    trajFiles <- list.files(telemtryFolder)
    sequences <- list()

        # Loop over the data
        for (mousetype in names(data)) { 

        for (mouse in names(data[[mousetype]])) { 

        for (type in names(data[[mousetype]][[mouse]])) { 

        for (state in names(data[[mousetype]][[mouse]][[type]])) { 

        for (day in names(data[[mousetype]][[mouse]][[type]][[state]])) { 

            # Get the experiment
            experiment <- data[[mousetype]][[mouse]][[type]][[state]][[day]][[dataType]]
            attribute <- h5readAttributes(datafile, paste0("/",mousetype,"/",mouse,"/",type,"/",state,"/",day,"/",dataType,"/"))
            MouseAttribute <- h5readAttributes(datafile, paste0("/",mousetype,"/",mouse,"/"))

            # Add on telemetry data from csv files
            index_string <- paste0(mouse,"_",state,"_",type,"_",day)
            which(str_detect(trajFiles,index_string)) -> ind

            telemFile <- paste0(telemtryFolder,"/",trajFiles[ind[1]])
            
            print("------------------")
            print(paste0("Started processing ",mousetype," ",mouse," ",type," ",state," ",day," ",length(experiment[["Bandit"]])))

            print(paste0("Looking for telemetry file: ",telemFile))

            if(file.exists(telemFile)) {
                    experiment[["LocationData"]] <- read.csv(telemFile,check.names = F)
            } else {
                    print(paste0("No Telemetry data for ",index_string))
                    missExcelCounter <- missExcelCounter + 1
                    experiment[["LocationData"]] <- NA
                    missingsTelemetryData <- c(missingsTelemetryData,index_string)

            }

            # Check if experment is Det or Prob
            if(type == "Det") {
                
                # HACK BEWARE
                1 -> attribute$A_Prob
                0.5 -> attribute$B_Prob
                0.25 -> attribute$C_Prob

                expDetStats <- mouseExpStats(experiment, attribute,"Deterministic Experiment")
                DetData <-c("MouseType"=mousetype,
                            "Mouse"=mouse,
                            "Sex" =MouseAttribute$Sex,
                            "Type"=type,
                            "State"=state,
                            "Day"=day,
                            expDetStats)
                            
                colNamesDet <- names(DetData)
                extractedFeaturesDet[[paste0(mousetype,"_",day,"_",state,dataType,type,runif(1))]] <- unname(DetData)
                extractedFeaturesRollingDet[[index_string]] <- RollingWindowStats(experiment,attribute, c("MouseType"=mousetype,
                            "Mouse"=mouse,
                            "Sex" =MouseAttribute$Sex,
                            "Type"=type,
                            "State"=state,
                            "Day"=day), 300, 150)

                detCounter <- detCounter + 1

            } else if (type == "Prob") {

                # Filter out determenistic part and compute stats for the entire probabilistic part
                expProbDetPart <- filterSeqByTime(experiment,0,attribute$Probabilistic_Switch)
                ProbPartDetStats <- mouseExpStats(expProbDetPart, attribute,"Probabilistic Experiment, Deterministic Part")

                expProb <- filterSeqByTime(experiment,
                                            attribute$Probabilistic_Switch,
                                            attribute$Duration)

                ProbPartEntireStats <- mouseExpStats(expProb, attribute,"Probabilistic Experiment, Probabilistic Part")

                # Get summary stats for the probabilistic part split into parts n
                ProbPartStatsSplit <- splitExperiment(experiment,attribute,numberOfSplits,TRUE)

                ProbDataInfo <- c("MouseType"=mousetype,
                                    "Mouse"=mouse,
                                    "Sex" =MouseAttribute$Sex,
                                    "Type"=type,
                                    "State"=state,
                                    "Day"=day)

                ProbDetPartA <- c(ProbDataInfo,ProbPartDetStats)
                ProbProbPartA <- c(ProbDataInfo,ProbPartEntireStats)

                extractedFeaturesProbDetPart[[paste0(mousetype,"_",day,"_",state,dataType,type,runif(1))]] <- unname(ProbDetPartA)
                extractedFeaturesProbProbPart[[paste0(mousetype,"_",day,"_",state,dataType,type,runif(1))]] <- unname(ProbProbPartA)

                for(i in 1:numberOfSplits) {
                    tempStr <- paste0(mousetype,"_",day,"_",state,dataType,type,runif(1))
                    splitExperiments[[i]][[tempStr]] <- unname(c(ProbDataInfo,ProbPartStatsSplit[[i]]))
                }

                extractedFeaturesRollingProb[[index_string]] <- RollingWindowStats(experiment,attribute, ProbDataInfo, 300, 150)

                probDetCounter <- probDetCounter + 1
                

            }
            print(paste0("Finished processing ",mousetype," ",mouse," ",type," ",state," ",day," ",length(experiment[["Bandit"]])))

            
        }

        }

        }
            
        }

        }


    # Convert the lists to dataframes and add column names

    extractedFeaturesDetDF <- do.call(rbind.data.frame, extractedFeaturesDet)
    colnames(extractedFeaturesDetDF) <- names(ProbProbPartA)
    extractedFeaturesDetDF <- readr::type_convert(extractedFeaturesDetDF)

    extractedFeaturesProbDetPartDF <- do.call(rbind.data.frame, extractedFeaturesProbDetPart)
    colnames(extractedFeaturesProbDetPartDF) <- names(ProbProbPartA)
    extractedFeaturesProbDetPartDF <- readr::type_convert(extractedFeaturesProbDetPartDF)

    extractedFeaturesProbProbPartDF <- do.call(rbind.data.frame, extractedFeaturesProbProbPart)
    colnames(extractedFeaturesProbProbPartDF) <- names(ProbProbPartA)
    extractedFeaturesProbProbPartDF <- readr::type_convert(extractedFeaturesProbProbPartDF)

    # Create a list of the dataframes for the split experiments
    splitExperimentsDF <- list()
    for(i in 1:numberOfSplits) {
        tempp <- do.call(rbind.data.frame, splitExperiments[[i]]) 
        colnames(tempp) <- names(ProbProbPartA)
        splitExperimentsDF[[i]] <- tempp %>% readr::type_convert()
    }
    names(splitExperimentsDF) <- names(ProbPartStatsSplit)

    print("Finished processing all experiments.")
    print(paste0(detCounter," deterministic experiments and ",probDetCounter," probabilistic experiments were processed."))
    print(paste0(missExcelCounter," experiments were missing telemetry data."))
    #print(missingsTelemetryData)
    
    return(list("DetExperiments" = extractedFeaturesDetDF,
                "ProbExperimentsDetPart" = extractedFeaturesProbDetPartDF,
                "ProbExperimentsProbPart" = extractedFeaturesProbProbPartDF,
                "SplitProbExperiments" = splitExperimentsDF,
                "MissingTelemetryData" = missingsTelemetryData,
                "RollingProbExperiments" = extractedFeaturesRollingProb,
                "RollingDetExperiments" = extractedFeaturesRollingDet))
                

    
}


RollingWindowStats <- function(experiment, attribute, attdata ,windowLength, overlap,Prob) {

    #if(attdata["Type"] == "Prob") {start <- 300} else {start <- 0}
    start <- 0
    # Generate the indices for the windows
    windowStarts <- seq(start,attribute$Duration-windowLength,windowLength-overlap)
    windowEnds <- windowStarts + windowLength

    # Create a list for storing the results
    res <- list()
    gotit <- FALSE

    for( i in 1:length(windowStarts) ) {s
        expTemp <- filterSeqByTime(experiment,windowStarts[i],windowEnds[i])
        temp <- c(attdata,mouseExpStats(expTemp, attribute,
                paste0("t_",windowStarts[i],"_",windowEnds[i])))
        if(!gotit & sum(is.na(temp)) < 2) {colnames <- names(temp)}
        res[[paste0("t_",windowStarts[i],"_",windowEnds[i])]] <- c(paste0("t_",windowStarts[i],"_",windowEnds[i]),unname(temp))
    }

    res <- do.call(rbind.data.frame, res) %>% readr::type_convert()
    colnames(res) <- c("Range",colnames)

    return(res)

}


gatherExperimentsRolling <- function(ListOfDataframes,variablesOfInterest) {
        # Hacky way to get the correct column names (hacks on hacks on hacks)
        nameVec <- (NA)
        i <- 1
        while(any(is.na(nameVec))) {
                nameVec <- names(ListOfDataframes[[i]])
                i <- i + 1
        }

        # Process the data frames, rename columns and interpolate
        for(i in 1:length(ListOfDataframes)) {
                # Rename the columns
                colnames(ListOfDataframes[[i]]) <- nameVec

                # Interpolate the data missing data by taking the previous non na value
                ListOfDataframes[[i]] <- ListOfDataframes[[i]] %>%
                        mutate_at(vars(variablesOfInterest),~zoo::na.approx(.,na.rm = FALSE))

                # Remove all rows that have any na values for the variables of interest
                ListOfDataframes[[i]] <- ListOfDataframes[[i]] %>%
                        filter_at(vars(variablesOfInterest),all_vars(!is.na(.)))


         }

       temp <- do.call(rbind.data.frame, ListOfDataframes)

        colnames(temp) <- nameVec

        temp <- temp %>% 
                mutate(TimeRange = str_remove(Range,"t_")) %>%
                tidyr::separate(TimeRange,c("Start","End"),sep = "_",remove = TRUE) %>%
                mutate(MidTimePoint = (as.numeric(Start) + as.numeric(End))/2) 

       return(select(temp,all_of(c("Range","MidTimePoint","MouseType","Mouse","Sex","State","Day","NTrials",variablesOfInterest))))

}
