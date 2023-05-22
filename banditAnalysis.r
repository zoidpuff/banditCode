
# THIS FILE CONTAINS A NUMBER OF FUNCTIONS THAT ARE USED IN THE BANDIT ANALYSIS
library(tidyverse)



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
          choiceCount <- table(factor(Choice) , levels = c("A","B","C"))
          # Depending on what the current bandit is, sample the choice with the lowest count
          if(Bandit[i] == "AB"){
            Choice[i] <- names(which.min(choiceCount[1:2]))
          } else if(Bandit == "AC"){
            Choice[i] <- names(which.min(choiceCount[c(1,3)]))
          } else if(Bandit == "BC"){
            Choice[i] <- names(which.min(choiceCount[2:3]))
          }
         }

        # Sample the time according to the bandit/choice combo and the banditTimeVec with std dev of 2
        Time[i] <- Time[i-1] + rnorm(1,timePerBandit[paste0(Bandit[i],">",Choice[i])],2)

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
addPrefix <- function(x,prefix){
    names(x) <- paste0(prefix,names(x))
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
    banditSplit <- str_split(bandit,">")
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
        return(paste0(isAt,">",pumbProb_A))
    }else if (pumpChoice == "B"){
        return(paste0(isAt,">",pumpProb_B))
    }else if (pumpChoice == "C"){
        return(paste0(isAt,">",pumpProb_C))
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
    banditChoices <- .banditChoice(experiment[["Choice"]],experiment[["Bandit"]],A_Prob,B_Prob,C_Prob)

    # If first string starts with NA, remove it
    if(str_starts(banditChoices[1],"NA")){
        banditChoices <- banditChoices[-1]
        } 

    if(!returnStats){
        return(banditChoices)
    }
    banditChoiceCount <- table(banditChoices)

    # Reorder the vectors and make sure they have the all all the bandits
    banditChoiceCount <- banditChoiceCount[c("0.25>0.5","0.25>1", "0.5>0.25","0.5>1","1>0.25","1>0.5")] 

    banditChoiceCount[is.na(banditChoiceCount)] <- 0

    # Compute the probability of choosing each bandit
    banditChoiceProb <- as_tibble(banditChoiceCount) %>%
                            separate("banditChoices",c("from","to"),sep = ">",remove = FALSE) %>%
                            group_by(from) %>%
                            mutate(Prob = n/sum(n))

    # Convert tibble to a named vector
    banditChoiceProb <- setNames(banditChoiceProb$Prob,banditChoiceProb$banditChoices)

    banditChoiceProb <- banditChoiceProb[c("0.25>0.5","0.25>1", "0.5>0.25","0.5>1","1>0.25","1>0.5")] %>%
                            addPrefix("banditChoiceProb_")

    banditChoiceProb[is.na(banditChoiceProb)] <- 0

    banditChoiceCount <- banditChoiceCount %>% addPrefix("banditChoiceCount_")

    # Compute the  of the bandit choices
        # First split the bandit choice three vectors based in bandit
        which(startsWith(banditChoices,"0.25")) -> inds_025
        which(startsWith(banditChoices,"0.5")) -> inds_05
        which(startsWith(banditChoices,"1")) -> inds_1

        # Compute the cramerV lagged for each bandit
        cramer_025 <- cramersVautoCorr(banditChoices[inds_025],1) %>% addPrefix("cramer_025")
        cramer_05 <- cramersVautoCorr(banditChoices[inds_05],1) %>% addPrefix("cramer_05")
        cramer_1 <- cramersVautoCorr(banditChoices[inds_1],1) %>% addPrefix("cramer_1")
        cramers <- c(cramer_025,cramer_05,cramer_1)

    # Compute stats around making the optimal choice
        # Check if the mouse made the optimal choice for each bandit in the sequence
        banditOptimalChoiceSeq <- optimalChoice(banditChoices) %>% as.logical()
        banditChoicesOptimalProb <- sum(banditOptimalChoiceSeq)/length(banditOptimalChoiceSeq) %>% 
                                    addPrefix("banditChoicesOptimalProb_")
        
        # Compute the cramerV for optimal choice
        cramerOptimalChoice <- cramersVautoCorr(as.character(banditOptimalChoiceSeq),1) %>% 
                                addPrefix("cramerVOptimalChoice")

    return(c(banditChoiceCount,
            banditChoiceProb,
            cramers,
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
                    addPrefix("ProbOfChoosing_")

    pumpCountRel[is.na(pumpCountRel)] <- 0
    pumpProbAbs[is.na(pumpProbAbs)] <- 0


    pumpProbAbs <- pumpProbAbs %>% addPrefix("ProbOfChoosing_")
    
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

    # print(paste0("Computing Circle Stats with ",length(directionChoices)," choices"))
    # Check if the sequence is atleast 3 long else just return na
    if(length(directionChoices) > 2) {
        countCircles <- 0
        for (i in 1:(length(directionChoices)-2)){
            if (directionChoices[i] == directionChoices[i+1] & directionChoices[i+1] == directionChoices[i+2]){
                countCircles <- countCircles + 1
            }
        }
        ProbOfCircle <- countCircles/(length(directionChoices)-2) %>% addPrefix("ProbOfCircle")


    
        #Count number of times the mouse changes direction
        countDirectionChanges <- 0
        for(i in 1:(length(directionChoices)-1)){
                if (directionChoices[i] != directionChoices[i+1]){
                    countDirectionChanges <- countDirectionChanges + 1
                }
            }
    
        probDirectionChanges <- countDirectionChanges/length(directionChoices) %>% addPrefix("probDirectionChanges")
        countCircles <- addPrefix(countCircles,"circleCount")
        countDirectionChanges <- addPrefix(countDirectionChanges,"countDirectionChanges")
    } else {
        countCircles <- c(NA) %>% addPrefix("countCircles")
        ProbOfCircle <- c(NA) %>% addPrefix("ProbOfCircle")
        probDirectionChanges <- c(NA) %>% addPrefix("probDirectionChanges")
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
                                separate("banditChoices",c("from","to"),sep = ">",remove = FALSE) %>%
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

    rewardProb <- sum(reward)/length(reward) %>% addPrefix("RewardProb")


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
    ObsRewardProbPerBandit <- ObsRewardProbPerBandit[c("0.25>0.5","0.25>1", "0.5>0.25","0.5>1","1>0.25","1>0.5")] %>% 
                                addPrefix("ObsRewardProbPerBandit")

    # Compute the observed reward per station
    ObsProbRewardPerStation <- data.frame("choice" = banditChoices,
                                    "reward" = reward) %>%
                        separate("choice",c("from","to"),sep = ">",remove = FALSE) %>%
                        group_by(to) %>%
                        summarise(rewardPerStation = mean(reward))


    # Convert tibble to a named vector and order
    ObsProbRewardPerStation <- setNames(ObsProbRewardPerStation$rewardPerStation,ObsProbRewardPerStation$to)
    ObsProbRewardPerStation <- ObsProbRewardPerStation[c("0.25","0.5","1")] %>% addPrefix("ObsProbRewardPerStation_")

    # Compute the probability of getting a reward after making the optimal choice
    optimalOrNot <- optimalChoice(banditChoices) %>% as.logical()
    rewardProbAfterOptimalChoic <- sum(reward[optimalOrNot])/sum(optimalOrNot) %>% addPrefix("rewardAfterOptimal")

    return(c(totalReward,
            rewardProb,
            ObsRewardProbPerBandit,
            ObsProbRewardPerStation,
            rewardProbAfterOptimalChoic))

}


# Function for determining some stats given a experiment object and its attribute
mouseExpStats <- function(experiment,attribute,inputDataType){

    #print(paste0("Computing experiment stats on :" ,inputDataType,". Number of Trials: ",length(experiment[["Choice"]])))
    # Get number of trials
    numberOfTrails <- length(experiment[["Choice"]]) %>% addPrefix("numberOfTrials")

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

    res <- c(numberOfTrails,banditChoices,biasStats,circleStats,mouseSpeeds,rewardAttainment)
    #print(res) debug
    return(res)
}
computeRollingPumpPreference <- function(experiment,attribute,inputDataType){
  # Get the estimate 
    banditChoices <- banditChoice(experiment, attribute, FALSE)

}



tabulateExperiments <- function(data,datafile,numberOfSplits) {
    # Create lists for storing the results
    extractedFeaturesDet <- list()
    extractedFeaturesProbDetPart <- list()
    extractedFeaturesProbProbPart <- list()
    # Create a nested list for the split probabilistic experiments
    splitExperiments <- list()
    for(i in 1:numberOfSplits) {
        splitExperiments[[i]] <- list()
    }
    detCounter <- 0
    probDetCounter <- 0
    # Loop over the data
    for (mousetype in names(data)) { # WT

    for (mouse in names(data[[mousetype]])) { # Mouse

            sex <- h5readAttributes(datafile, paste0("/",mousetype,"/",mouse,"/" ))$Sex

            for (type in names(data[[mousetype]][[mouse]])) { # Det vs prob

                for (state in names(data[[mousetype]][[mouse]][[type]])) { # Fasted vs fed

                    for (day in names(data[[mousetype]][[mouse]][[type]][[state]])) { # Day
                            
                            dataType = "Action"
                                # Get the experiment
                                experiment <- data[[mousetype]][[mouse]][[type]][[state]][[day]][[dataType]]
                                attribute <- h5readAttributes(datafile, paste0("/",mousetype,"/",mouse,"/",type,"/",state,"/",day,"/",dataType,"/"))
                                MouseAttribute <- h5readAttributes(datafile, paste0("/",mousetype,"/",mouse,"/"))
                                # Check if experment is Det or Prob
                                if(type == "Det") {
                                    print("------------------")
                                    print(paste0("Started processing ",mousetype," ",mouse," ",type," ",state," ",day," ",length(experiment[["Bandit"]])))
                                    # Get the bandit choices, direction choices, and pump choices for all trials

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
                                    detCounter <- detCounter + 1
                                    #print(str(experiment))
                                    #print(str(attribute))
                                    #stop("debug")



                                } else if (type == "Prob") {
                                    print("------------------")
                                    print(paste0("Started processing ",mousetype," ",mouse," ",type," ",state," ",day," ",length(experiment[["Bandit"]])))
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
                                    probDetCounter <- probDetCounter + 1

                                    ### Debugging
                                    #print(ProbProbPartA)
                                    #print(attribute)
                                    #View(data.frame(experiment[["Bandit"]],
                                    #                experiment[["Choice"]],
                                    #                experiment[["Reward"]],
                                    #                experiment[["Time"]],
                                    #                c(NA,banditChoice(experiment, attribute, FALSE)),
                                    #                c(NA,biasStats(experiment, attribute,TRUE))))
                                    #View(data.frame(expProb[["Bandit"]],
                                    #                expProb[["Choice"]],
                                    #                expProb[["Reward"]],
                                    #                expProb[["Time"]],
                                    #                c(banditChoice(expProb, attribute, FALSE)),
                                    #                c(biasStats(expProb, attribute,TRUE))))
                                    


                                }
                                print(paste0("Finished processing ",mousetype," ",mouse," ",type," ",state," ",day," ",length(experiment[["Bandit"]])))
                            #dataType = "Lick"
                            #  exp <- data[[mousetype]][[mouse]][[type]][[state]][[day]][[dataType]]
                            # Code_Lick <- paste0(exp[["Code_Lick"]],collapse = ",")
                            # Pump_Lick <- paste0(exp[["Pump_Lick"]],collapse = ",")
                            # Time_Lick <- paste0(exp[["Time_Lick"]],collapse = ",")
                        
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
    
    return(list("DetExperiments" = extractedFeaturesDetDF,
                "ProbExperimentsDetPart" = extractedFeaturesProbDetPartDF,
                "ProbExperimentsProbPart" = extractedFeaturesProbProbPartDF,
                "SplitProbExperiments" = splitExperimentsDF
                ))

    
}


panel.cor <- function(x, y, digits=2, prefix="", use="pairwise.complete.obs", method="pearson", cex.cor, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- cor(x, y, use=use, method=method) # MG: remove abs here
        txt <- format(c(r, 0.12345678910), digits=digits)[1]
        txt <- paste(prefix, txt, sep="")
        if(missing(cex.cor)) cex <- 0.6/strwidth(txt)

        test <- cor.test(as.numeric(x),as.numeric(y), method=method)
        # borrowed from printCoefmat
        Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                    symbols = c("***", "**", "*", ".", " "))
        # MG: add abs here and also include a 30% buffer for small numbers
        text(0.5, 0.5, txt, cex = cex * (abs(r) + .2) / 1.3)
        text(.8, .8, Signif, cex=cex, col=2)
    }

    hist.panel = function (x, ...=NULL ) {
        par(new = TRUE)
        hist(x,
             col = "light gray",
             probability = TRUE,
             axes = FALSE,
             main = "",
             breaks = "FD")
        lines(density(x, na.rm=TRUE),
              col = "red",
              lwd = 1)
        rug(x)
    }


# Function to avarage the data over replicated experiments
averageReplicates <- function(data) {
    data %>% 
        group_by(State,Mouse,Sex) %>%
        summarise_if(is.numeric,mean,na.rm=TRUE) %>%
        ungroup() %>%
        as.data.frame() %>%
        return()
}


# Create a PCA plotting function
plotPCAMouse <- function(data,columns,name,umapNeigbours,umapMinDist,CreateSDPlots) {
    # Remove any rows with NA values
    dataTrim <- data[complete.cases(data[,columns]),]
    # Print how many rows were removed
    print(paste0("Removed ",nrow(data)-nrow(data[complete.cases(data[,columns]),])," rows with NA values."))
    # Get the labels for the data
    dataLabs <- dataTrim[,1:6]
    # Select the columns to use for the PCA
    dataPCA <- dataTrim[,columns]
    # Remove columns with only one unique value
    dataPCA <- dataPCA[,sapply(dataPCA, function(x) length(unique(x)) > 1)]
    # Plot the Correlation matrix to see how the features correlate
    group <- as.factor(dataLabs$Mouse) %>% as.numeric()
    cols <- RColorBrewer::brewer.pal(n=length(unique(group)), name="Set3")[group]

    # Create folder if it doesnt exist
    dir.create((paste0("./plots/",name)))
    
    #dataPCA <- SecondHalfCollapsedFastedOnly[,TesmerFeatures]
    pdf(paste0("./plots/",name,"/_PairCorPlot.pdf"),width = ncol(dataPCA)/1.5,height = ncol(dataPCA)/1.5)
    #pdf(paste0("./plots/secondHalfFastedOnly.pdf"),width = ncol(dataPCA)/1.5,height = ncol(dataPCA)/1.5)
    graphics::pairs(dataPCA,
                    gap=0, 
                    lower.panel=graphics::panel.smooth, 
                    upper.panel=panel.cor, 
                    diag.panel=hist.panel,
                    labels=gsub("(.{12})", "\\1\n",colnames(dataPCA)),
                    pch=20,
                    col=cols,
                    font.labels=2,
                    cex.labels=0.7)
    dev.off()

    if(CreateSDPlots){
    # Plot a heatmap of SD for each var and each mouse stratified by state
    varianceMat <-  dataTrim %>%
                    group_by(Sex,Mouse,State) %>%
                    summarise_if(is.numeric,sd,na.rm = TRUE)
    
    varianceMatGrps <- varianceMat %>% ungroup() %>% select(Sex,Mouse,State)

    varianceMat <- varianceMat %>% 
                            ungroup() %>%
                            mutate_if(is.numeric,scales::rescale) %>%
                            as.data.frame()

    rownames(varianceMat) <- paste0(varianceMatGrps$Mouse,":",varianceMatGrps$Sex,"_",varianceMatGrps$State)

    annotationsRow <- data.frame(Mouse = as.factor(varianceMatGrps$Mouse),
                        State = as.factor(varianceMatGrps$State),
                        Sex = as.factor(varianceMatGrps$Sex),
                        row.names = rownames(varianceMat))

    pdf(paste0("./plots/",name,"/_SDPerMouseAndStateHeatMapPlot.pdf"),width = 12,height = 12)
    pheatmap::pheatmap(varianceMat[,columns],
                        annotation_row = annotationsRow,
                        show_rownames = FALSE,
                        main = "Scaled Standard Deviation Heatmap")
    dev.off()


    # Plot a heatmap of SD for each var and each mouse 
    varianceMat <-  dataTrim %>%
                    group_by(Sex,Mouse) %>%
                    summarise_if(is.numeric,sd,na.rm = TRUE)
    
    varianceMatGrps <- varianceMat %>% ungroup() %>% select(Sex,Mouse)

    varianceMat <- varianceMat %>% 
                            ungroup() %>%
                            mutate_if(is.numeric,scales::rescale) %>%
                            as.data.frame()

    rownames(varianceMat) <- paste0(varianceMatGrps$Mouse,":",varianceMatGrps$Sex)

    annotationsRow <- data.frame(Mouse = as.factor(varianceMatGrps$Mouse),
                        Sex = as.factor(varianceMatGrps$Sex),
                        row.names = rownames(varianceMat))

    pdf(paste0("./plots/",name,"/_SDperMouseHeatMapPlot.pdf"),width = 12,height = 12)
    pheatmap::pheatmap(varianceMat[,columns],
                        annotation_row = annotationsRow,
                        show_rownames = FALSE,
                        main = "Scaled Standard Deviation Heatmap")
    dev.off()
    }

    # Create a heatmap of features correlation
    pdf(paste0("./plots/",name,"/_FeatureCorrPlot.pdf"),width = 12,height = 12)
    pheatmap::pheatmap(cor(dataPCA,method = "pearson"),
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        main = "Features Correlation Heatmap")
    dev.off()

    # Plot a heatmap of samples correlation
        # Create annotation dataframe
        rownames(dataPCA) <- as.character(1:nrow(dataPCA))
        ann <- data.frame(Mouse = as.factor(dataLabs$Mouse),
                        State = as.factor(dataLabs$State),
                        Sex = as.factor(dataLabs$Sex),
                        row.names = rownames(dataPCA))
        

        # Create the heatmap with pheatmap
        pdf(paste0("./plots/",name,"/_ExperimentCorrPlot.pdf"),width = 12,height = 12)
        pheatmap::pheatmap(cor(t(dataPCA),method = "spearman"),
                                        annotation_col = ann,
                                        annotation_row = ann,
                                        show_rownames = FALSE,
                                        show_colnames = FALSE,
                                        main = "Experiment Correlation Heatmap")	
        dev.off()
    # Scale the data
    dataPCA <- scale(dataPCA)
    # Compute the PCA
    pca <- prcomp(dataPCA,center = FALSE) 
    var_explained = pca$sdev^2 / sum(pca$sdev^2)

    # Plot the scree plot with ggplot2
    screePlot <- ggplot(data.frame(var_explained = var_explained, PC = 1:length(var_explained)), 
                    aes(x = PC, y = var_explained)) +
                    geom_bar(stat = "identity", fill = "steelblue") +
                    geom_line(stat = "identity", color = "red", size = 1) +
                    scale_x_continuous(breaks = seq(1, length(var_explained), 1)) +
                    labs(x = "Principal Component",
                        y = "Variance Explained (%)",
                        title = "Scree Plot") +
                    theme_bw()	


    # Plot the PCA with ggplot2
    pcaPlot <-ggplot(data = dataLabs,
                            aes(x = pca$x[,1],
                                y = pca$x[,2],
                                color = as.factor(Mouse),
                                shape=as.factor(State))) +
                    geom_point(size = 4) +
                    labs(x = paste0( "PC1 (Variance Explained:", signif(var_explained[1]*100,4),"%)"),
                        y = paste0( "PC2 (Variance Explained:", signif(var_explained[2]*100,4),"%)"), color = "Mouse",shape = "State") +
                    theme_bw()


    # Plot a pairplot of the PCA
    pdf(paste0("./plots/",name,"/_PairPCAPlot.pdf"),width = ncol(dataPCA)/1.5,height = ncol(dataPCA)/1.5)
    graphics::pairs(pca$x,
                    gap=0, 
                    lower.panel=graphics::panel.smooth, 
                    upper.panel=panel.cor, 
                    diag.panel=hist.panel,
                    labels=gsub("(.{12})", "\\1\n",colnames(pca$x)),
                    pch=20,
                    col=cols,
                    font.labels=2,
                    cex.labels=0.7)
    dev.off()



    # Plot loadings PCA as a heatmap
      # Create color gradient function
      scaled_loadings <- sweep(pca$rotation, 2, var_explained, "*")
      colfunc <- colorRampPalette(c("blue", "white", "red"))

      # Get the maximum absolute value in the loadings matrix
      max_val <- max(abs(scaled_loadings))

      # Generate color breaks
      color_breaks <- seq(-max_val, max_val, length.out = 100)

    
      pdf(paste0("./plots/",name,"/_LoadingsHeatmapPlot.pdf"),width = 12,height = 12)
      pheatmap::pheatmap(scaled_loadings,
                          show_rownames = TRUE,
                          show_colnames = TRUE,
                          cluster_rows = FALSE,
                          cluster_cols = FALSE,
                          color = colfunc(100), 
                          breaks = color_breaks,
                          main = "Scaled Loadings Heatmap")
      dev.off()
                      

    # Plot the PCA with ggbiplot

    pcaBiPlotMice <-ggbiplot::ggbiplot(pca,ellipse=TRUE,
                        groups=as.factor(dataLabs$Mouse),
                        ellipse.prob = 0.5) + 
                        theme_bw() + 
                        ggtitle("BiPlot of PCA Grouped on Mice") +
                        scale_color_discrete(name = "Mouse")

    pcaBiPlotMice2 <-ggbiplot::ggbiplot(pca,ellipse=TRUE,choices = c(3,4),	
                        groups=as.factor(dataLabs$Mouse),
                        ellipse.prob = 0.5) + 
                        theme_bw() + 
                        ggtitle("BiPlot of PCA Grouped on Mice") +
                        scale_color_discrete(name = "Mouse")

    pcaBiPlotState <-ggbiplot::ggbiplot(pca,ellipse=TRUE,
                        groups=as.factor(dataLabs$State),
                        ellipse.prob = 0.5) + 
                        theme_bw() + 
                        ggtitle("BiPlot of PCA Grouped on State") +
                        scale_color_discrete(name = "Mouse")

    pcaBiPlotState2 <-ggbiplot::ggbiplot(pca,ellipse=TRUE,choices = c(3,4),
                        groups=as.factor(dataLabs$State),
                        ellipse.prob = 0.5) + 
                        theme_bw() + 
                        ggtitle("BiPlot of PCA Grouped on State") +
                        scale_color_discrete(name = "Mouse")

    # Compute a umap of the data

    umap <- umap::umap(dataPCA,
                        method = "umap-learn",
                        n_neighbors = umapNeigbours,
                        min_dist = umapMinDist,
                        metric = "correlation")

    # Plot the umap with ggplot2

    umapPlot <- ggplot(data = dataLabs,
                        aes(x = umap$layout[,1],
                            y = umap$layout[,2],
                            color = as.factor(Mouse),
                            shape=as.factor(State))) +
                    geom_point(size = 4) +
                    labs(x = "UMAP1",title = "UMAP",
                        y = "UMAP2", color = "Mouse",shape = "State") +
                    theme_bw()

    # Save the ggplots to pdf
    ggplot2::ggsave(paste0("./plots/",name,"/_ScreePlot.pdf"),device = "pdf",screePlot)
    ggplot2::ggsave(paste0("./plots/",name,"/_PCAPlot.pdf"),device = "pdf",pcaPlot)
    ggplot2::ggsave(paste0("./plots/",name,"/_PCABiPlotGroupedMicePC12.pdf"),device = "pdf",pcaBiPlotMice)
    ggplot2::ggsave(paste0("./plots/",name,"/_PCABiPlotGroupedMicePC34.pdf"),device = "pdf",pcaBiPlotMice2)
    ggplot2::ggsave(paste0("./plots/",name,"/_PCABiPlotGroupedStatePC12.pdf"),device = "pdf",pcaBiPlotState)
    ggplot2::ggsave(paste0("./plots/",name,"/_PCABiPlotGroupedStatePC34.pdf"),device = "pdf",pcaBiPlotState2)
    ggplot2::ggsave(paste0("./plots/",name,"/_UMAPPlot.pdf"),device = "pdf",umapPlot)


    return("Done")
}
     