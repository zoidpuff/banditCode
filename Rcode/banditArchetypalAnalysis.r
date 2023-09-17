# This file contains functions that used to perfom archtypal analysis on bandit mouse data



#############################################################################################
# # # # # # # # # # # # # # # # Manual Archetypal analysis # # # # # # # # # # # # # # # # # 
#############################################################################################

# Some helper functions 

## this is just cuz I dont wanna figure out how to feed the ^2 operator as a function argument, this is dumb 
squareVec <- function(x) {
    return(x^2)
}

# Find the value in the vector that is closest to some value
closeVal <- function(val,vec) {
    # remove any outliers in vec (kinda cursed)
    top <- quantile(vec,0.98)
    bottom <- quantile(vec,0.02)	
    vec <- vec[vec < top & vec > bottom]
    return(vec[which.min(abs(vec - val))])
}

archetypesNNLS <- function(archeMat,data,bigM = 200) {
        
        ##### COMPUTE ALPHAS WITH NNLS #####

    # Initialize the matrix of alphas
    alphas <- matrix(0,nrow = nrow(data), ncol = nrow(archeMat))

    # Scale the data matrix and store the scaling factors to use on archetypes matrix
    scaledData <- scale(data)

    # Add constraint dummy variable to the data matrix
    scaledDataDummy <- cbind(scaledData,bigM)

    # Scale the archetypes matrix with the same scaling factors as the data matrix
    scaledArchMat <- scale(archeMat,
                        center = attr(scaledData,"scaled:center"),
                        scale = attr(scaledData,"scaled:scale"))
    
    # Add big M constraint dummy variable to the archetypes matrix
    scaledArchMatDummy <- cbind(scaledArchMat,bigM)

    # For each data point compute the alphas
    for(i in 1:nrow(scaledData)) {
        alphas[i,] <- coef(nnls::nnls(t(scaledArchMatDummy),scaledDataDummy[i,]))
    }

    # Compute the RSS
    Approx <-  t(archeMat) %*% t(alphas)
    RSS <- sum((Approx - t(data))^2)

    # Convert the alphas to a dataframe
    alphas <- alphas %>% as.data.frame()

    colnames(alphas) <- rownames(archeMat)

    return(list(alphas,RSS))
}


##############
# OLD method #
##############

# Find  "optimal" alphas with gridsearch with some given archtypes matrix and data matrix. Specify a certain resolution and loss to use. 
archtypesGridSearch <- function(archeMat, resolution = 0.001, data, mouseData, lossFunction,scaling,returnHeatmap = FALSE,takeSqrt = TRUE,dataSetName = "NA") {

     ### OLD GRIDSEARCH CODE ####

    # Transpose the archetypes for lin alg reasons
    tarcheMat <- t(archeMat)

    # Initialize the matrix of alphas
    alphas <- matrix(0,nrow = ncol(tarcheMat), ncol = nrow(data))

    # Scale the data matrix and store the scaling factors to use on archetypes matrix
    ## This is so that features are of the same scale and contribute equally to the loss function
    if(scaling) {
        data <- scale(data, center = FALSE, scale = TRUE) # Dont think centering is necessary

        # Scale the archetypes matrix with the same scaling factors as the data matrix
        tarcheMat <- tarcheMat / attr(data, "scaled:scale")  
    }
    # Create a list of all possible combinations of the archetypes given resolution and store in matrix    
    matrixList <- list()
    for(i in seq(0, 1, by = resolution)) {
        if(i == 1) {j_seq <- c(0)} else {
            j_seq <- round(seq(0, 1-i, by = resolution), 7)
        }
        k_seq <- round((1 - i - j_seq), 7)
        matrixList[[as.character(i)]] <- matrix(c(rep(i,length(j_seq)), j_seq, k_seq), ncol = 3)
    }
    allcombs <- do.call(rbind, matrixList)


    if(nrow(archeMat)==4) {
        resForQuad <- 0.005
        matrixList <- list()
        for(i in seq(0, 1, by = resForQuad)) {
            if(i == 1) {j_seq <- c(0)} else {
                j_seq <- round(seq(0, 1-i, by = resForQuad), 7)
            }
            for(j in j_seq) {
                if(j == 1 | i == 1) {k_seq <- c(0)} else {
                    k_seq <- round(seq(0, max(c((1-i-j),0)), by = resForQuad), 7)
                }
                    l_seq <- round((1 - i - j - k_seq), 7)
                        matrixList[[paste(i,".",j)]] <- matrix(c(rep(i,length(k_seq)), rep(j,length(k_seq)), k_seq, l_seq), ncol = 4)
            }
        }
        allcombs2 <- do.call(rbind, matrixList)
        allcombs <- allcombs2
    }



    # Initialize the list of heatmaps
    heatPlots <- list()


    errorTotal <- 0
    for(d in 1:nrow(data) ){
        # print(d)
        dataPoint <-  data[d,]
        Aprox <-  tarcheMat %*% t(allcombs)
        ers <- colSums(lossFunction(sweep(Aprox, 1, t(dataPoint), "-")))
        if(takeSqrt) {
            ers <- sqrt(ers)
        }

        minError <- min(ers)
        errorTotal <- errorTotal + minError
        currBest <- allcombs[which.min(ers),]
        alphas[,d] <- currBest
        

        if(returnHeatmap) {
            heatmapPrim <- cbind(allcombs[,1:2],ers)
            colnames(heatmapPrim) <- c("Circle","GO","Error")
            heatmapPrim <- as.data.frame(heatmapPrim)
            #   heatmapPrim <- heatmapPrim %>% 
            #     pivot_wider(names_from = Circle, values_from = Error) %>% 
            #     column_to_rownames("GO")

            a <- ggplot(heatmapPrim,aes(x = Circle, y = GO, fill = Error)) +
                                geom_tile() + 
                                scale_fill_gradient(low = "#075AFF", high = "#FF0000") +
                                scale_x_continuous(breaks = NULL) +
                                scale_y_continuous(breaks = NULL) + 
                                theme_classic() +
                                coord_fixed() +
                                ggtitle(paste0("RMSE landscape for mouse: ",mouseData$Mouse[d], " ", mouseData$State[d]),
                                    subtitle = paste0("Minimum RMSE: ",round(minError,4))) +
                                annotate("point", x = currBest[1], y = currBest[2], size = 7, color = "#000000",shape = "x") + 
                                ylab("GO-Uncertainty") 
            ggplot2::ggsave(paste0("./plots/",dataSetName,"/",dataSetName,"_Mouse",mouseData$Mouse[d],"_State",mouseData$State[d],"_RMSE.jpeg"),device = "jpeg",a)
        }

    }

    alphas <- alphas %>% t() %>% as.data.frame()
    colnames(alphas) <- rownames(archeMat)


    return(list(alphas,errorTotal/nrow(data)))
    
}



#############################################################################################
# # # # # # # # # # # # # # # # Plotting Functions  # # # # # # # # # # # # # # # # # # # # #
#############################################################################################

# Plot the ternary plot with arrows
plotTernaryWithArrows <- function(alphasThree, Grouplabels, contrast,zoom,colVar,title=""){

    # Filter for the specific contrast

    filtInds <- Grouplabels$State %in% contrast

    varNames <- colnames(alphasThree)

    p <- ggtern::ggtern(data = alphasThree[filtInds,],
                aes_string(x = varNames[1], y = varNames[2], z = varNames[3],
                            color = paste0("as.character(Grouplabels$",colVar,")[filtInds]"),
                            group = "Grouplabels$Mouse[filtInds]")) + 
                geom_point(size=3,shape= Grouplabels$Sex[filtInds]) + 
                geom_line(arrow = arrow(length=unit(0.15,"cm"),
                            ends="both", type = "open"),
                            color = "#000000",
                            alpha=0.15) +
                ggtern::theme_rgbw() +
                guides(color = guide_legend(title = colVar)) +
                labs(x = varNames[1], y = varNames[2], z = varNames[3],) + 
                ggtern::theme_zoom_center(x = zoom) +
                ggtitle("",subtitle = title) +
                ggtern::theme_hidegrid()
    return(p)
  }

# Compute and plot deltas given an experimentla grouping and a alpha matrix 

create_diff_plots <- function(df, archetype_list,title,pairs_to_calculate_diff) {

  diff_df <- df %>%
    pivot_wider(names_from = State, values_from = archetype_list)

  
  for (archetype in archetype_list) {
    for (pair in pairs_to_calculate_diff) {
      diff_df <- diff_df %>%
        mutate(!!paste0("diff_", archetype, "_", paste(pair, collapse="Minus")) := 
                 .data[[paste0(archetype, "_", pair[1])]] - .data[[paste0(archetype, "_", pair[2])]])
    }
  }

  diff_df <- diff_df %>%
    select(starts_with("diff_") | starts_with("Se") | starts_with("Mous")) %>%
    pivot_longer(cols = starts_with("diff_"), names_to = "State", values_to = "Delta") %>%
    mutate(State = str_remove(State, "diff_")) %>%
    separate(State, into = c("Var", "Group"), sep = "_")

    levelsVec <- c()
    for(i in 1:length(pairs_to_calculate_diff)) {
        temp <- paste0(pairs_to_calculate_diff[[i]][2]," -> ",pairs_to_calculate_diff[[i]][1])
        diff_df$Group <- str_replace(diff_df$Group,
                        paste0(pairs_to_calculate_diff[[i]][1],"Minus",pairs_to_calculate_diff[[i]][2]), 
                        temp)
        levelsVec <- c(levelsVec,temp)	
    }

    diff_df$Group <- factor(diff_df$Group, levels = levelsVec)

    diff_df$Var <- factor(diff_df$Var, levels = archetype_list)

   # Filter rows with NA values
    diff_df <- diff_df[complete.cases(diff_df),]
    # Plot the delta of arachetypes across the groups with a ggplot box plot facted by group fill by state
    # Add a new variable to shift the x-values
    diff_df$X_pos <- ifelse(diff_df$Sex == "m", as.numeric(as.factor(diff_df$Var)) + 0.2, 
                            ifelse(diff_df$Sex == "f", as.numeric(as.factor(diff_df$Var)) - 0.2, NA))
    


    p <- ggplot(diff_df, aes(x = X_pos, y = Delta, fill = Var)) +
            geom_boxplot(aes(x = as.numeric(as.factor(Var))), outlier.alpha = 0,alpha=0.75) +
            geom_jitter(aes(shape = Sex), alpha = 0.4,size=3) +
            facet_wrap(~Group,ncol=3) +
            theme_classic() +
            labs(x = "Archetype", y = "Delta (Value at State 2 - Value at State 1)", fill = "State") +  
            geom_hline(yintercept = 0, linetype = "dotted", color = "red") + 
            ggtitle(title) + 
            scale_shape_manual(values=c("f","m")) +
            scale_fill_manual(values=c("darkblue", "darkred","darkgreen","orange"))+
            guides(shape = "none") +
            scale_x_discrete(breaks = c("","",""))

    # Create a second plot which shows the relationship between a starting location in terms of some archetypal value vs where it ends up after a change in state

    # Function to get the value of a variable from a row (lol)
    get_value <- function(row, var_vec) {
        var <- row[["Var"]]
        if (var %in% var_vec) {
            return(row[[var]])
        } else {
            return(NA_real_)
        }
    }
    
    diff_df$Group <- as.character(diff_df$Group)

    # Create a dataframe that shows that has columns for the value of each archetype for each mouse in the base state and the delta when moving to the other state
    diff_df %>% separate(Group, into = c("From", "To"), sep = " -> ",remove = FALSE) %>%
                left_join(df,by= c("Mouse" = "Mouse", "From" = "State")) %>%
                  rowwise() %>%
                  mutate(Value = get_value(cur_data(), archetype_list)) %>%
                    ungroup() -> joinedBoy
    

    joinedBoy$Group <- factor(joinedBoy$Group, levels = levelsVec)


    p2 <- ggplot(joinedBoy, aes(x=Value, y=Delta, color=Var)) +
                geom_point(aes(shape = Sex.y)) +
                geom_segment(aes(xend=Value, yend=0)) +
                geom_hline(yintercept=0, color="#494949d2", linetype="dashed") +
                stat_smooth(method="glm", se=TRUE,alpha=0.3) +
                facet_wrap(~Group+Var,ncol=length(unique(joinedBoy$Var))) +
                scale_colour_manual(values=c("darkblue", "darkred","darkgreen","orange"))+
                theme_classic() +
                labs(title="", x="Value at State 1", y="Delta (Value at State 2 - Value at State 1)",color="State") +
                guides(shape = guide_legend("Sex")) 
   
  return(list(p,p2))
}


lowDimPlottingFunc <- function(dataMatNumeric, dataMatMeta, archeTypeMat,archetypeName,dataSetName,name="", contrasts,archeTypeMatrix) {
    library(umap)

    # Create PCA embedding from data
    pca <- prcomp(dataMatNumeric, center = TRUE, scale. = TRUE)
    pca_df <- as.data.frame(pca$x)
    pca_df$State <- dataMatMeta$State

    # Map archetypes to PCA space using loading matrix
    archetypesPCA <- as.data.frame(scale(archeTypeMat,center = pca$center,scale = pca$scale) %*% pca$rotation)
    archetypesPCA$Archs <- archetypeName

    # Create plot of PCA with archetypes specially marked
    pcaPlot <- ggplot() +
                geom_point(data = pca_df,aes(x = PC1, y = PC2, color = State)) +
                geom_text(data = archetypesPCA,aes(x = PC1, y = PC2, label = Archs)) + 
                theme_bw() + 
                ggtitle("PCA with Archetypes") +
                guides(color = guide_legend(override.aes = list(label = NULL))) 

    # Create the UMAP embedding from the data
    umap <- umap::umap(dataMatNumeric, n_neighbors = 7, min_dist = 0.6, metric = "euclidean")
    umap_df <- as.data.frame(umap$layout)
    umap_df$State <- dataMatMeta$State
    
    # Map archetypes to UMAP space
    archetypesUMAP <- predict(umap,archeTypeMat) %>% as.data.frame()
    archetypesUMAP$Archs <- archetypeName


    # Create plot of UMAP with archetypes specially marked
    umapPlot <- ggplot() +
                geom_point(data = umap_df,aes(x = V1, y =V2, color = State)) +
                geom_text(data = archetypesUMAP,aes(x = V1, y = V2, label = Archs)) + 
                theme_bw() + 
                ggtitle("UMAP with Archetypes") + 
                xlab("UMAP1") +
                ylab("UMAP2") +
                guides(color = guide_legend(override.aes = list(label = NULL))) 


    # Plot the plots together and save in a pdf

    plotsTegether <- gridExtra::grid.arrange(pcaPlot, umapPlot, ncol = 2)
    ggsave(file=paste0("./plots/",dataSetName,"/",dataSetName,"_LowDimWithArchetypes",name,".pdf"),plotsTegether,width = 10,height = 5)


    # Plot the same things expcept filtering for different constrasts

    contrastsNames <- sapply(contrasts, function(x) paste0(rev(x), collapse = " vs. "))

    dfListPCA <- list()
    for(i in 1:length(contrasts)) {
        dfListPCA[[i]] <- pca_df[which(pca_df$State %in% contrasts[[i]]),]
        dfListPCA[[i]]$ExperimentGroup <- contrastsNames[i]
    }


    # Combine the dataframes
    masterDFpca <- do.call(rbind,dfListPCA)

    # Create plot of PCA with archetypes specially marked
    pcaPlotLarbe <- ggplot() +
                geom_point(data = masterDFpca,aes(x = PC1, y = PC2, color = State)) +
                stat_ellipse(data = masterDFpca,aes(x = PC1, y = PC2, color = State)) +
                facet_wrap(.~ExperimentGroup,ncol = 3) +
                geom_text(data = archetypesPCA,aes(x = PC1, y = PC2, label = Archs)) + 
                theme_bw() + 
                guides(color = guide_legend(override.aes = list(label = NULL))) 
    
    ggsave(file=paste0("./plots/",dataSetName,"/",dataSetName,"_LowDimWithArchetypesContrastsPCA",name,".pdf"),pcaPlotLarbe,width = 12,height = 7)

    # PLOT THE CONTRASTS FOR UMAP

    dfListUMAP <- list()

    for(i in 1:length(contrasts)) {
        dfListUMAP[[i]] <- umap_df[which(umap_df$State %in% contrasts[[i]]),]
        dfListUMAP[[i]]$ExperimentGroup <- contrastsNames[i]
    }


    # Combine the dataframes
    masterDFumap <- do.call(rbind,dfListUMAP)

    # Create plot of PCA with archetypes specially marked
    umapPlotLarbe <- ggplot() +
                geom_point(data = masterDFumap,aes(x = V1, y = V2, color = State)) +
                facet_wrap(.~ExperimentGroup,ncol = 3) +
                geom_text(data = archetypesUMAP,aes(x = V1, y = V2, label = Archs)) + 
                stat_ellipse(data = masterDFumap,aes(x = V1, y = V2, color = State)) +
                theme_bw() + 
                guides(color = guide_legend(override.aes = list(label = NULL))) +
                xlab("UMAP1") +
                ylab("UMAP2") 
    
    ggsave(file=paste0("./plots/",dataSetName,"/",dataSetName,"_LowDimWithArchetypesContrastsUMAP",name,".pdf"),umapPlotLarbe,width = 12,height = 7)

    # Plot the lows dims but colour each point by the value of each of the archetypes

    archtypesColnames <- colnames(archeTypeMatrix)[sapply(archeTypeMatrix, is.numeric)]


    newdf <- cbind(pca_df,archeTypeMatrix[,archtypesColnames]) 

    newdflonger <- newdf %>% 
                    pivot_longer(cols = all_of(archtypesColnames), names_to = "Archetypes", values_to = "Value")

    # Create plot of PCA with archetypes specially marked
    pcaPlotLarbe <- ggplot() +
                geom_point(data = newdflonger,aes(x = PC1, y = PC2, color = Value)) +
                scale_color_gradient(low = "lightgray", high = "red") +
                facet_wrap(.~Archetypes,ncol = 3) +
                geom_text(data = archetypesPCA,aes(x = PC1, y = PC2, label = Archs)) + 
                theme_bw() + 
                guides(color = guide_legend(override.aes = list(label = NULL))) 
    
    ggsave(file=paste0("./plots/",dataSetName,"/",dataSetName,"_LowDimArchetypesAssignmentPCA",name,".pdf"),pcaPlotLarbe,width = 12,height = 7)

    # Also do umap
    newdf <- cbind(umap_df,archeTypeMatrix[,archtypesColnames]) 

    newdflonger <- newdf %>% 
                    pivot_longer(cols = all_of(archtypesColnames), names_to = "Archetypes", values_to = "Value")


     # Create plot of PCA with archetypes specially marked
    pcaPlotLarbe <- ggplot() +
                geom_point(data = newdflonger,aes(x = V1, y = V2, color = Value)) +
                scale_color_gradient(low = "lightgray", high = "red") +
                facet_wrap(.~Archetypes,ncol = 3) +
                geom_text(data = archetypesPCA,aes(x = PC1, y = PC2, label = Archs)) + 
                theme_bw() + 
                guides(color = guide_legend(override.aes = list(label = NULL))) +
                xlab("UMAP1") +
                ylab("UMAP2") 
    
    ggsave(file=paste0("./plots/",dataSetName,"/",dataSetName,"_LowDimArchetypesAssignmentUMAP",name,".pdf"),pcaPlotLarbe,width = 12,height = 7)

    return(list(umap_df,archetypesUMAP,pca_df,archetypesPCA))

}





# A function that plot a heatmap of the archetypes
plotArchTypesMatrix <- function(archMats,name,features,dataSetName) {

                # Pivot the matrix to long format
                dftresz <- archMats %>% as.data.frame() %>% 
                    pivot_longer(cols = everything(), names_to = "x", values_to = "z") %>% 
                    mutate(y = rep(rownames(archMats),each=ncol(archMats)))

                # Plot the heatmap
                dftresz$x <- factor(dftresz$x, levels =features)
                ploa <- ggplot(dftresz, aes(x = x, y = y, fill = z)) + 
                geom_tile() +
                scale_fill_gradient(low = "white", high = "steelblue") +
                geom_text(aes(label = round(z, 2)), color = "black") + # Round to 2 decimal places
                labs(x = "", y = "") +	
                theme_bw() 

                # Save the heatmap
                ggsave(paste0("./plots/",dataSetName,"/",dataSetName,name,".pdf"), plot = ploa, device = "pdf", width = 10, height = 5)

}

# A function that plots a scatter plot matrix of the archetypes

plotScatterMatrixForArchetypes <- function(bestArchetypes,archMats = NA,dataMatNumeric,featurs,dataSetName,name) {
    pdf(paste0("./plots/",dataSetName,"/",dataSetName,name,".pdf"),width = 16,height = 16)
    par(mfrow=c(length(featurs),length(featurs)))
    for(i in 1:length(featurs)){
        for(j in 1:length(featurs)){
            # If diagonal plot text thats TesmerFeature[i]
            if(i == j) {
                plot.new()
                text(0.5,0.5,featurs[i])
            } else {
                if(is.matrix(archMats)) {
                    bestArchetypes$archetypes <- archMats[,c(i,j)]
                }
            xyplot(bestArchetypes, dataMatNumeric[,c(i,j)],xlim=c(0,1),ylim=c(0,1),cex = 0.5)
            }
        }
    }
    dev.off()
}

plotBoxPlotOfArchtypes <- function(ArchetypesWithGroups,contrasts,dataSetName,archNames) {

    contrastsNames <- sapply(contrasts, function(x) paste0(rev(x), collapse = " vs. "))

    # Create separate filtered dataframe containing only the relevant rows to each contrast and then combine them in one dataframe and do a long pivot
    
    # Create a list of dataframes
    dfList <- list()
    for(i in 1:length(contrasts)) {
        dfList[[i]] <- ArchetypesWithGroups[which(ArchetypesWithGroups$State %in% contrasts[[i]]),]
        dfList[[i]]$ExperimentGroup <- contrastsNames[i]
    }

    # Reverse the order of (cosmetic plotting reasons)
    contrastsRev <- lapply(contrasts, rev)

    # Combine the dataframes
    ArchetypesWithGroups <- do.call(rbind,dfList)

    pivotLongerArchetypes <- ArchetypesWithGroups %>%
                                pivot_longer(cols = archNames, names_to = "Archetype", values_to = "Value")
       
    
    pivotLongerArchetypes$State <- factor(pivotLongerArchetypes$State, levels = unique(unlist(contrastsRev)))
    pivotLongerArchetypes$ExperimentGroup <- factor(pivotLongerArchetypes$ExperimentGroup, levels = contrastsNames)
    pivotLongerArchetypes$Archetype <- factor(pivotLongerArchetypes$Archetype, levels = archNames )


    # Plot the boxplot
    archBoxPlot <- ggplot(pivotLongerArchetypes, aes(x = Archetype, y = Value, fill = State)) +
                geom_boxplot(outlier.alpha = 0, position = position_dodge(width = 0.75)) +
                geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size=1, alpha=0.3) +
                facet_wrap(.~ExperimentGroup,ncol = 3) +
                theme_classic()	+ 
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

    ggsave(paste0("./plots/",dataSetName,"/",dataSetName,"_ArchetypesBoxPlot.pdf"), plot = archBoxPlot, device = "pdf", width = 12, height = 8)

}

#############################################################################################
# # # # # # # # # # # # AutoArchtypes wrapper function # # # # # # # # # # # # # # # # # # # 
#############################################################################################

### Computes the autogenerated behavior archetypes from the observed data 
autoArchetypesAnaPlots <- function(dataMat,dataSetName,kInterest,featurs,contrasts,plotTern){
    library(archetypes)

    dataMatNumeric <- dataMat %>%
                        select(all_of(featurs))  

    dataMatMeta <- dataMat[,c("MouseType","Mouse","Sex","State")]

    # Compute the archetypes 50 reps over 10 k values	
    arhcetypal <- stepArchetypes(data = dataMatNumeric, 
                         k = 1:length(featurs), 
                         verbose = FALSE,
                         nrep = 25)

    bestArchetypes <- bestModel(arhcetypal)

    # Check if the plots directory exists if not create it
    if(!dir.exists(paste0("./plots/",dataSetName,"/"))){
        dir.create(paste0("./plots/",dataSetName,"/"))
    }

    # Plot the scree plot of the archetypes
    pdf(paste0("./plots/",dataSetName,"/",dataSetName,"_ArchetypesScreePlot.pdf"),width = 8,height = 8)
    screeplot(arhcetypal)
    dev.off()

    # Plot the bar plot of the k = kInterest archetypes
    pdf(paste0("./plots/",dataSetName,"/",dataSetName,"_ArchetypesBarPlot.pdf"),width = 8,height = 8)
    barplot(bestArchetypes[[kInterest]], dataMatNumeric, percentiles = TRUE)
    dev.off()

    # Plot the archtypes in a scatter plot
    plotScatterMatrixForArchetypes(bestArchetypes[[kInterest]],archMats = bestArchetypes[[kInterest]]$archetypes,dataMatNumeric,featurs,dataSetName,"AutoArchMatrixScatterPlot")

    # Plot the ternary plot of the k = 3
    # Get Alphas
    ternPlotData <- as.data.frame(coef(bestArchetypes[[kInterest]],'alphas'))
    colnames(ternPlotData) <- paste0("AT",1:kInterest)

    ternPlots <- list()
    if(plotTern & kInterest == 3) {
        # Create a List of Ternary Plots
        for(i in 1:length(contrasts)) {
                ternPlots[[i]] <- plotTernaryWithArrows(ternPlotData,
                                dataMatMeta,
                                contrasts[[i]],
                                1.1,
                                "State",
                                paste0(contrasts[[i]][2]," vs. ",contrasts[[i]][1]))	
        }
        plots <- ggtern::arrangeGrob(grobs = ternPlots,ncol = 2)
        ggsave(file=paste0("./plots/",dataSetName,"/",dataSetName,"_AutoArchetypesTernPlot.pdf"),plots,width = 16,height = 8)
    }

    # Plot the delta plots
    autoArchetypesWithGroups <- cbind(ternPlotData,dataMatMeta[,c("Mouse","State","Sex")])

    
    diffPlot <- create_diff_plots(autoArchetypesWithGroups, paste0("AT",1:kInterest), "",contrasts)
    ggplot2::ggsave(paste0("./plots/",dataSetName,"/",dataSetName,"AutoArchtypesDeltaPlots.pdf"),device = "pdf",diffPlot[[1]],width = 12,height = 5)
    ggplot2::ggsave(paste0("./plots/",dataSetName,"/",dataSetName,"AutoArchtypesDeltaPlotsPoint.pdf"),device = "pdf",diffPlot[[2]],width = 10,height = 10)



    ## PLOT A BOXPLOT
    boxplots <- plotBoxPlotOfArchtypes(autoArchetypesWithGroups,contrasts,dataSetName,paste0("AT",1:kInterest)) 

    # Plot the archetypes in low dim space
    representations <- lowDimPlottingFunc(dataMatNumeric, dataMatMeta,
                        bestArchetypes[[kInterest]]$archetypes,
                        paste0("AT",1:kInterest),
                        dataSetName,
                        "",
                        contrasts,
                        autoArchetypesWithGroups)

    return(list("ternPlot" = ternPlots,"diffPlot" = diffPlot,"boxplot" = boxplots,"alphasKinterst" = as.data.frame(coef(bestArchetypes[[kInterest]],'alphas'),
    "representations" = representations)))

}


#####################################################################
# # # # # # # # Predifined Archetypes Functions # # # # # # # # # # #
#####################################################################

ArchAnalysisWithPredefined <- function(archMats,
                                        dataMat,
                                        dataSetName,
                                        featurs,
                                        contrasts,
                                        useObsVals) {

    # Separate dataframes (cuz dumb)
    dataMatNumeric <- dataMat %>%
                        select(all_of(featurs))  

    dataMatMeta <- dataMat[,c("MouseType","Mouse","Sex","State")]

    # Check if the plots directory exists if not create it
    if(!dir.exists(paste0("./plots/",dataSetName,"/"))){
        dir.create(paste0("./plots/",dataSetName,"/"))
    }

    # Check if the user wants to adjust the archtypes to observed values (with quantile adjustment as to avoid outlier skew)
    if(useObsVals) {
        archMatObsVals <- matrix(0,nrow = nrow(archMats), ncol = ncol(archMats),dimnames = list( rownames(archMats),colnames(archMats)))
        for(i in 1:nrow(archMats)) {
            for(j in 1:ncol(archMats)) {
                archMatObsVals[i,j] <- closeVal(archMats[i,j],dataMatNumeric[,colnames(archMats)[j]])
            }
        }
        archMatsOld <- archMats
        archMats <- archMatObsVals
    }


    library(archetypes)
    archetypes(dataMatNumeric,k=nrow(archMats)) -> bestArchetypes

    plotScatterMatrixForArchetypes(bestArchetypes,archMats,dataMatNumeric,featurs,dataSetName,"PredefArchMatrixScatterPlot")

    ## Compute the alphas for each mouse with grid search func
    #gridSearchArches <- archtypesGridSearch(archeMat = archMats,
    #                            resolution = resolu,
    #                            data = dataMatNumeric,
    #                            mouseData = dataMatMeta,
    #                            lossFunction = lossFunc,
    #                            scaling = scale,
    #                            returnHeatmap = makeLossPlots,
    #                            takeSqrt = takeSq,
    #                            dataSetName = dataSetName)

    Alphas <- archetypesNNLS(archeMat = archMats, data = dataMatNumeric)

    AlphasWGroups <- cbind(Alphas[[1]],dataMatMeta[,c("MouseType","Mouse","State","Sex")])

    plotArchTypesMatrix(archMats,"ArchetypesHeatmap",featurs,dataSetName)

    if(useObsVals) {
            trash <-lowDimPlottingFunc(dataMatNumeric, dataMatMeta,
                    archMatsOld,
                    rownames(archMats),
                    dataSetName,"OrgVals",
                    contrasts,
                    AlphasWGroups)

        plotArchTypesMatrix(archMatsOld,"ArchetypesHeatmapIntededValues",featurs,dataSetName)	
    }

    # Plot the archetypes in low dim space
    representations <- lowDimPlottingFunc(dataMatNumeric, dataMatMeta,
                    archMats,
                    rownames(archMats),
                    dataSetName,"",
                    contrasts,
                    AlphasWGroups)



    return(list(AlphasWGroups,"error" = Alphas[[2]],"representations" = representations))
    }


plotArchetypalAnalysisPredefined <- function(plotTern = TRUE,
                                   contrasts,
                                   AlphasWGroups,
                                   dataSetName,
                                   archMats
                    ) {

     # If makes sense (i.e. there are three archetypes) plot the ternary plots
    ternPlots <- list()
    if(plotTern) {
        # Create a List of Ternary Plots
        for(i in 1:length(contrasts)) {
            ternPlots[[i]] <- plotTernaryWithArrows(AlphasWGroups[,c("Circle","GO","Uncertainty")],
                            AlphasWGroups[,c("MouseType","Mouse","State","Sex")],
                            contrasts[[i]],
                            1.1,
                            "State",
                            paste0(contrasts[[i]][2]," vs. ",contrasts[[i]][1]))	
        }
        plots <- ggtern::arrangeGrob(grobs = ternPlots,ncol = 2)
        ggsave(file=paste0("./plots/",dataSetName,"/",dataSetName,"_PredefArchetypesTernPlot.pdf"),plots,width = 16,height = 8)
    } 

    # PLot the delta plots (both delta between state and starting delta dependance on starting value)

    diffPlot <- create_diff_plots(AlphasWGroups, rownames(archMats), "",contrasts)
    ggplot2::ggsave(paste0("./plots/",dataSetName,"/",dataSetName,"PredefinedArchtypesDeltaPlots.pdf"),device = "pdf",diffPlot[[1]],width = 12,height = 8)
    ggplot2::ggsave(paste0("./plots/",dataSetName,"/",dataSetName,"PredefinedArchtypesDeltaPlotsPoint.pdf"),device = "pdf",diffPlot[[2]],width = 10,height = 12)

    # Plot boxplot of the archetypes
    archboxs <- plotBoxPlotOfArchtypes(AlphasWGroups,contrasts,dataSetName,rownames(archMats)) 
    

    return(list("TernPlots" = ternPlots, "DeltaPlots" = diffPlot, "BoxPlots" = archboxs))

}

