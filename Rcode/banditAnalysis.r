# This file contains a set of functions that are used to analyze and plot the bandit data

library(tidyverse)

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
    textPanel1 <- function(x = 0.5, y = 0.5, txt, cex, font) {   
		text(x, y, txt, cex = cex, font = font,col="#b300ca")
        }


# Function to avarage the data over replicated experiments
averageReplicates <- function(data) {
    data %>% 
        group_by(MouseType,State,Mouse,Sex) %>%
        summarise_if(is.numeric,mean,na.rm=TRUE) %>%
        ungroup() %>%
        as.data.frame() %>%
        return()
}


# Create a PCA plotting function
plotPCAMouse <- function(data,columns,name,umapNeigbours,umapMinDist,CreateSDPlots,addCols = c()) {
    # Remove any rows with NA values
    dataTrim <- data[complete.cases(data[,columns]),]
    # Print how many rows were removed
    print(paste0("Removed ",nrow(data)-nrow(data[complete.cases(data[,columns]),])," rows with NA values."))
    # Get the labels for the data
    dataLabs <- dataTrim[,1:6]
    dataLabs$Mouse <- as.character(dataLabs$Mouse)
    dataLabs <- dataLabs[,colnames(dataLabs) != "NTrials"]
    # Select the columns to use for the PCA
    dataPCA <- dataTrim[,columns]
    dataPCAallcols <- dataTrim
    # Remove columns with only one unique value
    dataPCA <- dataPCA[,sapply(dataPCA, function(x) length(unique(x)) > 1)]
    # Plot the Correlation matrix to see how the features correlate
    group <- as.factor(dataLabs$State) %>% as.numeric()
    cols <- RColorBrewer::brewer.pal(n=length(unique(group)), name="Set3")[group]

    # Additional Columns
    extraCols <- c(columns,addCols)

    # Create folder if it doesnt exist
    dir.create((paste0("./plots/",name)))

    # Create a legend
    # create a data frame with unique groups and their associated colors
    dftemp <- data.frame(group = dataLabs$State, color =cols) %>% distinct()

    pdf(paste0("./plots/",name,"/Legend.pdf"),width = 4,height = 6)
    # plot an empty plot just for the legend
    plot(1, 1, type="n", ann=FALSE, axes=FALSE)

    # add a legend
    legend("center", legend = dftemp$group, fill = dftemp$color)
    dev.off()



    
    #dataPCA <- SecondHalfCollapsedFastedOnly[,TesmerFeatures]
    pdf(paste0("./plots/",name,"/_PairCorPlot.pdf"),width = ncol(dataPCA)/1.5,height = ncol(dataPCA)/1.5)
    #pdf(paste0("./plots/secondHalfFastedOnly.pdf"),width = ncol(dataPCA)/1.5,height = ncol(dataPCA)/1.5)
    graphics::pairs(dataPCA,
                    gap=0, 
                    lower.panel=graphics::panel.smooth, 
                    upper.panel=panel.cor, 
                    diag.panel=hist.panel,
                    labels=paste0("\n ",gsub("(.{12})", "\\1\n",colnames(dataPCA))),
                    text.panel = textPanel1,
                    pch=20,
                    col=cols,
                    font.labels=2,
                    cex.labels=0.9,
                    label.pos=1
                    ) 
    dev.off()

    # Create a hetmap which shows corrlations of catigorical variables and numeric varibles

    # List of character variables
        # Get the names of all non-numeric variables and store in vector
        char_vars <- names(dataLabs)[sapply(dataLabs, is.character)]
        # Reomve the type var from the char_vars
        char_vars <- char_vars[char_vars != "Type" & char_vars != "MouseType"]

        # List of numeric variables
        num_vars <- colnames(dataPCA)

        num_vars <- c(num_vars,addCols)


        # Function to calculate -log(p-value) for anova test
        calculate_pvalues <- function(df, char_var, num_vars){
        pvalues <- sapply(num_vars, function(num_var){
            fit <- aov(df[[num_var]] ~ df[[char_var]])
            #print(summary(fit))
            pvalue <- summary(fit)[[1]][["Pr(>F)"]][1]
            if(is.na(pvalue)) {
            return(NA)
            } else {
            return(-log10(pvalue))
            }
        })
        return(pvalues)
        }

        # Apply the function to each character variable
        pvalues_list <- lapply(char_vars, calculate_pvalues, df = dataPCAallcols,
        num_vars = num_vars)

        # Combine the results into a data frame
        pvalues_df <- do.call(rbind, pvalues_list)
        row.names(pvalues_df) <- char_vars

        # Plot heatmap
        pvalues_df %>% 
        as_tibble(rownames = "char_var") %>% 
        gather(num_var, value, -char_var) -> temp

        temp$num_var <- factor(temp$num_var,levels = columns)
        
        temp %>% 
        ggplot(aes(x = num_var, y = char_var, fill = value)) +
        geom_tile() + 
        ylab("Categorical Variable") +
        xlab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        scale_fill_gradient(low = "white", high = "red", name = "-log10(p-value)") -> pval_heatmap
        ggplot2::ggsave(paste0("./plots/",name,"/_catVarAnovaHeatmap.pdf"),device = "pdf",pval_heatmap,width = 6,height = 6)




    if(CreateSDPlots){
    # Plot a heatmap of SD for each var and each mouse stratified by state
    #View(cbind(dataLabs,dataPCA))
    varianceMat <-  cbind(dataLabs,dataPCA) %>%
                    group_by(Sex,Mouse,State) %>%
                    summarise_if(is.numeric,sd,na.rm = TRUE)
    #print(varianceMat)
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

    pdf(paste0("./plots/",name,"/_SDPerMouseAndStateHeatMapPlot.pdf"),width = 10,height = 10)
    #print(varianceMat)
    pheatmap::pheatmap(varianceMat[,columns],
                        annotation_row = annotationsRow,
                        show_rownames = FALSE,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        main = "Scaled Standard Deviation Heatmap")
    dev.off()

    }

        # Plot a heatmap of SD for each var by state
    varianceMat <-  cbind(dataLabs,dataPCA) %>%
                    group_by(State) %>%
                    summarise_if(is.numeric,sd,na.rm = TRUE)
    
    varianceMatGrps <- varianceMat %>% ungroup() %>% select(State)

    varianceMat <- varianceMat %>% 
                            ungroup() %>%
                            mutate_if(is.numeric,scales::rescale) %>%
                            as.data.frame()

    rownames(varianceMat) <- paste0(varianceMatGrps$State)

    annotationsRow <- data.frame(State = as.factor(varianceMatGrps$State),
                                row.names = rownames(varianceMat))

    pdf(paste0("./plots/",name,"/_SDPerStateHeatMapPlot.pdf"),width = 6,height = 6)
    pheatmap::pheatmap(varianceMat[,columns],
                        annotation_row = annotationsRow,
                        show_rownames = TRUE,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        main = "Scaled Standard Deviation Heatmap")
    dev.off()



    # Plot a heatmap of SD for each var and each mouse 
    varianceMat <-  cbind(dataLabs,dataPCA) %>%
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

    pdf(paste0("./plots/",name,"/_SDperMouseHeatMapPlot.pdf"),width = 7,height = 7)
    pheatmap::pheatmap(varianceMat[,columns],
                        annotation_row = annotationsRow,
                        show_rownames = FALSE,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        main = "Scaled Standard Deviation Heatmap")
    dev.off()

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

    pcaPlot2 <-ggplot(data = dataLabs,
                        aes(x = pca$x[,1],
                            y = pca$x[,2],
                            color = as.factor(State))) +
                geom_point(size = 4) +
                labs(x = paste0( "PC1 (Variance Explained: ", signif(var_explained[1]*100,4),"%)"),
                    y = paste0( "PC2 (Variance Explained: ", signif(var_explained[2]*100,4),"%)"), color = "Mouse",shape = "State") +
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

    
      pdf(paste0("./plots/",name,"/_LoadingsHeatmapPlot.pdf"),width = 5,height = 5)
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

        umapPlot2 <- ggplot(data = dataLabs,
                        aes(x = umap$layout[,1],
                            y = umap$layout[,2],
                            color = as.factor(State))) +
                    geom_point(size = 4) +
                    labs(x = "UMAP1",title = "UMAP",
                        y = "UMAP2", color = "Mouse",shape = "State") +
                    theme_bw()

    # Save the ggplots to pdf
    ggplot2::ggsave(paste0("./plots/",name,"/_ScreePlot.pdf"),device = "pdf",screePlot, width = 6, height = 6)
    ggplot2::ggsave(paste0("./plots/",name,"/_PCAPlot.pdf"),device = "pdf",pcaPlot, width = 8, height = 6)
    ggplot2::ggsave(paste0("./plots/",name,"/_PCAPlotStateCol.pdf"),device = "pdf",pcaPlot2, width = 8, height = 6)
    ggplot2::ggsave(paste0("./plots/",name,"/_PCABiPlotGroupedMicePC12.pdf"),device = "pdf",pcaBiPlotMice, width = 8, height = 6)
    ggplot2::ggsave(paste0("./plots/",name,"/_PCABiPlotGroupedMicePC34.pdf"),device = "pdf",pcaBiPlotMice2, width = 8, height = 6)
    ggplot2::ggsave(paste0("./plots/",name,"/_PCABiPlotGroupedStatePC12.pdf"),device = "pdf",pcaBiPlotState, width = 8, height = 6)
    ggplot2::ggsave(paste0("./plots/",name,"/_PCABiPlotGroupedStatePC34.pdf"),device = "pdf",pcaBiPlotState2, width = 8, height = 6)
    ggplot2::ggsave(paste0("./plots/",name,"/_UMAPPlot.pdf"),device = "pdf",umapPlot, width = 8, height = 6)
    ggplot2::ggsave(paste0("./plots/",name,"/_UMAPPlotStateCol.pdf"),device = "pdf",umapPlot2, width = 8, height = 6)


    return("Done")
}
     

