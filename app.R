library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker) 
library(BiocManager)
library(SummarizedExperiment)
library(readr)
library(gridExtra)
library(tidyverse)
library(DT)
library(reshape2)
library('fgsea')
library(pheatmap) 

#options(shiny.maxRequestSize=30*1024^2)
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington’s Disease and neurologically normal individuals"),
    
    tabsetPanel(
      ###Sample Tab
      tabPanel("Samples",
               sidebarLayout(
                 sidebarPanel(fileInput("file1", "Choose CSV File", accept = ".csv")),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Samples", tableOutput("tablesummar")),
                     tabPanel("Table", dataTableOutput("tablesam")),
                     tabPanel( "Plot",plotOutput("samhisto")),
                     
                   )
                 )
                 
                 
               )
               
               
               ),
      
      ###Counts Tab
      tabPanel("Counts",sidebarLayout(
        sidebarPanel(fileInput("fileco", "Choose CSV File", accept = ".csv",),sliderInput("vari", "Select the percentile of variance", min = 0, max = 100, value = 15),
                     sliderInput("nonzero", "Select the number of nonzero samples", min = 0, max = 70, value = 15),width = 3
                     
                     ),
        
        mainPanel(
          
          column(width = 12,tabsetPanel(
            #summary
            tabPanel("Samples", tableOutput("tableco")),
            
            #Scatter plots
            tabPanel("Diagnostic Scatter", plotOutput("scatterco")),
            
            #HeatMap
            tabPanel("Heatmap", plotOutput("heatmapco")),
            
            
            #PCA
            tabPanel( "PCA", 
                      column(width = 3, radioButtons("pcx", "Choose PC for the x-axis",
                                                     choices = c("PC1","PC2","PC3","PC4","PC5","PC6"),
                                                     selected = "PC1"),
                             radioButtons("pcy", "Choose PC for the y-axis",
                                          choices = c("PC1","PC2","PC3","PC4","PC5","PC6"),
                                          selected = "PC2"),
                             #submitButton(text = "Select", icon = icon("sink")),
                             
                             
                             
                             ),
                     
                      column(width = 9,plotOutput("pcaplot")),
                      ),
          
          ),
          )
        )
      )
               ),
      
      ###DE tab
      tabPanel("DE", 
               sidebarLayout(
                 sidebarPanel(fileInput("filede", "Choose CSV File", accept = ".csv")),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Table", dataTableOutput("table")),
                     tabPanel("Plot",
                              sidebarLayout(
                                sidebarPanel(
                                             h4("A volcano plot can be generated with log2 fold-change on the x-axis and p-adjusted on the y-axis."),
                                             
                                             radioButtons("button1", "Choose the column for the x-axis",
                                                          choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"),
                                                          selected = "log2FoldChange"),
                                             radioButtons("button2", "Choose the column for the y-axis",
                                                          choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"),
                                                          selected = "padj"),
                                             colourInput("col1", "Base point color", "#22577A"),
                                             colourInput("col2", "High point color", "#FFCF56"),
                                             
                                             sliderInput("magnitude",
                                                         "Select the magnitude of the p adjusted coloring:",
                                                         min = -30,
                                                         max = 1,
                                                         value = -15),
                                            
                                             #submitButton(text = "Plot", icon = icon("sink")),
                                             #helpText("When you click the button above, you should see",
                                             #         "the output below update to reflect the value you",
                                             #         "entered at the top:"),
                                             #verbatimTextOutput("value")
                                ),
                                mainPanel(plotOutput("volcano")
                                )
                              )
                              ),
                   )
                 )
               )

               
               ),
      ###GESEA tab
      tabPanel("GSEA",
               sidebarLayout(
                 sidebarPanel(fileInput("filegsea", "Choose CSV File", accept = ".csv")),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Top Results", 
                              sidebarLayout(
                                sidebarPanel(sliderInput("pvalres",
                                                         "Top results adjusted by p-value, 10^(x)",
                                                         min = -20,
                                                         max = 1,
                                                         value = -5,
                                                         )),
                                mainPanel(plotOutput("bargsea"))
                              )

                              ),
                     tabPanel("Table", 
                              sidebarLayout(
                                sidebarPanel(
                                  
                                  sliderInput("pvalrestab","Adjust by p-value",min = -21,max = 1,value = -5),
                                  radioButtons("nespath", "Choose the NES Pathways",choices = c("all","positive", "negative"),selected = "all"),
                                  downloadButton('download',"Download the data"),
                                  ),
                                mainPanel(dataTableOutput("tablegsea"))
                              )
                              
                              ),
                     tabPanel( "Plot",
                               sidebarLayout(
                                 sidebarPanel(
                                   
                                   sliderInput("scatpval","Adjust by p-value",min = -21,max = 1,value = -5),
                                 ),
                                 mainPanel(plotOutput("scatgsea"))
                               )
                               
                               ),
                     
                   )
                 )
               )
               
               
      ),
      
      
    ),

    
    
)

##########################################################################################
# Define server logic required to draw a histogram
server <- function(input, output) {

  #observeEvent(input$pvalrestab, {print(paste0("You have chosen: ", input$pvalrestab))})
  ###Data loaded
  load_datasam <- reactive({
    validate(need(input$file1, "Please select a file"),
             )
    
    return(read.csv(input$file1$datapath, header = TRUE))
  })

  
  #Counts file
  load_dataco <- reactive({
    validate(need(input$fileco, "Please select a file"))
    return(read.csv(input$fileco$datapath, header = TRUE))
  })

  
  #DE file
  load_datade <- reactive({
    validate(need(input$filede, "Please select a file"))
    return(read.csv(input$filede$datapath, header = TRUE))
  })

  
  #GSEA file
  load_datagsea <- reactive({
    validate(need(input$filegsea, "Please select a file"))
    return(read.csv(input$filegsea$datapath, header = TRUE))
    
  })
  
  
  ######Functions
  
  ####Sample info
  
  ##Summary Table for IDs
  
  #helper that makes the last row
  helper <- function(elem){
    if (is.numeric(elem) == FALSE){dat <- paste0(unique(elem), sep = ", ")}
    else {
      dat <- list(round(mean(elem,na.rm=TRUE),1), " (+/-", round(sd(elem,na.rm=TRUE), 1), "), ")
    }
    return(dat)
  
  }
  draw_summary <- function(metadata2) {
    #create three columns
    metadata2 <- metadata2[, -c(0:1)]
    column_name <- base::colnames(metadata2)
    column_name<- gsub("_", " ", column_name)
    column_name<- gsub("[.]", " ", column_name)
    Type <- sapply(metadata2,type)
    md <- sapply(metadata2,helper)
    md <- sapply( md, paste0, collapse="")
    #create the table
    summar <- tibble(column_name,Type,md)
    summar <- as.data.frame(summar)
    summar$md = substr(summar$md, 1, nchar(summar$md)-2)

    return(summar)
  }
  
  #regular table
  draw_tablesam <- function(dataf) {
    return(dataf)
  }
  
  
  
  #makes histograms based on the parameters
  samp_histograms <- function(metadata2) {
    plot1 <- ggplot(metadata2, aes(x=age_of_death)) + geom_histogram()
    plot2 <- ggplot(metadata2, aes(x=pmi)) + geom_histogram()
    plot3 <- ggplot(metadata2, aes(x=RIN)) + geom_histogram()
    plot4 <- ggplot(metadata2, aes(x=age_of_onset)) + geom_histogram()
    plot5 <- ggplot(metadata2, aes(x=duration)) + geom_histogram()
    plot6 <- ggplot(metadata2, aes(x=mrna.seq_reads)) + geom_histogram()
    plot7 <- ggplot(metadata2, aes(x=vonsattel_grade)) + geom_histogram(bins=2)
    plot8 <- ggplot(metadata2, aes(x=h.v_cortical_score)) + geom_histogram()
    plot9 <- ggplot(metadata2, aes(x=h.v_striatal_score)) + geom_histogram()
    #makes a grid of the plots above
    x <- grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9)
    return(x)
    
  }
  
  
  
  
  
  ####DE
  #Volcano plot
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    #construct the dataframe based on the parameters
    colvec = c(color1,color2)
    titlestr <- "P-adj < e"
    titlestr <- paste(titlestr,slider, sep = "")
    
    dataf1 <- dataf
    dataf1$difex <- "FALSE"
    dataf1$difex[dataf1$pvalue < 1*10^slider] <- "TRUE"
    
    #construct the plot
    plot1 <- ggplot(dataf1, aes(x=!!sym(x_name), y= -log10(!!sym(y_name)), col = difex  )  ) + 
      geom_point()+
      theme_light()+
      scale_colour_manual(values = colvec)+
      labs(color = titlestr) 
    
    
    return(plot1)
  }
  
  #Table
  draw_table <- function(dataf) {
    names(dataf)[names(dataf) == 'X'] <- 'gene'
    
    
    return(dataf)
  }
  
  
  
  ####Counts functions
  filter_var <- function(verse_counts, vari) {
    #remove columns that are not numbers
    dataf <- verse_counts[, -c(0:3)]
    x <- as_tibble(apply(dataf,1,var))
    #reconstruct the df with variance included
    dataf <- add_column(dataf,x)
    dataf <-cbind(verse_counts["...2"],dataf)
    #filter by var
    dataf<- dplyr::filter(dataf, quantile(value, (vari/100))>value)
    #dataf <-select(dataf,-value) 
    
    return(dataf)
  }
  
  
  filter_zero <- function(verse_counts, zero) {
    #dataf<-verse_counts
    dataf <- verse_counts[, -c(0:1)]
    dataf$nonzeros<-(ncol(dataf)-rowSums(dataf==0))
    dataf <-cbind(verse_counts["...2"],dataf)
    dataf<- dplyr::filter(dataf, zero < nonzeros)
    #dataf <-select(dataf,-zeros)
    
    return(dataf)
    
  }
  #counts table
  draw_tableco <- function(verse_counts,vari,zero) {
    #return(verse_counts)
    #dataf<-dataf[!is.na(foldvec)]
    
    dataf <- filter_var(verse_counts,vari)
    dataf <- filter_zero(dataf,zero)
    #return(dataf) 
    colnum <-ncol(verse_counts)-3
    Title <- c("Number of Samples", "Number of Genes","Genes passing filter","Percent genes passing filter", "Genes not passing filter","Percent genes not passing filter")
    Data <- c(colnum,nrow(verse_counts),nrow(dataf),nrow(dataf)/nrow(verse_counts)*100,(nrow(verse_counts)-nrow(dataf)),(nrow(verse_counts)-nrow(dataf))/nrow(verse_counts)*100 )
    Data <- as.integer(Data)
    summar <- tibble(Title,Data)
    return(summar)
  }
  #counts scatters
  plot_scatterco <- function(data,vari,zero) {
    #getting variance
    dataf <- data[, -c(0:3)]
    median <- apply(dataf, 1, median, na.rm=T)
    
    xval <- as_tibble(apply(dataf,1,var))
    dataf <- add_column(dataf,xval)
    #different than function because this is zeros, not non zeros
    dataf$zeros<-(rowSums(dataf==0))
    dataf$filt <- "Passed Filter"
    dataf$filt[quantile(dataf$value, 1-(vari/100))>dataf$value | zero>(ncol(dataf)-dataf$zeros-3)] <- "Not passed Filter"
    
    #-1 for the added columns
    #dataf$filt[zero<(ncol(dataf)-dataf$zeros-1)] <- "Not passed Filter"
    dataf <- add_column(dataf, median)
    
    
    colvec = c("#e0d9d1","#21130d" )
    #plot
    plot1 <- ggplot(dataf, aes(x=median,y=value,col = filt)) + geom_point() +ggtitle('Variance Scatter Plot') + scale_colour_manual(values = colvec)
    plot2 <- ggplot(dataf, aes(x=median,y=zeros,col = filt)) + geom_point() + ggtitle('Zero counts Scatter Plot') +scale_colour_manual(values = colvec)
    
    x <- grid.arrange(plot1,plot2)
    return(x)
    
  }
  # counts heatmap
  plot_heatmap <- function(data,vari,zero) {
    #dataf <- filter_var(data,vari)
    #dataf <- filter_zero(dataf,zero)
    #dataf <- melt(dataf)
    #dataf <- tibble(dataf["...2"],dataf["variable"])
    # ggplot(dataf, aes(x = ...2, y = variable, fill = value)) +
    #   geom_tile() +
    #   labs(title = "Correlation Heatmap",
    #        x = "Genes",
    #        y = "Sample")
    
    #filter data and remove the extra columns
    dataf <- filter_var(data,vari)
    dataf <- filter_zero(dataf,zero)
    dataf <-select(dataf,-nonzeros)
    dataf <-select(dataf,-value)
    #make first column into rownames
    dataf2 <- dataf[,-1]
    rownames(dataf2) <- dataf[,1]
    #scale
    dataf <- scale(dataf2)
    #plot
    x<-pheatmap(dataf, main = "Clustered Heatmap")
    return(x)
  }
  
  #counts PCA
  plot_pca <- function(data,pcx,pcy,vari,zero) {
    #cleaning df for pca
    dataf <- filter_var(data,vari)
    dataf <- filter_zero(dataf,zero)
    dataf <-select(dataf,-value)
    dataf <-select(dataf,-nonzeros)
    dataf <- dataf[, -c(0:1)]
    #pca
    pcar <- prcomp(t(dataf))
    newtib <- as_tibble(pcar$x)
    
    #variance
    summ <- summary(pcar)
    xlabel <- paste(pcx, ": ",(summ$importance[2,pcx]*100), "% variance", sep = "" )
    ylabel <- paste(pcy, ": ",(summ$importance[2,pcy]*100), "% variance", sep = "" )

    plot1 <- ggplot(newtib, aes(x=!!sym(pcx),y=!!sym(pcy))) +
      geom_point()+
      ggtitle('PCA plot') +
      xlab(xlabel) +
      ylab(ylabel)
    
    return(plot1)
    
  }
  
  
  
  
  ###GSEA ##########
  run_gsea <- function(dataf){
    dataf<-arrange(dataf,desc(log2FoldChange))
    xvec <- pull(dataf, symbol )
    logvec <- pull(dataf, log2FoldChange)
    foldvec <- setNames(logvec, xvec)
    foldvec<-foldvec[!is.na(foldvec)]
    hallmark_pathways_fgsea <- fgsea::gmtPathways('/Users/zachderse/Documents/591 project/R shiny project/h.all.v2023.2.Hs.symbols.gmt')
    fgsea_results <- fgsea(hallmark_pathways_fgsea, foldvec, minSize = 15, maxSize=500)
    fgsea_results <- fgsea_results %>% as_tibble()
    return(fgsea_results)
    
  }
  
  plot_bargsea <- function(dataf,pvalres){
    fgsea_results <- run_gsea(dataf)
    
    
    fgsea_results <- filter(fgsea_results, padj < 1*10^pvalres)
    #fgsea_results$color <- "Positive"
    #fgsea_results$color[fgsea_results$NES<0] <- "Negative"
    
    
    fgsea_results %>% arrange(padj)
    fgsea_results %>%
      mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
      ggplot() +
      geom_bar(aes(x=pathway, y=NES, fill = NES <  0), stat='identity') +
      scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
      theme_minimal() +
      theme(legend.position = "none")+
      ggtitle('FGSEA results for Hallmark MSigDB gene sets') +
      ylab('Normalized Enrichment Score (NES)') +
      xlab('') +
      coord_flip()
  }
  
  gsea_table <- function(dataf,pvalres,nespath){
    
    fgsea_results <- run_gsea(dataf)
    
    #fgsea_results <- unnest(fgsea_results, leadingEdge)
    
    #fgsea_results %>%
    #    mutate(new column = paste(unlist(leadingEdge), collapse=',')) %>%
    #     ungroup()
      
    #fgsea_results$leadingEdge <- as.character(unlist(fgsea_results$leadingEdge))
    fgsea_results$leadingEdge <- lapply(fgsea_results$leadingEdge,paste0,collapse=",")
    fgsea_results$leadingEdge <- as.character(fgsea_results$leadingEdge)
    
    
    
    fgsea_results <- filter(fgsea_results, padj < 1*10^pvalres)
    if (nespath == "positive") {
      fgsea_results <- filter(fgsea_results, NES > 0)
    }
    if (nespath == "negative") {
      fgsea_results <- filter(fgsea_results, NES < 0)
    }
    
    #fgsea_results %>%
    #  mutate(leadingEdge = map_chr(leadingEdge, ~ .[[1]] %>% str_c(collapse = ", ")))
    
    return(fgsea_results)
  }
  
  
  plot_scatgsea <- function(dataf,scatpval){
    
    fgsea_results <- run_gsea(dataf)
    #fgsea_results <- unnest(fgsea_results, leadingEdge)
    #fgsea_results <- filter(fgsea_results, padj < 1*10^pvalres)
    
    colvec = c("#e0d9d1","#21130d" )
    fgsea_results$filt <- "Not Passed Filter"
    fgsea_results$filt[fgsea_results$padj<1*10^scatpval] <- "Passed Filter"
    plot1 <- ggplot(fgsea_results, aes(x=NES,y=-log10(padj),col = filt)) + geom_point()  + ggtitle('FGSEA scatter plot') + scale_colour_manual(values = colvec)
    return(plot1)
  }
  
  
  reactgsea <- reactive({gsea_table(load_datagsea(),input$pvalrestab,input$nespath)})
  reactscatg <- reactive({plot_scatgsea(load_datagsea(),input$scatpval)})
  ######Output
  
  ##samples
  output$tablesummar <- renderTable({draw_summary(load_datasam())})
  output$tablesam <- DT::renderDataTable({draw_tablesam(load_datasam())})
  output$samhisto <- renderPlot({samp_histograms(load_datasam())})
  
  ##counts
  output$tableco <- renderTable({draw_tableco(load_dataco(),input$vari,input$nonzero)})
  output$scatterco <- renderPlot({ plot_scatterco(load_dataco(),input$vari,input$nonzero)})
  output$heatmapco <- renderPlot({ plot_heatmap(load_dataco(),input$vari,input$nonzero)})
  output$pcaplot <- renderPlot({ plot_pca(load_dataco(),input$pcx,input$pcy,input$vari,input$nonzero)})
  
  ##DE
  output$table <- DT::renderDataTable({draw_table(load_datade())})
  
  output$volcano <- renderPlot({ volcano_plot(load_datade(),input$button1,input$button2,
                                              input$magnitude, input$col1, input$col2)})
  
  #GESEA
  
  
  
  output$bargsea <- renderPlot({plot_bargsea(load_datagsea(),input$pvalres)})
  
  
  output$tablegsea <- DT::renderDataTable({reactgsea()})
  
  output$scatgsea <- renderPlot({reactscatg()})
  
  
  
  #download file
  
  output$download <- downloadHandler(
    filename <- function(){"filteredgsea.csv"},
    content <- function(fname){
      write.csv(reactgsea(), fname)
    }
  )
  
}
# Run the application 
shinyApp(ui = ui, server = server)
