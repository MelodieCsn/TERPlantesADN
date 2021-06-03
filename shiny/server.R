
# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  

    output$distPlot <- renderPlot({

        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')

    })
    
    # el <- reactive({
    #   if(input$element == "Na"){
    #     el <- "Na"
    #   }
    #   return (el)
    # })
    
    # output$selected_var <- renderImage({ 
    #   paste("You have selected", input$element)
    # })
    
        output$myImage <- renderImage({
          # A temp file to save the output.
          # This file will be removed later by renderImage
          outfile <- tempfile(fileext = '.png')
          
          # Generate the PNG
          png(outfile, width = 400, height = 300)
  
        
        if (is.null(input$element))
          {return(NULL)}
        
        if (input$element == "Na") {
          return(list(
            src = "images/var.PNG",
            contentType = "image/png",
            alt = "Face"
          ))
        }
      })
        # Loading files
        statgen_0kb <- read.table("data/genes_0kb_maf02.csv")
        statgen_1kb <- read.table("data/genes_1kb_maf02.csv")
        mlmm_0kb <- read.table("data/MLMM_genes_0kb.csv")
        mlmm_1kb <- read.table('data/MLMM_genes_1kb.csv')
        statgen_5kb <- read.table("data/genes_5kb_maf02.csv")
        statgen_10kb <- read.table("data/genes_10kb_maf02.csv")
        mlmm_5kb <- read.table("data/MLMM_genes_5kb.csv")
        mlmm_10kb <- read.table('data/MLMM_genes_10kb.csv')
        
        # Creating function
        removena <- function(df, first_row = 5){
          filtered_df <- df[rowSums(is.na(df[,first_row:ncol(df)])) != ncol(df[,first_row:ncol(df)]), ]
          return(filtered_df)
        }
        
        overval <- function(df, first_row = 5,val = 1, button = FALSE){
          if(button == T){
            temp_df <- df
            temp_df$sum1 <- rowSums(df[,first_row:ncol(df)], na.rm = T)
            temp_df <- filter(temp_df,sum1 > val)
            return(temp_df[,1:ncol(temp_df)-1])}
          else {return(df)}
        }
       # NewDF <- MAINDF[,sort(c(col.num, col.num - 1))]
        # Rending table
        datasetInput <- reactive({
          if(input$distance == "0kb"){if(input$package == "statgen"){statgen_0kb%>%
              select(c("chr","annot","type","name",unlist(input$element_candidat)))%>%
              removena(first_row = 5)%>%
              overval(first_row = 5,val = 1, button = input$commun)}
            else if(input$package == "mlmm"){mlmm_0kb%>%
                select(c("chr","annot","type","name",unlist(input$element_candidat)))%>%
                removena(first_row = 5)%>%
                overval(first_row = 5,val = 1, button = input$commun)}}
          else if(input$distance == "1kb"){if(input$package == "statgen"){statgen_1kb%>%
              select(c("chr","annot","type","name",unlist(input$element_candidat)))%>%
              removena(first_row = 5)%>%
              overval(first_row = 5,val = 1, button = input$commun)}
            else if(input$package == "mlmm"){mlmm_1kb%>%
                select(c("chr","annot","type","name",unlist(input$element_candidat)))%>%
                removena(first_row = 5)%>%
                overval(first_row = 5,val = 1, button = input$commun)}}
          else if(input$distance == "5kb"){if(input$package == "statgen"){statgen_5kb%>%
              select(c("chr","annot","type","name",unlist(input$element_candidat)))%>%
              removena(first_row = 5)%>%
              overval(first_row = 5,val = 1, button = input$commun)}
            else if(input$package == "mlmm"){mlmm_5kb%>%
                select(c("chr","annot","type","name",unlist(input$element_candidat)))%>%
                removena(first_row = 5)%>%
                overval(first_row = 5,val = 1, button = input$commun)}}
          else if(input$distance == "10kb"){if(input$package == "statgen"){statgen_10kb%>%
              select(c("chr","annot","type","name",unlist(input$element_candidat)))%>%
              removena(first_row = 5)%>%
              overval(first_row = 5,val = 1, button = input$commun)}
            else if(input$package == "mlmm"){mlmm_10kb%>%
                select(c("chr","annot","type","name",unlist(input$element_candidat)))%>%
                removena(first_row = 5)%>%
                overval(first_row = 5,val = 1, button = input$commun)}}
          
        })
        output$candidat <- renderTable({
          datasetInput()
          })
        output$downloadData <- downloadHandler(
          filename = function() {
            paste(input$package,input$distance, input$format, sep = "")
          },
          if(input$format == ".csv"){
            content = function(file) {
              write.csv(datasetInput(), file, row.names = FALSE)
            }}
          else if(input$format == ".txt"){
            content = function(file) {
              write.xlsx(datasetInput(), file, asTable = T)
            }}
        )
        
        ####CARTE####
        access <- read.csv("data/accessions_regMap.csv")
        cluster <- read.csv2("data/3cluster_regmap.csv")
        
        # Merging the dataframe
        
        cluster$accession <- cluster$genotype
        total <- merge(cluster, access,by="accession")
        
        # By original name
        
        access$genotype <- access$original_n..
        total2 <- merge(cluster, access, by="genotype")
        
        # Add popup information
        
        total2 <- total2%>%mutate(popup_inf = paste(accession.y,"<br/>", original_n.., "<br/>", cluster,"<br/>", country))
        total2[,'cluster']<-factor(total2[,'cluster'])
        
        # Render map
        output$mymap <- renderLeaflet({
        pal = colorFactor(as.character(wes_palette("Darjeeling1")[c(1,2,3)]), 
                          total2$cluster)
        leaflet()%>%addProviderTiles(providers$Stamen.TonerLite)%>%
          addCircleMarkers(data = total2, lat = ~Latitude, lng = ~Longitude,opacity = 0.7, radius = ~3, popup = ~popup_inf, color = ~pal(cluster))
        })
        
        ######PostGWAS############
        
        #DataTransfo
        common0 <- merge(statgen_0kb, mlmm_0kb, by.x='name', by.y='name')
        common1 <- merge(statgen_1kb, mlmm_1kb, by.x='name', by.y='name')
        common5 <- merge(statgen_5kb, mlmm_5kb, by.x='name', by.y='name')
        common10 <- merge(statgen_10kb, mlmm_10kb, by.x='name', by.y='name')
        
        #Rendering Venn
        output$venn <- renderPlot({
          if(input$post_distance == "0kb"){
            draw.pairwise.venn(area1 = nrow(statgen_0kb),
                               area2 = nrow(mlmm_0kb),
                               cross.area  = nrow(common0),
                               fill = c("#fc1505", "#096be3"),
                               category = c("statgen", "MLMM"),
                               cat.cex =2,
                               cex = 3)
          }
          else if(input$post_distance == "1kb"){
            draw.pairwise.venn(area1 = nrow(statgen_1kb),
                               area2 = nrow(mlmm_1kb),
                               cross.area  = nrow(common1),
                               fill = c("#fc1505", "#096be3"),
                               category = c("statgen", "MLMM"),
                               cat.cex =2,
                               cex = 3)
          }
          else if(input$post_distance == "5kb"){
            draw.pairwise.venn(area1 = nrow(statgen_5kb),
                               area2 = nrow(mlmm_5kb),
                               cross.area  = nrow(common5),
                               fill = c("#fc1505", "#096be3"),
                               category = c("statgen", "MLMM"),
                               cat.cex =2,
                               cex = 3)
          }
          else if(input$post_distance == "10kb"){
            draw.pairwise.venn(area1 = nrow(statgen_10kb),
                               area2 = nrow(mlmm_10kb),
                               cross.area  = nrow(common10),
                               fill = c("#fc1505", "#096be3"),
                               category = c("statgen", "MLMM"),
                               cat.cex =2,
                               cex = 3)
          }
          
        })
        ### only for Zn Need to be adapted
        snpZn <- readRDS("data/snpZn.RDS")
        regmap <-readRDS("data/regmap.RDS")
        g<-list()
        for (row in (1:nrow(snpZn))){
          name <- snpZn$`snpZn$newColName`[row]
          liste <- strsplit(snpZn$nonzeros[row], " ")
          #liste = list(snpZn$nonzero2[row])
          print(name)
          reg <- regmap
          reg$ecotype <- as.character(regmap$ecotype)
          #reg$Statut <- ifelse(is.element((regmap$ecotype),liste, "with", "without"))
          reg$Statut <- ifelse((regmap$ecotype) %in% liste[[1]], "Present", "Absent")
          g[[row]] <- ggplot(reg, aes(x=as.factor(Statut), y=Zn_change) ) + 
            geom_point(aes(color=as.factor(Statut)), size=2, alpha=0.8, position = position_jitterdodge(), show.legend = FALSE) +
            geom_boxplot(aes(color=as.factor(Statut)), alpha=0.2, lwd = 1.2,show.legend = FALSE) + 
            xlab(paste0("SNP name : ",name))+
            ylab("Variation de la concentration en Zn sous fort Co2 (en %)")+
            ggtitle("Effet de la présence d'un SNP sur le phénotype considéré")+
            theme_bw()+
            theme(text = element_text(size=10))
          #plot(g)
          output$snpeffect <- renderPlot(multiplot(g[[1]],g[[2]],g[[3]], g[[4]], g[[5]], g[[6]],
                                                   #g[[7]],g[[8]],g[[9]],g[[10]],
                                                   cols=3))
        }
})


