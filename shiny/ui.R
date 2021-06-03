
library("openxlsx")
library("shiny")
library("dplyr")
library("ggplot2")
library("tidyr")
library("stats")
library("plyr")
library("plotrix")
library("shinydashboard")
library("plotly")
library("hrbrthemes")
library("viridis")
library("RColorBrewer")
library("ggridges")
library("VennDiagram")
library("leaflet")
library("rgdal")
library("wesanderson")
library("tidyr")
library("rlist")
library("shinycssloaders")
library("shinyWidgets")
library("tidyverse")
library(VennDiagram)
library(Rmisc)



dashboardPage(
  
  #Header
  dashboardHeader(title = "TER ADN"),
  #Sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem("Présentation", tabName = "Presentation"),
      menuItem("Statistiques descriptives", tabName = "Stats"),
      menuItem("ACP", tabName = "ACP"),
      menuItem("Clustering", tabName = "Clusterng"),
      menuItem("GWAS statgenGWAS", tabName = "StatgenGWAS"),
      menuItem("GWAS mlmm", tabName = "Mlmm"),
      menuItem("Gènes candidats", tabName = "Genes"),
      menuItem("Visualisations post GWAS", tabName = "postGWAS"),
      menuItem("Carte", tabName = "Carte")
      
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
                                /* logo */
                                .skin-blue .main-header .logo {
                                background-color: #7de193;
                                }

                                /* logo when hovered */
                                .skin-blue .main-header .logo:hover {
                                background-color: #61c878;
                                }

                                /* navbar (rest of the header) */
                                .skin-blue .main-header .navbar {
                                background-color: #7de193;
                                }

                                /* main sidebar */
                                .skin-blue .main-sidebar {
                                background-color: #7de193;
                                }

                                /* active selected tab in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                                background-color: #61c878;
                                }

                                /* other links in the sidebarmenu */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                                background-color: #7de193;
                                color: #000000;
                                }

                                /* other links in the sidebarmenu when hovered */
                                .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                                background-color: #61c878;
                                }
                                /* toggle button when hovered  */
                                .skin-blue .main-header .navbar .sidebar-toggle:hover{
                                background-color: #61c878;
                                }

                                /* body */
                                .content-wrapper, .right-side {
                                background-color: #FFFFF7;
                                }

                                '))),
    tags$style(type = "text/css", "#mymap {height: calc(100vh - 80px) !important;}"),
    tabItems(
      tabItem(tabName = "Presentation",
              
              sidebarLayout(
                sidebarPanel(
                  sliderInput("bins",
                              "Number of bins:",
                              min = 1,
                              max = 50,
                              value = 30)
                ),
                
                # Show a plot of the generated distribution
                mainPanel(
                  plotOutput("distPlot")
                )
              )
              
              
              
      ),
      tabItem(tabName = "Stats",
              #titlePanel("Répartition géographique"),
              # Sidebar panel for inputs ----
              radioButtons("element", "Choisissez l'élément à afficher : ",
                           choices= c("Carbone" = "C",
                                      "Zinc" = "Zn",
                                      "Fer" = "Fe" ,
                                      "Magnésium" = "Mg",
                                      "Sodium"="Na",
                                      "Cuivre" = "Cu",
                                      "Manganèse" = "Mn",
                                      "Azote" = "N"
                                      )),
              imageOutput("myImage")
              
      ),
      tabItem(tabName = "Carte",
              leafletOutput(outputId = "mymap")
              #titlePanel("Répartition géographique"),
              # Sidebar panel for inputs ----
              
              
      ),
      tabItem(tabName = "postGWAS",
        sidebarPanel(
          sidebarLayout(
            prettyRadioButtons(
              inputId = "post_distance",
              label = "Distance maximum entre le SNP et le gène", 
              choices = c("0kb","1kb","5kb","10kb"),
              icon = icon("check"), 
              bigger = TRUE,
              status = "info",
              animation = "jelly"
                  ),
            prettyRadioButtons(
              inputId = "post_element",
              label = "Phénotype considéré (Variation de la composition de l'élèment en eCO2/aCO2)", 
              choices = c("Zn", "C", "Cu", "Mg", "Mn","Na", "Fe","Composante 1 de l'ACP"="Comp1"),
              icon = icon("check"), 
              bigger = TRUE,
              status = "info",
              animation = "jelly"
                  )
                  
          )),
        mainPanel(plotOutput("venn"), plotOutput("snpeffect"))
              ),
      
      tabItem(tabName="Genes",
              p("Sélectionnez et télécharger une liste de gènes candidats selon plusieurs critères (distance au SNP
                identifiés, package utilisé, phénotype considéré"),
              sidebarLayout(
                sidebarPanel(
                  sliderTextInput(
                    inputId = "distance",
                    label = "Distance maximum entre le SNP identifié et le gène:", 
                    choices = c("0kb", "1kb", "5kb", "10kb")
                  ),
                  awesomeRadio(
                    inputId = "package",
                    label = "Méthode de GWAS utilisée:", 
                    choices = c("StatgenGWAS" = "statgen", "MLMM" ="mlmm", "Gènes en commun entre les deux méthodes"="both"),
                    selected = "statgen"
                  ),
                  pickerInput(
                    inputId = "element_candidat",
                    label = "Phénotype(s) considéré(s) (Variation de la composition de l'élèment en eCO2/aCO2)", 
                    choices = c("Zn", "C", "Cu", "Mg", "Mn","Na", "Fe","Composante 1 de l'ACP"="Comp1"),
                    selected = c("Zn", "C", "Cu", "Mg", "Mn","Na", "Fe","Composante 1 de l'ACP"="Comp1"),
                    options = list(
                      `actions-box` = TRUE), 
                    multiple = TRUE
                  ),
                  awesomeCheckbox(
                    inputId = "commun",
                    label = "Gènes identifiés pour au moins deux élèments (parmi ceux séléctionnés au dessus)", 
                    value = FALSE,
                    status = "info"
                  ),
                  radioGroupButtons(
                    inputId = "format",
                    label = "Format de sortie",
                    choices = c("Tableau (csv)" =".csv", 
                                "Tableau(xlsx)" = ".xlsx", "Liste (txt)"=".txt"),
                    justified = TRUE,
                    checkIcon = list(
                      yes = icon("ok", 
                                 lib = "glyphicon"))
                  ),
                  downloadButton("downloadData", "Download")),
              mainPanel(
                tableOutput('candidat')
              ))
              
      )
      
    )
    
  )
)

