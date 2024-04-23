library(shiny)
library(shinydashboard)
library(leaflet)
library(leaflegend)
library(maps)
library(mapproj)
library(tidyverse)
library(leaflegend)
library(bslib)
library(rsconnect)
library(readr)
library(ggplot2)

# Load data ----
edna_data = suppressWarnings(read.csv("Data/winter_spring_historical - 3.csv", sep = ";", fileEncoding = "UTF-8"))
# jitter (add small variability) on lat and lon to avoid oevrlapping
edna_data$lat <- jitter(edna_data$lat, factor = 0.0001)
edna_data$lon <- jitter(edna_data$lon, factor = 0.0001)




shinyApp(
  ui = navbarPage("Fjord eDNA Biodiversity",
                  theme = bs_theme(
    # Controls the default grayscale palette
    bg = "#FFFFFF",
    fg = "#060606", # font color
    "navbar-bg" = "#428BCA",
    # Controls the accent (e.g., hyperlink, button, etc) colors
    primary = "#008BBC",
    secondary = "#F0F0F0",
    "input-border-color" = "#EA80FC"),
                  
                  
                  
                  #############################
                  # Setup the information panel
                  #############################
                  tabPanel(title = "Project information",
                           # style of the text
                           # padding-left control left margin
                               tags$style(type='text/css', 'body { overflow-y: scroll;}'), # make scrolling possible
                               tags$style(type="text/css", "h1 {width: 100%; text-align:center; font-size: 30px; padding-left:30px; padding-right:30px; color: #428BCA}"),
                               tags$style(type="text/css", "h2 {width: 100%; text-align:center; font-size: 25px; padding-left:30px; padding-right:30px; color: #428BCA}"),
                               tags$style(type="text/css", "h3 {width: 100%; text-align:center; font-size: 20px; padding-left:30px; padding-right:30px}"),
                               h1("Welcome to the Fjord eDNA Biodiversity shiny app!"),
                           
                               h2("Why are fjords so important?"),
                               h3("Fjords play a key role in coastal ecosystem processes serving as important habitats for marine species,
                               including sensitive species like rockfish,eulachons, salmons and glass sponges. They provide critical feeding grounds for anadromous
                               fish during their early life stages and recent hypotheses suggest that fjords can also serve as thermal refuges
                               for temperature-sensitive fish, such as caplin. Despite these unique characteristics, British Columbia fjords
                               remain chronically understudied by western science and continue to receive limited conservation protection. 
                               This project aims to use environmental DNA (eDNA) to assess the distribution of fish in BC fjords to guide marine conservation strategies."),
                           
                               img(src = "fjords_important.jpg", height = 300, width = 440, style="display: block; margin-left: auto; margin-right: auto;"),
                           
                               h2("What is eDNA?"),
                               h3("Marine eDNA is organism DNA in the water from microbial cells, organisms' tissues, skin and scales, metabolic waste,
                                  or dissolved molecules. A sample of eDNA simply requires a collection of fixed volume of water and collection of
                                  associated eDNA onto a filter. Extraction and sequencing of this DNA can be used to investigate the taxonomic composition
                                  of whole marine communities ranging from invertebrates to fish and marine mammals."),
                           
                               h2("How to interpret the interactive map?"),
                               h3("The map shows the eDNA index value associated with each species detected by eDNA during our sampling events. This species-specific
                                  index ranges from 0 (no DNA detected for the species) to 100 (maximum proportion of DNA detected for the species across
                                  all samples). The index reflects changes in the relative biomass (or proportion) of the species
                                  but do not reflect absolute changes in biomass. Because of biais associated
                                  with DNA Amplification, the index values cannot be compared among different species, but can be used to assess spatial variations
                                  of biomass for a given species"),
                           
                               h2("Sampling plan"),
                               h3("We collected eDNA from water samples in spring, summer, fall and winter of 2023, in Indian Arm, Atl'ka7tsem/Howe Sound, Toba Inlet, Bute Inlet and Rivers Inlet.
                               The interactive map also include samples collected in summer 2019 and 2022 by the Hakai Institute.
                                  At each site, eDNA replicates were colllected at 4 depths: surface, 10m, 100m and near bottom."),
                               img(src = "Sampling_map_fjords.jpg", height = 336, width = 600, style="horizontal-align:left; margin-left: 50px; margin-right: auto;"),
                               img(src = "sampling_vertical.JPG", height = 336, width = 650, style="horizontal-align:right; margin-left: 200px; margin-right: auto;"),
                           
                           
                               h2("Do you want to learn more?"),
                               h3("Please don't hesitate to contact us for more informations:
                                  Loic Jacquemot (l.jacquemot@oceans.ubc.ca); Brian Hunt (b.hunt@oceans.ubc.ca)")
                  ),
                  
                  ##############################
                  # Setup the map panel
                  ##############################
                  tabPanel("Interactive map",
                           tags$style(type = "text/css", "#map {height: calc(100vh - 53px) !important;}"), # control height of the map
                           fillPage(
                           leafletOutput('map')),
                           tags$style(type = "text/css", ".container-fluid {padding-left:0px; padding-right:0px;}"), #control left and right margins of the map
                           tags$style(type = "text/css", ".navbar {margin-bottom: -20px}"), # control space between the navbar and the map
                           
                           ##########################################
                           # add a right panel with vertical profiles
                           ##########################################
                           absolutePanel(top = 80, right = 50,
                                         style = "background-color: white;
                         opacity: 0.90;
  	                     padding: 1px 10px 1px 10px; # control size of the plot within the window - haut droite bas gauche
  	                     margin: 100;
  	                     box-shadow: 0pt 0pt 6pt 0px rgba(61,59,61,0.48);
  	                     border-radius: 10pt" # control angle of border
                                         ,
                                         
                                         # add some text          
                                         h3("Section of selected Fjord"), # text level 2
                                         
                                         # Fjord checkbox panel
                                         selectInput("fjord", 
                                                     label = "Fjord:",
                                                     choices = unique(edna_data$fjord),
                                                     selected = c("Bute")),
                                         # add vertical plot
                                         plotOutput("vertical", height = 200, width = 450)
                           ),
                           
                           ###################################################
                           # add a left panel with logo and selective panels
                           ###################################################
                           absolutePanel(bottom = 8, left = 20,
                                         # add the PEL and Hakai logo
                                         img(src = "PEL_LOGO_single.png", height = 140, width = 140),
                                         img(src = "Hakai_logo.png", height = 60, width = 140),  
                                         style = "color: white ; padding: 1px 10px 1px 10px",

                                         # select species panel          
                                         selectInput("species", 
                                                     label = "Choose a species to display",
                                                     choices = unique(edna_data$species),
                                                     selected = "Pacific herring (Clupea pallasii)"),
                                         
                                         # Year checkbox panel
                                         #selectInput("Year", 
                                         #                   label = "Year:",
                                         #                   choices = unique(edna_data$Year),
                                         #                   selected = c("2023")),
                                         
                                         # Season checkbox panel
                                         checkboxGroupInput("season", 
                                                            label = "Season:",
                                                            choices = unique(edna_data$season),
                                                            selected = c("May 2023", "February 2023")),
                                         
                                         # Depth slider panel
                                         sliderInput("depth", 
                                                     label = "Depth range (m)",
                                                     min = min(edna_data$depth),  
                                                     max = max(edna_data$depth), 
                                                     value = range(edna_data$depth))
                           )
                  )
  ),
  
  
######################################################################################################################
# Server logic ----
######################################################################################################################
  
  
  server = function(input, output, session) { 
    
    
    Sys.setlocale("LC_ALL", "C")
    
    #############################
    # create the reactives
    #############################
    
    # create the reactive for the map
    map_filter = reactive({
      edna_data %>%
        filter(species %in% input$species)%>%
        filter(depth >= input$depth[1] & depth <= input$depth[2])%>%
        filter(season %in% input$season)
      #filter(Year %in% input$Year)
    })
    
    # create the reactive for the plot
    plot_filter = reactive({
      edna_data %>%
        filter(species %in% input$species)%>%
        filter(fjord %in% input$fjord)%>%
        filter(season %in% input$season)
      #filter(Year %in% input$Year)
    })
    
    # create a color gradient
    pal <- colorNumeric(palette = c("white", "orange", "red")
                        , domain = c(0,100)
    )
    
    factorPal <- colorFactor(c("white", "orange", "red")
                             , domain = c(0,100)
    )
    #################
    # create the map
    #################
    output$map <- renderLeaflet({
      leaflet(edna_data) %>%
        addTiles() %>%
        addProviderTiles(providers$Esri.WorldImagery)%>% # suppose to add depth
        fitBounds(~min(lon)-0.5, ~min(lat)-0.5, ~max(lon)+0.5, ~max(lat)+0.5)%>%
        
        # control circles
        addCircles(data = map_filter(),
                   lat =  ~lat,
                   lng =  ~lon,
                   color = ~pal(indexvalue_100),
                   weight = ~indexvalue_100/2.5, # 2.5 to fit with legend size
                   label = ~station_depth,
                   popup = ~paste(indexvalue_popup, Date, sep = "   "))%>%
        
        # control circles legend
        addLegendSize(pal = pal,
                      shape = 'circle',
                      fillOpacity = .5,
                      opacity = 0,
                      values = c(-1,100),
                      baseSize = 20,
                      breaks = 10,
                      orientation = "horizontal",
                      title = 'eDNA Index value',
                      position = 'bottomright',
                      data = map_filter())
      #addLegend("bottomright",
      #          pal = pal,
      #          values = edna_data$indexvalue_100,
      #          title = "eDNA index value",
      #          opacity = 1)
    })
    
    ###################################
    # create the vertical section plot
    ###################################
    output$vertical <- renderPlot({
      ggplot(data = plot_filter(), aes(x = lat, y = depth)) + #prepare le mapping, sur quelles variables on travaille
        geom_point(aes(size = indexvalue_100, fill = indexvalue_100), shape = 21, limits = c(0,100))+
        scale_size_area(limits = c(0, 100), name = "eDNA Index Value")+
         scale_fill_gradient2(low = "white", mid = "orange",  high = "red" , midpoint = 50, limits = c(0,100))+
        labs(x="Latitude", y = "Depth (m)")+
        scale_y_reverse() +
        scale_x_continuous(breaks = plot_filter()$lat,
                           labels = plot_filter()$Site_ID,
                           minor_breaks = NULL)+
       theme_bw()+
        theme(
          legend.title = element_text(size=1),
          #panel.grid.minor = element_blank(),
          #panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          legend.position="bottom")

    })
    
    
    }
)


# see here to deploy on a html page : https://www.shinyapps.io/admin/#/dashboard
# NB: does not work well when integrated into OneDrive
# rsconnect::deployApp('C:/Users/loicj/Documents/02_Shiny_apps/census-app1 - winter-spring-historical - 3')