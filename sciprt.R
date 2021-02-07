#Hacklytics
#Time Series forecasting of COVID19 new infection in GA
#
#
#
install.packages("MLmetrics")
install.packages("epitools")
install.packages("SimBIID")
library(epitools)
library(SimBIID)
library(MLmetrics)
library(coda)
install.packages("sf")
install.packages("raster")
install.packages("spData")
library(tmap)    # for static and interactive maps
library(leaflet) # for interactive maps
library(ggplot2) # tidyverse data visualization packag
library(EpiEstim)
library(shiny)
library(RcppXPtrUtils)
X <- read.csv(url("https://raw.githubusercontent.com/JieYingWu/COVID-19_US_County-level_Summaries/master/data/counties.csv"))
#Daily infections time series
Y <- read.csv(url("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"))
#daily deaths time series
#Z<- read.csv(url("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv"))
#merge daily infections and deaths with demographic data
XY <- merge(X,Y)
#XZ <- merge(X,Z)


XyFilt <- XY[,c(2,3,5,7,40,110,54,120,330,331,274,302,313,672:ncol(XY))]
#XzFilt <- XZ[,c(3,5,7,40,110,54,120,330,331,274,302,313,359:ncol(XZ))]
rm(X,Y,XY)
# remove rows with na (252 out of 3142)
XyFilt2 <- na.omit(XyFilt)
#XzFilt2 <- na.omit(XzFilt)

#datesThruSunDeath <- epidate(x = colnames(XzFilt[13:ncol(XzFilt)]), format = "X%m.%d.%y")$dates
datesThruThursInfec <- epidate(x = colnames(XyFilt[14:ncol(XyFilt)]), format = "X%m.%d.%y")$dates

XyFilt2 <- subset(XyFilt2, State == "NC")


GASocDist <- read.csv("~/Downloads/Unacast_Social_Distancing_Score_GA.csv")
library(outbreaker2)
#library(pomp)
library(projections)

#county ROs per week predictor
countyR0sWeek <- vector(mode="numeric",length = nrow(XyFilt2))
for (i in 1:nrow(XyFilt2)){
a <- as.numeric(as.matrix(as.vector(t(XyFilt2[i,14:ncol(XyFilt2)]))))
ex <- estimate_R(a,method = "parametric_si",
           config = make_config(list(mean_si = 3.96, std_si = 4.75)))
ex$dates <- datesThruThursInfec
ex$R$t_start <- datesThruThursInfec[ex$R$t_start]
ex$R$t_end <- datesThruThursInfec[ex$R$t_end]
countyR0sWeek[i] <- round(ex$R[nrow(ex$R),3],2)}

XYts <- ts(t(XyFilt2[,14:ncol(XyFilt)]))
colnames(XYts) <- XyFilt2$Area_Name
rownames(XYts) <- colnames(XyFilt2)[14:ncol(XyFilt2)]


datesThruThursInfec <- epidate(x = colnames(XyFilt[14:ncol(XyFilt)]), format = "X%m.%d.%y")$dates

NCcounties <- list()
for (cas in 1:100){
  casescounty <- vector(mode = "numeric", length = sum(XYts[,cas]))
  NCcounties[[cas]] <- rep(rownames(XYts) , XYts[,cas])
  names(NCcounties)[[cas]] <- colnames(XYts)[cas]
}
rownames(XYts) <- as.character(datesThruThursInfec)
NCcountyIncidence <- list()
for (cas in 1:100){
  i2 <- incidence(rep(rownames(XYts) , XYts[,cas]))
  NCcountyIncidence[[cas]] <- i2
  print(cas)
}

library(distcrete)

si <- distcrete("gamma", interval = 1L,
                shape = 1.5,
                scale = 2, w = 0)
NCprojection <- list()
NCprojNegBi <- list()
for (cas in 1:100){
  P <- project(NCcountyIncidence[[cas]], 
               R = countyR0sWeek[cas], si = si, n_days = 5)
  NCprojection[[cas]] <- P
  #names(NCprojection[[cas]]) <- colnames(XYts)[cas]

  #NP <- project(NCcountyIncidence[[cas]], 
  #               R = countyR0sWeek[cas], si = si,
  #               n_days = 5, model = "negbin")
  # NCprojNegBi[[cas]] <- NP
  # names(NCprojNegBi[[cas]]) <- colnames(XYts)[cas]
  print(cas)
    }
names(NCprojection) <- gsub(pattern = " County", replacement = "", x = names(NCprojection))



ui <- fluidPage(
  
  # Application title
  titlePanel("NC covid predictions"),
  
  # Top panel with county name
  verticalLayout(
    
    wellPanel(textOutput("cnty")),
    plotOutput("plot"),
    # the map itself
    mainPanel(
      leafletOutput("map")
    )
  )
)
server <- function(input, output) {
  
  output$map <- renderLeaflet({
    leaflet() %>% 
      addProviderTiles("Stamen.Toner") %>% 
      addPolygons(data = shape, 
                  fillColor = "aliceblue", 
                  color = "grey",
                  layerId = ~CNTY_ID)
  })
  
  # this is the fun part!!! :)
  observe({ 
    event <- input$map_shape_click
    nameo <-shape$NAME[shape$CNTY_ID == event$id]
    output$cnty <- renderText(nameo)
    output$plot <- renderPlot(plot(cumulate(NCprojection[[match(nameo,names(NCprojection))]])))
                                    
                                
  })
}

# Run the application 
shinyApp(ui = ui, server = server)






