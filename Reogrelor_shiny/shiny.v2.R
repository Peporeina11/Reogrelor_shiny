library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinycustomloader)
library(deSolve)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(MASS)
library(plyr)
library(compiler)
library(doParallel)

cl <- makePSOCKcluster(detectCores() - 1)

clusterEvalQ(cl, list(library(deSolve), library(foreach)))

registerDoParallel(cl)	

# UI ---- 
ui <- dashboardPage(
  title = "Rheogrelor Case Study",
  
  # Header
  dashboardHeader(
    title = h4("Rheogrelor Case Study", style = "font-family:'helvetica'; color:'aqua';"),
    titleWidth = 200
  ),
  
  # Sidebar
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Graphs", icon = icon("line-chart"), selected = TRUE, tabName = "pkpd_sim1")
    ),
    fluidPage(
      h5('', style = "padding:350px;"),
      hr(style = "border-color:gray;"),
      h5('Models by: A. Di Deo, P. Laddomada, U. Villani, O. Della Pasqua', style = "color:gray;"),
      h5('Created by: S. D\'Agate and A. Di Deo', style = "color:gray;")
    )
  ),
  
  # Body
  dashboardBody(
    tabItems(
      tabItem(tabName = "pkpd_sim1",
              fluidRow(
                # Plot 1 Card
                box(
                  title = "Rheogrelor Concentration",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 6,
                  withLoader(plotOutput("plotCONC", height = 400, width = 600),
                             type = "image",
                             loader = 'https://img.pikbest.com/png-images/20190918/cartoon-snail-loading-loading-gif-animation_2734139.png!bw700')
                ),
                
                # Plot 2 Card
                box(
                  title = "Platelet Aggregation",
                  status = "warning",
                  solidHeader = TRUE,
                  width = 6,
                  withLoader(plotOutput("plotPHT", height = 400, width = 600),
                             type = "image",
                             loader = 'https://img.pikbest.com/png-images/20190918/cartoon-snail-loading-loading-gif-animation_2734139.png!bw700')
                )
              ),
              
              # Inputs Card
              fluidRow(
                box(
                  title = "Simulation Parameters",
                  status = "info",
                  solidHeader = TRUE,
                  width = 12,
                  
                  sliderInput("n", HTML('Number of subjects'), min = 0, max = 100, value = 50, step = 1, ticks = F),
                  fluidRow(
                    column(6, numericInput("min_WT", HTML("Minimum weight [kg]"), value = 50)),
                    column(6, numericInput("max_WT", "Maximum weight [kg]", value = 120))
                  ),
                  sliderInput("interval_1", "Infusion duration [hours]", min = 1, max = 6, value = 1, step = 0.5, ticks = F),
                  fluidRow(
                    column(6, numericInput("dose_CNG_1", HTML("Initial Rheogrelor Bolo [μg/kg]"), value = 30)),
                    column(6, numericInput("dose_Infusion", HTML("Infusion Dose [μg/kg/min]"), value = 5))
                  ),
                  hr(style = "border-color:lightgray;"),
                  h5('Clinical Trial Design', style = "color:gray;"),
                  sliderInput("Samples", "Number of Samples", min = 1, max = 24, value = 1, step = 1, ticks = F),
                  submitButton("Simulate")
                )
              )
      )
    )
  )
)


# SERVER ----
server <- function(input, output, session) {
  
  
  mod1 <- reactive({
    
    sumfuncx <- function(x) {
      stat1 <-  median(x)
      stat2 <-  quantile(x, probs=0.05, names=F) 
      stat3 <-  quantile(x, probs=0.95, names=F)
      stat4 <-  length(x)
      result <- c("median"=stat1, "low"=stat2, "hi"=stat3, "n"=stat4)
      result
    }
    
    perc_resp <- function(x) {
      perc <- c('median' = sum(x)/length(x)*100, 'low' = NA, 'hi' = NA, 'n' = length(x))
      perc
      perc
    }
    
    #Create a dataframe with ID and parameter values for each individual	
    n <- input$n
    par.data <- seq(from = 1, to = n, by = 1)
    par.data <- data.frame(par.data)
    set.seed(02072024)
    WT <- runif(n, min = input$min_WT, max = input$max_WT)
    names(par.data)[1] <- "ID"
    
    #Define population values
    
    #PK Cangrelor
    TVCL_CNG = 39.4
    TVV1_CNG = 3.3
    
    #PD Platelet Inhibition
    TVKIN = 1
    TVKOUT = 0
    TVIC50 = 15.0
    
    
    #Define population parameter variability
    
    #PK Cangrelor
    ETACL_CNG <- rnorm(n, mean = 0, sd = sqrt(0.03))
    ETAV1_CNG <- rnorm(n, mean = 0, sd = sqrt(0.04))
    
    #PD Platelet Inhibition
    ETAIC50 = rnorm(n, mean = 0, sd = sqrt(0.22))
    
    #Simulate individual values	
    #PK Cangrelor
    par.data$CL_CNG = TVCL_CNG * (WT/70)^1 * exp(ETACL_CNG)
    par.data$V1_CNG = TVV1_CNG * (WT/70) * exp(ETAV1_CNG)
    
    #PD Platelet Inhibition
    par.data$KIN = TVKIN 
    par.data$KOUT = TVKOUT
    par.data$IC50 = TVIC50 * exp(ETAIC50)
    par.data$WT = WT
    
    TIME = sort(unique(seq(0,3,0.1)))
    
    platelet_mod <- function(t, S, parameters) {
      with(as.list(c(S, parameters)), {
        
        ### ODE system
        dSdt=vector(len=3)
        
        ### PK Cangrelor
        Ke_CNG = CL_CNG/V1_CNG
        
        # Infusion rate (0 before, during T_inf it's Dose_inf / T_inf, 0 after)
        R_inf <- ifelse(t <= interval_1, dose_Infusion*WT*60/interval_1, 0)
        
        # Central compartment
        dSdt[1] = R_inf - Ke_CNG * S[1]
        
        # AUC
        dSdt[2] = S[1]/V1_CNG
        
        # C/Css
        C_CNG = S[1]/V1_CNG
        Cavg_CNG = ifelse(S[2] == 0, yes = 1e-10, no = S[2]/t)
        
        ### PD Platelet Inhibition
        INH = C_CNG/(C_CNG+IC50)
        
        # Platelet compartment
        dSdt[3] = KIN * (1 - INH) - KOUT*S[3]
        
        list(dSdt)
      })
    }
    
    
    th.cmpf <- cmpfun(platelet_mod)
    #----------------------------------------------------------------------------------------	
    #Function for simulating concentrations for the ith patient
    simulate.conc <- function(par.data, dose_CNG_1, interval_1, dose_Infusion, treat_dur) {
      
      par <- c(par.data, dose_Infusion = dose_Infusion, interval_1 = interval_1)
      
      #Set initial conditions in each compartment
      S = c(
        S1 = dose_CNG_1*par.data$WT,
        S2 = 0,
        S3 = 1
      )
      
      events_dataset <- rbind(data.frame(var = c("S1"), time=0, value = dose_CNG_1*par.data$WT, method=c("add")),
                              data.frame(var = c("S1"), time=0, value = dose_Infusion*par.data$WT, method=c("add")))
      
      sim.data <- lsoda(y = S, 
                        times = TIME, 
                        func = th.cmpf, 
                        parms = par,
                        atol = 1e-3,
                        ynames = F,
                        rtol	= 1e-3,
                        events	= list(data=events_dataset))
      
      sim.data <- as.data.frame(sim.data)	
    }
    
    #Compile simulate.conc function	
    simulate.conc.cmpf <- cmpfun(simulate.conc)
    
    #Apply simulate.conc.cmpf function to each individual in par.data
    #ID is provided as a variable in which each value has simulate.conc.cmpf applied to
    #V1 and V2 are required to be here to preserve their values, so concentrations can be
    #calculated later
    suppressWarnings({
      sim.data <- ddply(par.data, .(ID, V1_CNG), simulate.conc.cmpf, 
                        dose_CNG_1 = input$dose_CNG_1,
                        dose_Infusion = input$dose_Infusion,
                        interval_1 = input$interval_1,
                        treat_dur = 12,
                        .parallel = TRUE)
    })
    
    sim.data$Cavg_CNG <- sim.data$S2/(sim.data$time+1e-10)
    sim.data$CONC_CNG <- sim.data$S1/sim.data$V1_CNG
    sim.data$PLAT_INH <- 1 * (1-(sim.data$CONC_CNG/(sim.data$CONC_CNG+23.2)))
    
    #sim.data$resp_DFX <- ifelse(test = ((sim.data$FER_DFX < 2500) | ((sim.data$FER_DFX-TVFERBASE)/TVFERBASE <= -0.2)), yes = 1, no = 0)
    
    
    suppressWarnings({
      statsCONC_DFX <- ddply(sim.data, 
                             .(time), function(sim.data) sumfuncx(sim.data$CONC_CNG), .parallel = TRUE)
      statsCONC_DFX$time <- (statsCONC_DFX$time)
      statsFER_DFX <- ddply(sim.data, 
                            .(time), function(sim.data) sumfuncx(1-sim.data$PLAT_INH), .parallel = TRUE)
      statsPERC_DFX <- ddply(sim.data, 
                             .(time), function(sim.data) perc_resp(sim.data$resp_DFX), .parallel = TRUE)
    })
    
    #Combine datasets
    list(
      conc = statsCONC_DFX,
      PLT  = statsFER_DFX
    )
    
    
  })
  
  output$plotCONC <- renderPlot({
    
    df <- mod1()$conc
    
    ggplot(df, aes(x = time*60, y = median)) +
      geom_line(size = 1.2, colour = "steelblue") +
      geom_ribbon(aes(ymin = low, ymax = hi),
                  fill = "steelblue", alpha = 0.3) +
      scale_x_continuous(breaks = seq(0, 180, 30)) +
      labs(y = "Reogrelor Plasma Concentration [ng/mL]", x = "Time Since Start of Infusion [min]") +
      theme_bw() +
      theme(
        plot.background = element_rect(fill = '#ecf0f4', colour = '#ecf0f4'),
        plot.title = element_text(color = "#0099f8", size = 16, face = "bold", hjust = 0.5),
      )
    
  })
  
  output$plotPHT <- renderPlot({
    
    df <- mod1()$PLT
    
    ggplot(df, aes(x = time*60, y = median)) +
      geom_line(size = 1.2, colour = "steelblue") +
      geom_ribbon(aes(ymin = low, ymax = hi),
                  fill = "steelblue", alpha = 0.3) +
      scale_x_continuous("\nTime [min]") +
      ylab("Rheogrelor concentration [ng/mL]") +
      geom_hline(yintercept = c(0.8, 0.9), linetype =2) +
      theme_bw() +
      theme(
        plot.background = element_rect(fill = '#ecf0f4', colour = '#ecf0f4')
      )
    
  })
  
}

### EXECUTE APP
shinyApp(ui, server)




