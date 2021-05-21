#############################################################
# R Shiny app for visualizing and exploring how temperature #
# influences insect dynamics under climate change           #
#############################################################

# Attach libraries
library(shiny)
library(shinythemes)
library(tidyverse)
library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

##################
# USER INTERFACE #
##################


# Define IU
ui <- shinyUI(fluidPage(
    
    # Application title
    titlePanel("Temperature effects on insect population dynamics"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            h4("Current climate:"),
            sliderInput("MeanTemp",
                        "Mean habitat temperature (°C):",
                        min = 0,
                        max = 40,
                        value = 20,
                        step = 0.1),
            sliderInput("TempFluc",
                        "Seasonal temperature fluctuations (difference between warmest and coldest month):",
                        min = 0,
                        max = 20,
                        value = 4,
                        step = 0.1),
            h4("Climate change:"),
            sliderInput("MeanT.incr",
                        "Increase in mean temperature over the next 100 years:",
                        min = 0,
                        max = 10,
                        value = 0,
                        step = 0.1),
            sliderInput("Fluc.incr",
                        "Increase in seasonal temperature fluctuations over the next 100 years:",
                        min = 0,
                        max = 20,
                        value = 0,
                        step = 0.1),
            h4("Starting population densities:"),
            numericInput("startA",
                         "Adults:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step = 0.01),
            numericInput("startJ",
                         "Juveniles:",
                         min = 0,
                         max = 1,
                         value = 0.01,
                         step = 0.01),
            h4("Plot options:"),
            numericInput("xmax",
                         "Length of simulation (years):",
                         min = 0,
                         max = 10e6,
                         value = 100,
                         step = 1),
            numericInput("ymax",
                         "Y-axis maximum:",
                         min = 0.1,
                         max = 3,
                         value = 1.5,
                         step = 0.1)
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            #tabsetPanel(
            #    tabPanel("Model dynamics",
                         fluidRow(
                             column(4, plotOutput("vital_repro")),
                             column(4, plotOutput("vital_devel")),
                             column(4, plotOutput("vital_mort"))
                         ),
                         fluidRow(
                             column(12, plotOutput("POPDYN"))
                #         )
                #),
                #tabPanel("Model description",
                #         withMathJax(includeMarkdown("ModelDescription.md")))
            )
        )
    )
))


##########
# SERVER #
##########


# define temperature conversion functions
C2K <- function(x){
    x+273.15
}
K2C <- function(x){
    x-273.15
}
# set color scheme
J_color <- "dodgerblue2"
A_color <- "firebrick2"
# define fixed model parameters
params <- matrix(c(1,293,5,0.1,283,15000,5000,70000,278,295,0.1,5000,0.1,2500),
                 dimnames=list(c("r_Topt","T_opt","s","m_R","T_R","A_m","A_L","A_H","T_L","T_H","d_JR","A_dJ","d_AR","A_dA"),
                               c("value")))

# calculate vital rate plot data, make plot templates
vitalrates <- data.frame(TempC = seq(0, 40, by = 0.1)) %>%
    mutate(TempK = C2K(TempC)) %>%
    mutate(Repro = params["r_Topt",]*exp(-((TempK-params["T_opt",])^2)/(2*(params["s",])^2)),
           Devel = params["m_R",]*(TempK/params["T_R",])*
               ((exp(params["A_m",]*((1/params["T_R",])-(1/TempK))))/
                    (1+exp(params["A_L",]*((1/params["T_L",])-(1/TempK)))+exp(params["A_H",]*((1/params["T_H",])-(1/TempK))))),
           Mort_J = params["d_JR",]*exp(params["A_dJ",]*((1/params["T_R",])-(1/TempK))),
           Mort_A = params["d_AR",]*exp(params["A_dA",]*((1/params["T_R",])-(1/TempK))))
# reproduction
t_vital_repro <- ggplot(vitalrates, aes(x = TempC, y = Repro))+
    geom_line(size = 2, color = J_color)+
    scale_x_continuous(name = "Temperature (°C)", limits = c(0, 40))+
    scale_y_continuous(name = "Reproduction", limits = c(0, 1))+
    scale_color_hue()+
    theme_bw()+
    theme(axis.title = element_text(size = 24, face = "bold"),
          axis.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16),
          panel.grid = element_blank())
# development
t_vital_devel <- ggplot(vitalrates, aes(x = TempC, y = Devel))+
    geom_line(size = 2, color = J_color)+
    scale_x_continuous(name = "Temperature (°C)", limits = c(0, 40))+
    scale_y_continuous(name = "Development", limits = c(0, 0.2))+
    scale_color_hue()+
    theme_bw()+
    theme(axis.title = element_text(size = 24, face = "bold"),
          axis.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16),
          panel.grid = element_blank())
# mortality
t_vital_mort <- ggplot(vitalrates, aes(x = TempC))+
    geom_line(aes(y = Mort_J), size = 2, color = J_color)+
    geom_line(aes(y = Mort_A), size = 2, color = A_color)+
    scale_x_continuous(name = "Temperature (°C)", limits = c(0, 40))+
    scale_y_continuous(name = "Mortality", limits = c(0, 0.6))+
    scale_color_hue()+
    theme_bw()+
    theme(axis.title = element_text(size = 24, face = "bold"),
          axis.text = element_text(size = 18),
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16),
          panel.grid = element_blank())

#################
# Define server #
#################

server <- shinyServer(function(input, output) {

    #############
    # REACTIVES #
    #############
    
    # fecundity plot
    output$vital_repro <- renderPlot({
        t_vital_repro+
            geom_vline(xintercept = input$MeanTemp)
    })
    # development plot
    output$vital_devel <- renderPlot({
        t_vital_devel+
            geom_vline(xintercept = input$MeanTemp)
    })
    # mortality plot
    output$vital_mort <- renderPlot({
        t_vital_mort+
            geom_vline(xintercept = input$MeanTemp)
    })
    
    # population dynamics
    output$POPDYN <- renderPlot({
        time <- seq(0, input$xmax*365, by = 1)
        state <- c(J = input$startJ, A = input$startA)
        
        # calculate habitat temperature
        HabitatTemp <- 273.15 + (input$MeanTemp + input$MeanT.incr*time/(365*100)) + (input$TempFluc/2 + input$Fluc.incr/2*time/(365*100))*sin(2*pi*time/365 + 0)
        # time-series of temperature values
        signal <- as.data.frame(list(times = time, temp = rep(0, length(time))))
        signal$temp <- HabitatTemp
        # create interpolating function
        inter.func <- approxfun(signal, rule = 2 ) # rule = 2 sets any extrapolated point to the closest data extreme
        
        # define population dynamic model
        TempMod <- function (time, state, Tparams) {
            with(as.list(c(state, Tparams)), {
                # import interpolated temperature function
                TempK <-  inter.func(time)
                
                # calculate dynamic population parameters
                r <- params["r_Topt",]*exp(-((TempK-params["T_opt",])^2)/(2*(params["s",])^2))
                m <- params["m_R",]*(TempK/params["T_R",])*((exp(params["A_m",]*((1/params["T_R",])-(1/TempK))))/
                                                                (1+exp(params["A_L",]*((1/params["T_L",])-(1/TempK)))+exp(params["A_H",]*((1/params["T_H",])-(1/TempK)))))
                d_J <- params["d_JR",]*exp(params["A_dJ",]*((1/params["T_R",])-(1/TempK)))
                d_A <- params["d_AR",]*exp(params["A_dA",]*((1/params["T_R",])-(1/TempK)))
                
                # population dynamics
                dJdt = r*A*exp(-A)-m*J-d_J*J
                dAdt = m*J-d_A*A
                return(list(c(dJdt, dAdt), signal = TempK))
            })
        }
        
        # project population dynamics
        popdyn <- as.data.frame(ode(func = TempMod, y = state, parms = params, times = time)) %>%
            gather("stage","popdens",2:3)
        # gather model variables (time, temp, life stage, density)
        model.output = popdyn[seq(0, dim(popdyn)[1], by=1), ]
        model.output = gather(model.output, key=Variable, value=Output, -time)
        
        # draw the population dynamics figure
        ggplot(popdyn, aes(x = time, y = popdens, group = stage, color = stage))+
            geom_line(size = 3)+
            scale_x_continuous(name = "Time (years)", limits = c(0, input$xmax))+
            scale_y_continuous(name = "Population Density", limits = c(0, input$ymax))+
            scale_color_manual(name = "Stage", breaks = c("J", "A"), labels = c("Juveniles", "Adults"), values = c(J_color, A_color))+
            theme_bw()+
            theme(axis.title = element_text(size = 24, face = "bold"),
                  axis.text = element_text(size = 18),
                  legend.title = element_text(size = 18, face = "bold"),
                  legend.text = element_text(size = 16),
                  legend.position = "top",
                  panel.grid = element_blank())
        
        # draw habitat temperature figure
        ggplot(model.output[model.output$Variable %in% c("signal"), ], aes(x=time, y=Output, color=Variable)) + 
            geom_line(size = 3) +
            scale_color_manual(values=c("signal"="#d1495b")) + 
            scale_x_continuous(name = "Time (years)", limits = c(0, input$xmax))+
            scale_y_continuous(name = "Temperature (C)", limits = c(input$MeanTemp - input$Fluc.incr/2, input$MeanTemp + input$TempFluc/2 + input$Temp.incr + input$Fluc.incr/2)) +
            theme_bw()+
            theme(axis.title = element_text(size = 24, face = "bold"),
                  axis.text = element_text(size = 18),
                  legend.title = element_text(size = 18, face = "bold"),
                  legend.text = element_text(size = 16),
                  legend.position = "top",
                  panel.grid = element_blank())
        
    })
    
})


###########
# RUN APP #
###########

shinyApp(ui = ui, server = server)
