#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

source("core.R")

evidence <- list(prev = c(43L, 457L), sn = c(41L, 2L), sp = c(147, 310))

# Define UI for application that draws a histogram
ui <- fluidPage(
  
    # Application title
    titlePanel("ENBS calculator"),
    tabsetPanel(id="input_tabset",
      tabPanel("Introduction",HTML("
               <H3><B>Welcome to the Expected Net Benefit of Sampling (ENBS) calculator</B></H3><BR/>
               <P>This Web app helps you determine the optimal sample size for the external validation of your risk prediction model based on uncertainty around its net benefit (NB).</P>
               <P><B>To use this app, you need the following categories of information: </B><BR/>
                  1. Two 'exchange rates': one is the exchange rate between true and false positives, implied by the risk threshold, and another one between the sampling efforts and clinical benefit.<BR/>
                    
                  2. The performance of the test in terms of its sensitivity and specificity (alongside outcome prevalence) at the risk threshold of interest.<BR/>
                  
                  3. An estimate of the total expected usage of the model for decision making in the target population.<BR/>
                  
                  4. General setup for calculations (range of target sample size, number of Monte Carlo simulations).<BR/></P>
                  
                  <P> The main goal is to calculate ENBS. This is done via the following equation:</P>
                  
                  <P align='center'><B>ENBS = EVSI(z)*N-n/lambda</B>, where:</P><BR/>
                  
                  <P>
                  - EVSI is the Expected Value of Sample Information for the decision at the chosen risk threshold <BR/>
                  
                  - z is the risk threshold <BR/>
                  
                  - N is the total number of times the model will be used (EVSI*N is the 'Population EVSI') <BR/>
                  
                  - n is the sample size of the future validation study <BR/>
                  
                  - lambda is the sampling effort in true positive units <BR/>
                  </P>
                  
                  <P>ENBS can be used to determine the optimal sample size from a decision-theoretic perspective.</P>
                  <P>The nominal value of ENBS, which is in true positive units by default, can also be used to judge whether the external validation study of a given sample size is worth taking</P>
      ")),
      tabPanel("Exchange rates",
               fluidRow(
                 column(4,
                        sliderInput("z", label="Risk threshold (%)", min=0, max=100, value=2, width="100%")),
                 column(8, textOutput("z_desc"))
               ),
               hr(),
               fluidRow(
                 column(4, 
                        sliderInput("lambda", "Desired 'Number Needed to Study': The efforts erequired to procure how many observations is equal to the benefit of one true positive detection?", min=0, max=10000, value=100, step=1, width="100%")),
                 column(8, textOutput("lambda_desc"))
               ),
               actionButton("infer_lambda","Infer minimum NNS from the risk threshold"),
               textOutput("infer_lambda_desc")
      ),
      tabPanel("Model performance",
        HTML("Here we solicit your assessment of the model performance and your uncertainty around this asessment. 
             Currently, this is based on trivariate specification (prevalence, sensitivity, and specificity), with Beta distribution for modeling uncertainty for each component. In the future, other methods of uncertainty assessment will be added"),
        hr(),
        fluidRow(
          column(4,sliderInput("prev",label="Prevalence (%)", min=0, max=100, step=0.1, value=evidence$prev[1]/sum(evidence$prev)*100, width="100%")),
          column(4,numericInput("prev_n",label="Your estimate of prevalence is based on how many observations?", value=sum(evidence$prev))),
          column(4, textOutput("prev_dist"))
        ),hr(),
        fluidRow(
          column(4,sliderInput("sn",label="Sensitivity (%)", min=0, max=100, step=0.1, value=evidence$sn[1]/sum(evidence$sn)*100, width="100%")),
          column(4,numericInput("sn_n",label="Your estimate of sensitivity is based on how many observations?", value=sum(evidence$sn))),
          column(4, textOutput("sn_dist"))
        ),hr(),
        fluidRow(
          column(4,sliderInput("sp",label="Specificity (%)", min=0, max=100, step=0.1, value=evidence$sp[1]/sum(evidence$sp)*100, width="100%")),
          column(4,numericInput("sp_n",label="Your estimate of sensitivity is based on how many observations?", value=sum(evidence$sp))),
          column(4, textOutput("sp_dist"))
        )
      ),
      tabPanel("Population parameters",
        fluidRow(
          column(4,
            numericInput("N", label="Expected number of times the medical decision is to be made.", value=800000)),
          column(4, "Think of the prevalence of the condition, the number of years the model will remain in user (e.g., before the next revision), and the 'market penetration' of the model (e.g., which proportion of practitioners will consider using the model)")
        )
      ),
      tabPanel("Analysis setup",
        numericInput("n_min","Starting sample size", value=500),
        numericInput("n_max","Ending sample size", value=16000),
        sliderInput("n_step", "Sample size multiplier at each step", value=2, min=1, max=10, step=0.25),
        numericInput("n_sim", "Number of simulations", min=100, max=10^67, value=10^5)
      ),
      tabPanel("Results",
        actionButton("run", "Run"), actionButton("clear_results", "Clear results"),
        uiOutput("results")
      )
    )
)




# Define server logic required to draw a histogram
server <- function(input, output) 
{
  construct_evidence <- reactive(
  {
    evidence <- list(
      prev=c(round(input$prev/100*input$prev_n), round(input$prev_n-input$prev/100*input$prev_n)),
      sn=c(round(input$sn/100*input$sn_n), round(input$sn_n-input$sn/100*input$sn_n)),
      sp=c(round(input$sp/100*input$sp_n), round(input$sp_n-input$sp/100*input$sp_n))
    )
    evidence
  })
  
  evidence_changed <- reactive({
    list(input$prev,input$prev_n, input$sn, input$sn_n, input$sp, input$sp_n)
  })
  observeEvent(evidence_changed(),{
    evidence <- construct_evidence()
    
    make95CrI <- function(parms)
    {
      paste0("95%CI ", 
            round(qbeta(0.025,parms[1],parms[2])*100,1),
            "%-",
            round(qbeta(0.975,parms[1],parms[2])*100,1),"%"
            )
    }
    
    output$prev_dist <- renderText(paste0("prev~Beta(",paste0(evidence$prev,collapse=","),")\n", make95CrI(evidence$prev)))
    output$sn_dist <- renderText(paste0("sn~Beta(",paste0(evidence$sn,collapse=","),")\n", make95CrI(evidence$sn)))
    output$sp_dist <- renderText(paste0("sp~Beta(",paste0(evidence$sp,collapse=","),")\n", make95CrI(evidence$sp)))
  })
  
  observeEvent(input$z, {
    z <- input$z/100
    output$z_desc <- renderText(paste0("The risk threshold is the threshold on predicted risk at which point the care provider / patient is indifferent between using or not using the tool.\n
          A risk threhsold of ",input$z, "% indicates that the benefit of a true positive case is equal (in the opposite direction) to the harm of ",(1-z)/z," false positive cases."))
  })
  
  observeEvent(input$lambda, {
    output$lambda_desc <- renderText(paste("The Number Needed to Study (NNS) represents the trade-off between sampling efforts and clinical utility. For example, a NNS of 100 means the investigator believes the efforts required to procure 100 more samples is equal to the benefit of one true positive diagnosis.\n
          To choose this value, think of the context: How important the clinical event is? For example, for a catstrophic event like a stroke, one might justify procuring a sample in hundreds, or thoushands, to prevent one more event\n
          This is also related to the ease by which samples can be obtained. Will this external validation study be based on primary or secondary data collection? The former might require significantly more efforts. \n
          In thinking about sampling efforts, do not consider one-time costs and efforts required for conducting a validaiton study. Those one-time costs are not affected by sample size. Instead, think of  'incremental effort of procuring samples.
          "))
  })
  
  observeEvent(input$infer_lambda, {
    z <- input$z/100
    updateSliderInput(inputId="lambda", value=max(1, (1-z)/(z)))
    output$infer_lambda_desc <- renderText(paste("The inferred minimum NNS is based on the notion that in most contexts, it is justifiable to obtain one more observation insofar as that observation is expected to prevent one incorrect medical decision.\n 
                                                 Incorrect decisions can be in terms of not detecting an individual who will experience the event (missinge a true positive), or unncecesarily treating an individuals who will not experience the event (causing a false positive).\n"
                                                 ,ifelse(z<0.5, 
                                                    paste0("Because at risk threshold of ", input$z, "% each true positive detection is equal to", (1-z)/z, " false positives, and because the EVSI is in true positive units, the desired NNS should be at least", (1-z)/z,"."),
                                                    ""
                                                    ),
                                                  "In general, this reasoning results in NNS of max(1,(1-z)/z) where z is the clinical threshold."))
  })
  
  observeEvent(input$clear_results, {
    output$results=renderUI(HTML(""))
  })
  
  observeEvent(input$run, {
    evidence <- construct_evidence()
      
    z <- input$z/100
    n_sim <- input$n_sim
    
    n_stars <- c(0,as.integer(round(exp(seq(from=log(input$n_min), to=log(input$n_max), by=log(input$n_step))))))
    VoI <- EVSI_ag(evidence, z, n_sim=n_sim, future_sample_sizes=n_stars[-1])
    EVPI <- VoI$EVPI
    p_best <- VoI$p_best
    EVSIs <- c(0,VoI$EVSI)
    lambda <- input$lambda
    N <- input$N
    ENBS <- EVSIs*N-n_stars/lambda
    
    require("knitr")
    output$results <- renderUI(list(
      renderText(paste("EVPI=",EVPI)),
      hr(),
      renderTable(p_best),
      HTML(ifelse(max(p_best$p_best)>0.99, "<B style='color:red; font-weight:bold;'>In more than 99% of simulations the same stratgy had the same NB. This indicates there is not much uncertainty around this decision. VoI analysis might not be informative and might be degenrate.</B>","")),
      hr(),
      renderTable((data.frame("sample size"=as.integer(n_stars),
                                                      "EVSI"=format(EVSIs, nsmall=5),
                                                      "Population EVSI"=as.double(EVSIs*N),
                                                    "ENBS"=as.double(ENBS)))),
      renderPlot({
        plot(n_stars, EVSIs, type='l', xlab="Sample size of the future study", ylab="EVSI")
        y2 <- pretty(c(0,N*EVPI))
        axis(4, at=y2/N, labels=y2)
        mtext("Population EVSI", side = 4)
        lines(c(0,max(n_stars)), rep(EVPI,2), col="gray")
      }),
    
      renderPlot({
        plot(n_stars, ENBS, type='l', xlab="Sample size of the future study", ylab="Expected Net Benefit of Sampling")
        winner <- which.max(ENBS)
        if(winner!=1 & winner!=length(ENBS))
        {
          lines(c(n_stars[winner], n_stars[winner]), c(0, max(ENBS)), col='red')
          text(n_stars[winner]*1.1,max(ENBS)/2,paste0("Optimal sample size:",n_stars[winner]), col='red')
        }
        else
        {
          text(0,max(ENBS)/2,"Edge case!", col='red')
        }
      })
    ))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
