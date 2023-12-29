#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(evsiexval)

evidence <- list(prev = c(43L, 457L), se = c(41L, 2L), sp = c(147, 310))

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("ENBS calculator (CONFIDENTIAL - please do not share the link)"),
    tabsetPanel(id="input_tabset",
      tabPanel("Introduction",HTML("
               <H3><B>Welcome to the Expected Net Benefit of Sampling (ENBS) calculator for external validation studies of risk prediction models</B></H3><BR/>
               <P>This Web app helps you determine the optimal sample size for the external validation of your risk prediction model based on uncertainty around its net benefit (NB).</P>
               <P><B>To use this app, you need the following categories of information: </B><BR/>
                  1. Two 'exchange rates': one is the exchange rate between true and false positives, implied by the risk threshold, and another one between the sampling efforts and clinical benefit.<BR/>

                  2. The performance of the model in terms of its sensitivity and specificity at the risk threshold of interest (alongside outcome prevalence).<BR/>

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
                  <HR/>
                  <DIV style='display: flex;	 align-items: flex-end;  align-text: center'>App version 2023.12.20. For questions and bug reports contact msafavi@mail.ubc.ca</DIV>
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
             Currently, this is based on trivariate specification (outcome prevalence, sensitivity, and specificity), with Beta distribution for modeling uncertainty for each component. In the future, other methods of uncertainty assessment will be added"),
        hr(),
        selectInput("evidence_type", "What type of evidence you have?", choices=c("0.PLEASE SELECT",
                                                                                  "1.Independent beta performance",
                                                                                  "2.Summary performance measures",
                                                                                  "3.Draws from the posteroir distirbution")
                    ),
        uiOutput("evidence_inputs")
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
        sliderInput("n_step", "Sample size multiplier at each step", value=2, min=1.25, max=10, step=0.25),
        fluidRow(
          column(4,
            numericInput("n_sim_outer", "Number of outer-level simulations", min=100, max=10^6, value=10^5)),
          column(4,
            textOutput("n_sim_outer_note"))
        ),
        numericInput("n_sim_inner", "Number of inner-level simulations", min=100, max=10^6, value=10^5),
        textOutput("n_sim_inner_note"),
        p("Note: the computational difficulty of VoI methods differ, and the number of simulations should be chosen accordingly"),
        br(),
        p("For any non-trivial computations, please use the R package or the local version of this app on your computer")
      ),
      tabPanel("Results",
        actionButton("prelim_run", "Preliminary analysis"), actionButton("clear_results", "Clear results"),
        uiOutput("prelim_results"),
        uiOutput("results")
      )
    ),
)




# Define server logic required to draw a histogram
server <- function(input, output)
{
  construct_evidence <- reactive(
  {
    evidece_type <- substring(input$evidence_type,1,1)

    if(evidece_type=="1")
    {
      evidence <- list(
        type=1,
        params=list(
          prev=c(round(input$prev/100*input$prev_n), round(input$prev_n-input$prev/100*input$prev_n)),
          se=c(round(input$se/100*input$se_n), round(input$se_n-input$se/100*input$se_n)),
          sp=c(round(input$sp/100*input$sp_n), round(input$sp_n-input$sp/100*input$sp_n))
        )
      )
    }
    if(evidece_type=="2")
    {
      evidence <- list(
        type=2,
        params=list(
          prev=c(mean=input$prev[1]/100, upper_ci=input$prev[2]/100),
          cs=c(mean=input$cstat[1], upper_ci=input$cstat[2]),
          A=c(mean=input$cal_intercept/1, sd=input$cal_intercept_sd/1),
          B=c(mean=input$cal_slope/1, sd=input$cal_slope_sd/1)
        )
      )
    }
    evidence
  })


  make95CrIFromBeta <- function(a,b)
  {
    paste0("95%CI ",
           round(qbeta(0.025,a,b)*100,1),
           "%-",
           round(qbeta(0.975,a,b)*100,1),"%"
    )
  }

  observeEvent(input$z, {
    z <- input$z/100
    output$z_desc <- renderText(paste0("The risk threshold is the threshold on predicted risk at which point the care provider / patient is indifferent between treatment or no treatment.\n
          A risk threhsold of ",input$z, "% indicates that the benefit of a true positive case is equal (in the opposite direction) to the harm of ",(1-z)/z," false positive cases."))
  })

  observeEvent(input$lambda, {
    output$lambda_desc <- renderText(paste("The Number Needed to Study (NNS) represents the trade-off between sampling efforts and clinical utility. For example, a NNS of 100 means the investigator believes the efforts required to procure 100 more samples is equal to the benefit of one true positive diagnosis.\n
          To choose this value, think of the context: How important the clinical event is? For example, for a catstrophic event like a stroke, one might justify procuring a sample in hundreds, or thoushands, to prevent one more event\n
          This is also related to the ease at which samples can be obtained. Will this external validation study be based on primary or secondary data collection? The former might require significantly more efforts. \n
          In thinking about sampling efforts, do not consider one-time costs and efforts required for setting up a validaiton study. Those one-time costs are not affected by sample size. Instead, think of  'incremental' effort of procuring samples.
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

  observeEvent(input$evidence_type, {
    choice <- substring(input$evidence_type,1,1)
    output$n_sim_outer_note <- renderText("")
    output$evidence_inputs <- renderUI("")

    if(choice=="1")
    {
      output$evidence_inputs <- renderUI(list(
        fluidRow(
          column(4,sliderInput("prev",label="Expected outcome risk (outcome prevalence) (%)", min=0, max=100, step=0.1, value=evidence$prev[1]/sum(evidence$prev)*100, width="100%")),
          column(4,numericInput("prev_n",label="Your estimate of the average risk of the clinical outcome is based on how many observations?", value=sum(evidence$prev))),
          column(4, textOutput("prev_dist"))
        ),hr(),
        fluidRow(
          column(4,sliderInput("se",label="Sensitivity of the model at the chosen risk threshold (%)", min=0, max=100, step=0.1, value=evidence$se[1]/sum(evidence$se)*100, width="100%")),
          column(4,numericInput("se_n",label="Your estimate of sensitivity is based on how many observations?", value=sum(evidence$se))),
          column(4, textOutput("se_dist"))
        ),hr(),
        fluidRow(
          column(4,sliderInput("sp",label="Specificity the model at the chosen risk threshold (%)", min=0, max=100, step=0.1, value=evidence$sp[1]/sum(evidence$sp)*100, width="100%")),
          column(4,numericInput("sp_n",label="Your estimate of sensitivity is based on how many observations?", value=sum(evidence$sp))),
          column(4, textOutput("sp_dist"))
        )
      ))
      observeEvent(input$prev, {
        a <- input$prev/100*input$prev_n
        b <- input$prev_n-input$prev/100*input$prev_n
        output$prev_dist <- renderText(paste0("prev~Beta(", a,",", b ,") | ", make95CrIFromBeta(a,b)))
      }, autoDestroy=T)
      observeEvent(input$se, {
        a <- input$se/100*input$se_n
        b <- input$se_n-input$se/100*input$se_n
        output$se_dist <- renderText(paste0("se~Beta(", a,",", b ,") | ", make95CrIFromBeta(a,b)))
      })
      observeEvent(input$sp, {
        a <- input$sp/100*input$sp_n
        b <- input$sp_n-input$sp/100*input$prev_n
        output$sp_dist <- renderText(paste0("sp~Beta(", a,",", b ,") | ", make95CrIFromBeta(a,b)))
      })
      updateNumericInput(inputId="n_sim_outer",value="")
      shinyjs::disable("n_sim_outer")
      output$n_sim_outer_note <- renderText("Outer simulation is not applicable to the selected type of evidence specification on model performance because of conjugate probability distirbutions.")
    }
    if(choice=="2")
    {
      output$evidence_inputs <- renderUI(list(
        fluidRow(
          column(4, sliderInput("prev",label="Expected outcome risk (outcome prevalence) (%)", min=0, max=100, step=0.1, value=c(40,60), width="100%")),
          column(4, textOutput("prev_dist"))
        ),hr(),
        fluidRow(
          column(4, sliderInput("cstat",label="Expected value and upper 95%CI bound of the c-statistic", min=0.51, max=0.99, step=0.001, value=c(0.70,0.80), width="100%")),
          column(4, textOutput("cstat_dist"))
        ),hr(),
        fluidRow(
          column(4, sliderInput("cal_intercept",label="Expected calibration intercept", min=-1, max=1, step=0.01, value=0, width="100%")),
          column(4, numericInput("cal_intercept_sd",label="SE of the intercept", value=0.1)),
          column(4, textOutput("cal_intercept_dist"))
        ),
        fluidRow(
          column(4, sliderInput("cal_slope",label="Expected calibration slope", min=-2, max=2, step=0.01, value=1, width="100%")),
          column(4, numericInput("cal_slope_sd",label="SE of the slope", value=0.1)),
          column(4, textOutput("cal_slope_dist"))
        )
      ))
      observeEvent(input$cstat, {
        res <- evsiexval::solve_beta_given_mean_upper_ci(input$cstat[1],input$cstat[2])
        output$cstat_dist <- renderText(paste0("c~Beta(", round(res$a,2),",",round(res$b,2),") | ", make95CrIFromBeta(res$a,res$b)))
      })
      observeEvent(input$prev, {
        res <- evsiexval::solve_beta_given_mean_upper_ci(input$prev[1]/100,input$prev[2]/100)
        output$prev_dist <- renderText(paste0("prev~Beta(", round(res$a,2),",",round(res$b,2),") | ", make95CrIFromBeta(res$a,res$b)))
      })
      shinyjs::enable("n_sim_outer")
      output$n_sim_outer_note <- renderText("This value cannot be set to >100 for the web app. For larger values use the R package directly")
    }
    if(choice=="3")
    {

    }
  })

  observeEvent(input$prelim_run, {
    evidence <- construct_evidence()

    z <- input$z/100
    N <- input$N
    n_sim_outer <- input$n_sim_outer*1
    n_sim_inner <- input$n_sim_inner*1
    n_stars <- c(0,as.integer(round(exp(seq(from=log(input$n_min), to=log(input$n_max), by=log(input$n_step))))))

    if(evidence$type==1)
    {
      VoI <- EVSI_ag(evidence$params, z, n_sim=1, future_sample_sizes=c())
    }
    if(evidence$type==2)
    {
      samples <<- gen_triplets(n_sim_outer, z=z, prev=evidence$params$prev, cs=evidence$params$cs, A=evidence$params$A, B=evidence$params$B)
      VoI <- EVSI_g(samples[,c('prev','se','sp')], z, n_sim=1, future_sample_sizes=c())
    }

    EVPI <- VoI$EVPI
    p_best <- VoI$p_best

    require("knitr")
    output$prelim_results <- renderUI(list(
      renderText(paste("EVPI=",format(EVPI, nsmall=5))),
      renderText(paste("Population EVPI=",format(EVPI*N, nsmall=2))),
      hr(),
      renderTable(p_best),
      HTML(ifelse(max(p_best$p_best)>0.99, "<B style='color:red; font-weight:bold;'>In more than 99% of simulations the same stratgy had the highest NB. This indicates there is not much uncertainty around this decision. VoI analysis might be degenrate and non-informative.</B>","")),
      actionButton("evsi_run","Run EVSI analysis")
      ))
  })

  observeEvent(input$evsi_run, {
    evidence <- construct_evidence()

    N <- input$N
    z <- input$z/100
    n_sim_outer <- input$n_sim_outer*1
    n_sim_inner <- input$n_sim_inner*1
    n_stars <- c(0,as.integer(round(exp(seq(from=log(input$n_min), to=log(input$n_max), by=log(input$n_step))))))

    if(evidence$type==1)
    {
      VoI <- EVSI_ag(evidence$params, z, n_sim=n_sim_inner, future_sample_sizes=n_stars[-1])
    }
    if(evidence$type==2)
    {
      #samples <- gen_triplets(n_sim_outer, z=z, prev=evidence$params$prev, cs=evidence$params$cs, A=evidence$params$A, B=evidence$params$B)
      VoI <- EVSI_gf(samples[,c('prev','se','sp')], z, n_sim=n_sim_inner, future_sample_sizes=n_stars[-1])
    }

    EVPI <- VoI$EVPI
    EVSIs <- c(0,VoI$EVSI)
    lambda <- input$lambda
    ENBS <- EVSIs*N-n_stars/lambda


    require("knitr")
    output$results <- renderUI(list(
      renderTable(data.frame("sample size"=as.integer(n_stars),
                                                      "EVSI"=format(EVSIs, nsmall=5),
                                                      "Population EVSI"=as.double(EVSIs*N),
                                                      "ENBS"=as.double(ENBS)
                                                      )),
      renderPlot({
        plot(n_stars, EVSIs, type='l', xlab="Sample size of the future study", ylab="EVSI")
        title("Expected Value of Sample Information (in true positive units)")
        y2 <- pretty(c(0,N*EVPI))
        axis(4, at=y2/N, labels=y2)
        mtext("Population EVSI", side = 4)
        lines(c(0,max(n_stars)), rep(EVPI,2), col="gray")
      }),

      renderPlot({
        plot(n_stars, ENBS, type='l', xlab="Sample size of the future study", ylab="ENBS")
        title("Expected Net benefit of Sampling (in true positive units)")
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


shinyApp(ui = ui, server = server)
