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

source("process_input.R")

make95CrIFromBeta <- function(a,b)
{
  paste0("95%CI ",
         round(qbeta(0.025,a,b)*100,1),
         "%-",
         round(qbeta(0.975,a,b)*100,1),"%"
  )
}


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
      tabPanel("Evidence on model performance",
        HTML("Here we solicit your assessment of the model performance and your uncertainty around this asessment. Your knowledge of model performance can be specified in different ways."),
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
        actionButton("run1", "Run the analysis"), actionButton("clear_results", "Clear results"),
        uiOutput("results1")
      )
    ),
)





















# Define server logic required to draw a histogram
server <- function(input, output)
{
  global_vars <<- list(
    evidence_type="", #ind_beta, prev_cs_A_B, sample
    evidence=list(),
    default_values=list(
      exchange_rates=list(z=0.02, lambda=100),
      evidence=list(
        ind_beta=list(prev = c(43L, 457L), se = c(41L, 2L), sp = c(147, 310)),
        prev_cs_A_B=list(prev=c(point=0.086000, upper_ci=0.110659), cs=c(point=0.8474887, upper_ci=0.9147262), A=c(point=0.5379061, sd=0.3824605), B=c(point=1.1793535, sd=0.1692672))
      )
    ),
    result_level=0,
    sample=data.frame(),
    str_val="", #Will contain text generated by process_input()
    str_desc ="" #General textual description of the parameters
  )

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
    global_vars$result_level<<-0
    output$results1=renderUI(HTML(""))
  })

  observeEvent(input$evidence_type, {
    choice <- substring(input$evidence_type,1,1)
    output$n_sim_outer_note <- renderText("")
    output$evidence_inputs <- renderUI("")
    if(choice=="0")
    {
      global_vars$evidence_type <<- ""
    }
    if(choice=="1")
    {
      global_vars$evidence_type <<- "ind_beta"
      output$evidence_inputs <- renderUI(list(
        renderText("This method is based on trivariate specification (outcome prevalence, sensitivity, and specificity), with Beta distribution for modeling uncertainty for each component."),
        hr(),
        fluidRow(
          column(4,sliderInput("prev",label="Expected outcome risk (outcome prevalence) (%)", min=0, max=100, step=0.1, value=global_vars$default_values$evidence$ind_beta$prev[1]/sum(global_vars$default_values$evidence$ind_beta$prev)*100, width="100%")),
          column(4,numericInput("prev_n",label="Your estimate of the average risk of the clinical outcome is based on how many observations?", value=sum(global_vars$default_values$evidence$ind_beta$prev))),
          column(4, textOutput("prev_dist"))
        ),hr(),
        fluidRow(
          column(4,sliderInput("se",label="Sensitivity of the model at the chosen risk threshold (%)", min=0, max=100, step=0.1, value=global_vars$default_values$evidence$ind_beta$se[1]/sum(global_vars$default_values$evidence$ind_beta$se)*100, width="100%")),
          column(4,numericInput("se_n",label="Your estimate of sensitivity is based on how many observations?", value=sum(global_vars$default_values$evidence$ind_beta$se))),
          column(4, textOutput("se_dist"))
        ),hr(),
        fluidRow(
          column(4,sliderInput("sp",label="Specificity the model at the chosen risk threshold (%)", min=0, max=100, step=0.1, value=global_vars$default_values$evidence$ind_beta$sp[1]/sum(global_vars$default_values$evidence$ind_beta$sp)*100, width="100%")),
          column(4,numericInput("sp_n",label="Your estimate of sensitivity is based on how many observations?", value=sum(global_vars$default_values$evidence$ind_beta$sp))),
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
      shinyjs::disable("n_sim_outer")
      output$n_sim_outer_note <- renderText("Outer simulation is not applicable to the selected type of evidence specification on model performance because of conjugate probability distirbutions.")
    }
    if(choice=="2")
    {
      global_vars$evidence_type <<- "prev_cs_A_B"
      output$evidence_inputs <- renderUI(list(
        renderText("This method is based on the specification of your informaiton on typical performance metrics of a model (c-statistic, and calibration intercept and slope)."),
        hr(),
        fluidRow(
          column(4, sliderInput("prev",label="Expected outcome risk (outcome prevalence) (%)", min=0, max=100, step=0.1, value=global_vars$default_values$evidence$prev_cs_A_B$prev*100, width="100%")),
          column(4, textOutput("prev_dist"))
        ),hr(),
        fluidRow(
          column(4, sliderInput("cstat",label="Expected value and upper 95%CI bound of the c-statistic", min=0.51, max=0.99, step=0.001, value=global_vars$default_values$evidence$prev_cs_A_B$cs, width="100%")),
          column(4, textOutput("cstat_dist"))
        ),hr(),
        fluidRow(
          column(4, sliderInput("cal_intercept",label="Expected calibration intercept", min=-1, max=1, step=0.01, value=global_vars$default_values$evidence$prev_cs_A_B$A[1], width="100%")),
          column(4, numericInput("cal_intercept_sd",label="SE of the intercept", value=global_vars$default_values$evidence$prev_cs_A_B$A[2])),
          column(4, textOutput("cal_intercept_dist"))
        ),
        fluidRow(
          column(4, sliderInput("cal_slope",label="Expected calibration slope", min=-2, max=2, step=0.01, value=global_vars$default_values$evidence$prev_cs_A_B$B[1], width="100%")),
          column(4, numericInput("cal_slope_sd",label="SE of the slope", value=global_vars$default_values$evidence$prev_cs_A_B$B[2])),
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
      global_vars$evidence_type <<- "sample"
      output$evidence_inputs <- renderUI(list(
        renderText("To be implemented. You will be able to upload a CSV file that contains triplets of prevalence(prev), sensitivity (se), and specificity(sp)."),
        hr()
      ))
    }
  })


  observeEvent(input$run1, {
    res <- process_input(input)

    stuff_to_render <- list(renderUI(pre(global_vars$str_desc)))

    if(res==T)
    {
      if(global_vars$evidence_type=="ind_beta")
      {
        VoI <- EVSI_ag(global_vars$evidence, global_vars$z, n_sim=global_vars$n_sim_inner, future_sample_sizes=c())
      }

      if(global_vars$evidence_type=="prev_cs_A_B")
      {
        global_vars$sample <<- gen_triplets(global_vars$n_sim_outer, z=global_vars$z, prev=global_vars$evidence$prev, cs=global_vars$evidence$cs, A=global_vars$evidence$A, B=global_vars$evidence$B)
        stuff_to_render <- c(stuff_to_render, renderText(paste("Results are based on a sample of", nrow(global_vars$sample), "draws.")))
        VoI <- evsiexval::EVSI_g(global_vars$sample[,c('prev','se','sp')], global_vars$z, n_sim=1, future_sample_sizes=c())
      }

      EVPI <- VoI$EVPI
      summary <- VoI$summary

      require("knitr")
      stuff_to_render <- c(stuff_to_render,list(
        renderTable(summary, rownames=T, digits=6),
        HTML(ifelse(max(summary['P_best',])>0.99, "<B style='color:red; font-weight:bold;'>In more than 99% of simulations the same stratgy had the highest NB. This indicates there is not much uncertainty around this decision. VoI analysis might be degenrate and non-informative.</B>","")),
        hr(),
        renderText(paste("EVPI=",format(EVPI, nsmall=5))),
        renderText(paste("Population EVPI=",format(EVPI*global_vars$N, nsmall=2))),
        hr(),
        actionButton("evsi_run","Run EVSI analysis"),
        uiOutput("results2")
        ))
      global_vars$result_level <<- 1
    }
    else
    {
      global_vars$result_level <<- 0
      stuff_to_render <- c(stuff_to_render, renderText(paste("Error:", global_vars$str_val)))
    }
    output$results1 <- renderUI(stuff_to_render)
    output$results2 <- renderUI("")
  })



  observeEvent(input$evsi_run, {
    if(global_vars$evidence_type=="ind_beta")
    {
      VoI <- EVSI_ag(global_vars$evidence, global_vars$z, n_sim=global_vars$n_sim_inner, future_sample_sizes=global_vars$n_stars[-1])
    }
    if(global_vars$evidence_type=="prev_cs_A_B")
    {
      withProgress({
        VoI <- EVSI_gf(global_vars$sample[,c('prev','se','sp')], global_vars$z, n_sim=global_vars$n_sim_inner, future_sample_sizes=global_vars$n_stars[-1])
      }, message="Computing")
    }

    EVPI <- VoI$EVPI
    EVSIs <- c(0,VoI$EVSI)
    ENBS <- EVSIs*global_vars$N-global_vars$n_stars/global_vars$lambda

    require("knitr")
    output$results2 <- renderUI(list(
      renderTable(data.frame("sample size"=as.integer(global_vars$n_stars),
                                                      "EVSI"=format(EVSIs, nsmall=5),
                                                      "Population EVSI"=as.double(EVSIs*global_vars$N),
                                                      "ENBS"=as.double(ENBS)
                                                      )),
      renderPlot({
        plot(global_vars$n_stars, EVSIs, type='l', xlab="Sample size of the future study", ylab="EVSI")
        title("Expected Value of Sample Information (in true positive units)")
        y2 <- pretty(c(0,global_vars$N*EVPI))
        axis(4, at=y2/global_vars$N, labels=y2)
        mtext("Population EVSI", side = 4)
        lines(c(0,max(global_vars$n_stars)), rep(EVPI,2), col="gray")
      }),

      renderPlot({
        plot(global_vars$n_stars, ENBS, type='l', xlab="Sample size of the future study", ylab="ENBS")
        title("Expected Net benefit of Sampling (in true positive units)")
        winner <- which.max(ENBS)
        if(winner!=1 & winner!=length(ENBS))
        {
          lines(c(global_vars$n_stars[winner], global_vars$n_stars[winner]), c(0, max(ENBS)), col='red')
          text(global_vars$n_stars[winner]*1.1,max(ENBS)/2,paste0("Optimal sample size:",global_vars$n_stars[winner]), col='red')
        }
        else
        {
          text(0,max(ENBS)/2,"Edge case!", col='red')
        }
      })
    ))
    global_vars$result_level <<- 2
  })
}


shinyApp(ui = ui, server = server)
