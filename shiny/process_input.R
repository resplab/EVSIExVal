process_input <- function(input)
{
  global_vars$str_val <<- ""
  global_vars$str_desc <<- ""

  res <- T #returns true if all Ok, otherwise false. Materials is dumped into global_vars

  if(is.null(global_vars$evidence_type) | global_vars$evidence_type=="")
  {
    global_vars$str_val <<- c(global_vars$str_val, "You have not characterized model performance.")
    return(F)
  }

  global_vars$z <<- input$z/100
  global_vars$str_desc <<- paste(global_vars$str_desc, "Risk threshold (z):", global_vars$z, "\n")

  global_vars$lambda<<-input$lambda
  global_vars$str_desc <<- paste(global_vars$str_desc, "Number needed to study (NNS):", global_vars$lambda, "\n")

  global_vars$N <<- input$N
  global_vars$str_desc <<- paste(global_vars$str_desc, "Population scaling factor (N):", global_vars$N, "\n")

  global_vars$n_sim_outer <<- input$n_sim_outer*1 #Description for this variable will be provided only if evidence type is relevent so it is relegated there

  global_vars$n_sim_inner <<- input$n_sim_inner*1
  global_vars$str_desc <<- paste(global_vars$str_desc, "Number of inner Monte Carlo simulations:", global_vars$n_sim_inner, "\n")

  global_vars$n_stars <<-NULL
  try(
    {global_vars$n_stars <<- c(0,as.integer(round(exp(seq(from=log(input$n_min), to=log(input$n_max), by=log(input$n_step))))))},
  )
  if(is.null(global_vars$n_stars))
  {
    global_vars$str_val <<- paste(global_vars$str_val, "Sample size specification is invalid");
    return(F)
  }
  global_vars$str_desc <<- paste(global_vars$str_desc, "Sample sizes evaluated (n_stars):", paste0(global_vars$n_stars,collapse=","), "\n")



  if(global_vars$evidence_type=="ind_beta")
  {
    global_vars$str_desc <<- paste0(global_vars$str_desc, "Evidence type: independent Beta distributions (",global_vars$evidence_type,").\n")

    if(is.null(global_vars$n_sim_inner) | global_vars$n_sim_inner=="")
    {
      global_vars$str_val <<- paste(global_vars$str_val,"The number of inner simulations need to be specified.")
      return(F)
    }
    if(global_vars$n_sim_inner>10^6)
    {
      global_vars$str_val <<- paste(global_vars$str_val,"Maximum number of inner simulations should be no more than 10^6. For higher numbers use the R package locally.")
      return(F)
    }
    if(global_vars$n_sim_inner<0)
    {
      global_vars$str_val <<- paste(global_vars$str_val,"Invalid value for the number of inner simulations.")
      return(F)
    }
    global_vars$str_desc <<- paste(global_vars$str_desc, "Number of simulations (n_sim_inner):", global_vars$n_sim_inner, "\n")

    evidence <- list(
      prev=c(round(input$prev/100*input$prev_n), round(input$prev_n-input$prev/100*input$prev_n)),
      se=c(round(input$se/100*input$se_n), round(input$se_n-input$se/100*input$se_n)),
      sp=c(round(input$sp/100*input$sp_n), round(input$sp_n-input$sp/100*input$sp_n))
    )

  }
  if(global_vars$evidence_type=="prev_cs_A_B")
  {
    global_vars$str_desc <<- paste0(global_vars$str_desc, "Evidence type: Prevalence, c-statistic, and calibration (",global_vars$evidence_type,").\n")

    if(is.null(global_vars$n_sim_outer) | global_vars$n_sim_outer=="")
    {
      global_vars$str_val <<- paste(global_vars$str_val,"The number of outer simulations need to be specified.")
      global_vars$str_desc <<- paste(global_vars$str_desc, "Number of outer simulations (n_sim_outer):", global_vars$n_sim_outer, "\n")
      return(F)
    }
    if(global_vars$n_sim_outer>1000)
    {
      global_vars$str_val <<- paste(global_vars$str_val,"Number of outer simulations for this type of evidence cannot be more than 1000. For larger values please use the R package locally.")
      return(F)
    }
    if(is.null(global_vars$n_sim_inner) | global_vars$n_sim_inner=="")
    {
      global_vars$str_val <<- paste(global_vars$str_val,"The number of inner simulations need to be specified.")
      return(F)
    }


    if(global_vars$n_sim_inner>1000)
    {
      global_vars$str_val <<- paste(global_vars$str_val,"Number of inner simulations for this type cannot be more than 1000. For larger values please use the R package locally.")
      return(F)
    }

    evidence <- list(
      prev=c(mean=input$prev[1]/100, upper_ci=input$prev[2]/100),
      cs=c(mean=input$cstat[1], upper_ci=input$cstat[2]),
      A=c(mean=input$cal_intercept/1, sd=input$cal_intercept_sd/1),
      B=c(mean=input$cal_slope/1, sd=input$cal_slope_sd/1)
    )
  }
  if(global_vars$evidence_type=="sample")
  {
    global_vars$str_val <<- "This type of evidence synthesis is not yet implemented"
    res <- F
    evidence <- list()
  }

  global_vars$evidence<<-evidence

  res
}

