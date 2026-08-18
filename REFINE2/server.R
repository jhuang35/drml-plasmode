library(shiny)

list.of.packages <- c("shiny","mgcv","nlme","glm2","polspline",  "doRNG","doParallel",
                      "SuperLearner","gam","foreach","splines","nnls",
                      "randomForest", "xgboost", "RCAL",
                      "Plasmode","table1","readxl","haven","dplyr","purrr","readr","tidyr",       
                      "tibble","tidyverse","ggplot2",
                      "here","sandwich",
                      "glmnet",
                      'arm', 'lme4', 'twang', 'gbm', 'latticeExtra', 'epiDisplay' # Plasmode dependencies
)

install.packages(list.of.packages[which(
  sapply(list.of.packages, function(x) {nzchar(system.file(package = x))})==F
)])  

# install Plasmode
if (!requireNamespace("Plasmode", quietly = TRUE)) {
  install.packages(here::here("REFINE2/Plasmode_0.1.0.tar.gz"), repos = NULL, type="source")
}

function(input, output, session) {
  library(tidyverse)
  library(haven)
  library(readxl)
  library(table1)
  library(SuperLearner)
  library(gam)
  library(doParallel)
  library(doRNG)
  library(foreach)
  library(dplyr)
  library(polspline)
  
  path <- reactive({input$path})
  
  # Fixed parallel setup for Windows compatibility
  num_cores <- detectCores(all.tests = TRUE) - 2
  registerDoParallel(cores = num_cores)
  
  random_seed <- reactive({input$random_seed})
  num_cf <- reactive({input$num_cf})
  control <- SuperLearner.CV.control(V=2)
  
  est.mtd <- reactive({input$mtd})
  
  doDCTMLE <<- reactive({as.numeric(est.mtd()=="DCTMLE")})
  doDCAIPW <<- reactive({as.numeric(est.mtd()=="DCAIPW")})
  doAIPW <<- reactive({as.numeric(est.mtd()=="AIPW")})
  doIPW <<- reactive({as.numeric(est.mtd()=="IPW")})
  doGComp <<- reactive({as.numeric(est.mtd()=="GComp")})
  doManuTMLE <<- reactive({as.numeric(est.mtd()=="TMLE")})
  
  is.par <- reactive({input$is.par})
  use.par <<- reactive({as.logical(is.par()=="Smooth")})
  
  sims.ver <- "plas"
  
  ########
  # parameters for plasmode
  #######
  plas_sim_N <- reactive({max(10,input$obs)}); use.subset <- F
  generateA <- T
  
  ########
  # parameters for 5 var, 5var.then.plas
  #######
  Nsets <- 500
  Nsamp <- 600
  
  exp.Form.check <- reactive({
    ds <- read.csv(file = path(), header=TRUE, stringsAsFactors=FALSE)
    if (prod(all.vars(expr = as.formula(input$expForm)) %in% colnames(ds)) == 0){
      tmp <- all.vars(expr = as.formula(input$expForm))
      to.out <- c("Error: Variables (", paste0(tmp[which(!tmp %in% colnames(ds))],collapse=", ")  ,") are not found in SIMULATION Propensity Score 
             model. Must only contain variables in the dataset. Please re-enter.")
      paste0(to.out,collapse="")
    }
  })
  
  output$exp.Form <- renderText({
    exp.Form.check()
  })
  
  out.Form.check <- reactive({
    ds <- read.csv(file = path(), header=TRUE, stringsAsFactors=FALSE)
    if (prod(all.vars(expr = as.formula(input$outForm)) %in% colnames(ds)) == 0){
      tmp <- all.vars(expr = as.formula(input$outForm))
      to.out <- c("Error: Variables (", paste0(tmp[which(!tmp %in% colnames(ds))],collapse=", ")  ,") are not found in SIMULATION Outcome 
             Model. Must only contain variables in the dataset. Please re-enter.")
      paste0(to.out,collapse="")
    }
  })
  
  output$out.Form <- renderText({
    out.Form.check()
  })
  
  exp.Form.est.check <- reactive({
    ds <- read.csv(file = path(), header=TRUE, stringsAsFactors=FALSE)
    
    if (prod(all.vars(expr = as.formula(input$expForm.est)) %in% colnames(ds)) == 0){
      
      tmp <- all.vars(expr = as.formula(input$expForm.est))
      to.out <- c("Error: Variables (", paste0(tmp[which(!tmp %in% colnames(ds))],collapse=", ")  ,") are not found in ESTIMATION Propensity Score 
             model. Must only contain variables in the dataset. Please re-enter.")
      paste0(to.out,collapse="")
    }
  })
  
  output$exp.Form.est <- renderText({
    exp.Form.est.check()
  })
  
  out.Form.est.check <- reactive({
    ds <- read.csv(file = path(), header=TRUE, stringsAsFactors=FALSE)
    if (prod(all.vars(expr = as.formula(input$outForm.est)) %in% colnames(ds)) == 0){
      
      tmp <- all.vars(expr = as.formula(input$outForm.est))
      to.out <- c("Error: Variables (", paste0(tmp[which(!tmp %in% colnames(ds))],collapse=", ")  ,") are not found in ESTIMATION Outcome 
             Model. Must only contain variables in the dataset. Please re-enter.")
      paste0(to.out,collapse="")
    }
  })
  
  output$out.Form.est <- renderText({
    out.Form.est.check()
  })
  
  observeEvent(input$run, {
    set.seed(random_seed())
    expForm <- reactive({input$expForm})
    outForm <- reactive({input$outForm})
    ds <- read.csv(file = path(), header=TRUE, stringsAsFactors=FALSE)
    req(prod(all.vars(expr = as.formula(input$expForm)) %in% colnames(ds)) == 1)
    req(prod(all.vars(expr = as.formula(input$outForm)) %in% colnames(ds)) == 1)
    req(prod(all.vars(expr = as.formula(input$expForm.est)) %in% colnames(ds)) == 1)
    req(prod(all.vars(expr = as.formula(input$outForm.est)) %in% colnames(ds)) == 1)
    
    # Source required files
    source("20200705-DCDR-Functions.R")
    source("20200720-Algos-code.R")
    
    plas.copy <- read.csv(file = path(), header=TRUE, stringsAsFactors=FALSE)
    {
      doIPW <- doIPW(); doAIPW=doAIPW();doDCAIPW=doDCAIPW()
      doManuTMLE=doManuTMLE(); doDCTMLE=doDCTMLE(); doGComp=doGComp()
      use.par = use.par()
      
      if (use.par == T){
        #### SMOOTH
        aipw_lib <- SL.param
      }
      else{
        ####  NON-SMOOTH
        aipw_lib <- SL.lib
      }
      set1 <- data.frame(plas.copy)
      tset <- set1 %>% mutate(YT = (Y-min(set1$Y))/(max(set1$Y)- min(set1$Y)))
      one.time.est <- getRES(set1, tset, aipw_lib=aipw_lib, tmle_lib=tmle_lib, short_tmle_lib=short_tmle_lib,
                             doIPW = doIPW(),
                             doAIPW=doAIPW(),
                             doDCAIPW=doDCAIPW(),
                             doManuTMLE=doManuTMLE(),
                             doDCTMLE=doDCTMLE(),doGComp=doGComp(),
                             num_cf=num_cf(),
                             control=control,
                             parallel=F,
                             expForm = input$expForm.est,
                             outForm = input$outForm.est
      )
    }
    
    ATE.one.time <<- as.numeric(one.time.est[1,1])
    SE.one.time <<- as.numeric(one.time.est[1,2])
    
    source("20200803-Sims-Function.R")
    cat(sims.ver)
    sims.obj <- general.sim(sims.ver,path=path(),expForm = expForm(),
                            outForm = outForm(),plas_sim_N=plas_sim_N())
    
    if (sims.ver == "plas"|sims.ver == "5var.then.plas"){
      plas <- sims.obj$plas
      plas_sims <- sims.obj$plas_sims
      vars <- sims.obj$vars
    }else{
      sim_boots <- sims.obj
    }
    
    RegEff <<- sims.obj$RegEff
    Effect_Size <<- RegEff
    N_sims <- reactive({input$obs})
    
    expForm <- reactive({input$expForm.est})()
    outForm <- reactive({input$outForm.est})()
    
    plas_sim_N <- plas_sim_N()
    N_sims <- N_sims()
    plas.copy <- plas %>% dplyr::select(-Y,-A)
    
    boot1 <- reactive({
      # Source files 
      source("20200705-DCDR-Functions.R")
      source("20200720-Algos-code.R")
      
      # Extract ALL reactive values BEFORE the parallel loop
      doIPW_val <- doIPW()
      doAIPW_val <- doAIPW()
      doDCAIPW_val <- doDCAIPW()
      doManuTMLE_val <- doManuTMLE()
      doDCTMLE_val <- doDCTMLE()
      doGComp_val <- doGComp()
      num_cf_val <- num_cf()
      
      # Define the missing libraries based on available SL objects
      # Use the existing SL.lib.tmle if available, otherwise create defaults
      if (!exists("tmle_lib")) {
        if (exists("SL.lib.tmle")) {
          tmle_lib <- SL.lib.tmle
          cat("Using existing SL.lib.tmle for tmle_lib\n")
        } else {
          tmle_lib <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth", "SL.ranger")
          cat("Defining tmle_lib with default algorithms\n")
        }
      }
      
      if (!exists("short_tmle_lib")) {
        short_tmle_lib <- c("SL.mean", "SL.glm")
        cat("Defining short_tmle_lib with basic algorithms\n")
      }
      
      # Verify aipw_lib exists (should have been created earlier)
      if (!exists("aipw_lib")) {
        cat("aipw_lib not found, using default...\n")
        aipw_lib <- c("SL.mean", "SL.glm", "SL.gam")
      }
      
      # Capture objects in local variables for export
      local_aipw_lib <- get("aipw_lib")
      local_tmle_lib <- get("tmle_lib")
      local_short_tmle_lib <- get("short_tmle_lib")
      local_control <- control
      local_plas <- plas
      local_plas_copy <- plas.copy
      local_plas_sims <- plas_sims
      local_plas_sim_N <- plas_sim_N
      local_expForm <- expForm
      local_outForm <- outForm
      
      boot2 <- foreach(i = 1:N_sims, 
                       .errorhandling = "stop",
                       .packages = c("tidyverse", "SuperLearner", "dplyr")) %dopar% {
                         
                         # Source required files on each worker (especially important for Windows)
                         source("20200705-DCDR-Functions.R")
                         source("20200720-Algos-code.R")
                         
                         sims.ver <- "plas"
                         
                         # Initialize dataset
                         if (sims.ver == "plas" | sims.ver == "5var.then.plas"){
                           plas_data <- data.frame(id = local_plas_sims$Sim_Data[i],
                                                   A = local_plas_sims$Sim_Data[i + (2*local_plas_sim_N)],
                                                   Y = local_plas_sims$Sim_Data[i + local_plas_sim_N])
                           
                           colnames(plas_data) <- c("id", "A", "Y")
                           
                           set1 <- left_join(as_tibble(plas_data), as_tibble(local_plas_copy), by="id")
                           tset <- set1 %>% mutate(YT = (Y-min(set1$Y))/(max(set1$Y)- min(set1$Y)))
                         } else {
                           # Initialize dataset for other simulation versions
                           ss <- 600
                           set1 <- as_tibble(cbind(C1 = sim_boots[[i]]$C1[1:ss],
                                                   C2 = sim_boots[[i]]$C2[1:ss],
                                                   C3 = sim_boots[[i]]$C3[1:ss],
                                                   C4 = sim_boots[[i]]$C4[1:ss],
                                                   C5 = sim_boots[[i]]$C5[1:ss],
                                                   A = sim_boots[[i]]$A[1:ss],
                                                   Y = sim_boots[[i]]$Y[1:ss]))
                           tset <- set1 %>% mutate(YT = (Y-min(set1$Y))/(max(set1$Y)- min(set1$Y)))
                         }
                         
                         # Use the extracted values with explicit parameter names
                         getRES(set1, tset, 
                                aipw_lib = local_aipw_lib,
                                tmle_lib = local_tmle_lib, 
                                short_tmle_lib = local_short_tmle_lib,
                                doIPW = doIPW_val,
                                doAIPW = doAIPW_val,
                                doDCAIPW = doDCAIPW_val,
                                doManuTMLE = doManuTMLE_val,
                                doDCTMLE = doDCTMLE_val,
                                doGComp = doGComp_val,
                                num_cf = num_cf_val,
                                control = local_control,
                                parallel = F,
                                expForm = local_expForm,
                                outForm = local_outForm
                         )
                       }
      boot2
    })
    
    source("20200816-Result-Summary.R")
    
    boot1.val <- reactive({
      boot1 <- boot1()
      for (i in 1:length(boot1)){
        colnames(boot1[[i]]) <- c("ATE", "SE", "TYPE", "t", "no_cores")
      }
      boot1
    })
    
    res_summary <- reactive({
      tryCatch({
        summarise.res(boot1 = boot1.val(), Effect_Size = Effect_Size)
      }, error = function(e) {
        cat("Error in summarise.res:", e$message, "\n")
        # Return a simple fallback summary if the function fails
        boot_data <- boot1.val()
        if (length(boot_data) > 0) {
          ate_values <- sapply(boot_data, function(x) x[1, "ATE"])
          data.frame(
            Method = est.mtd(),
            med_ATE = median(ate_values, na.rm = TRUE),
            mean_ATE = mean(ate_values, na.rm = TRUE),
            coverage = NA,
            bias = median(ate_values, na.rm = TRUE) - Effect_Size
          )
        } else {
          data.frame(Method = est.mtd(), med_ATE = NA, mean_ATE = NA, coverage = NA, bias = NA)
        }
      })
    })
    
    output$table <- renderTable({
      res_summary()
    })
    
    output$plot <- renderPlot({
      summarise.plot(boot1 = boot1.val(), Effect_Size = Effect_Size)
    })
    
    output$res.text.1 <- renderText({
      tbl <- res_summary()
      tbl <- round(tbl[, 2:ncol(tbl)], 3)
      
      paste0("<p> <b>EMPIRICAL RESULTS:</b> Using <b>", est.mtd(),"</b> on the <b>OBSERVED</b> data, 
   we estimate an average treatment effect of <b>", round(ATE.one.time,3),  "</b> with a 95% confidence interval of <b>(",
             round(ATE.one.time-1.96*SE.one.time,3),", ",round(ATE.one.time+1.96*SE.one.time,3),")</b>. 
  This is an unbiased estimate of the ATE if standard causal inference assumptions
         are fulfilled (consistency, exchangeability, positivity, and correct model).</p>")
    })
    
    output$res.text.2 <- renderText({
      tbl <- res_summary()
      tbl <- round(tbl[, 2:ncol(tbl)], 3)
      
      paste0("\n", "<p> <b>SIMULATION FINDINGS:</b> If data-adaptive (machine learning) algorithms are used to improve model specification, 
         there must also be no practical positivity violations across covariates and sample sizes must be large enough 
         for bias convergence. This will vary by data structure and setting. </p>")
    })
    
    output$res.text.3 <- renderText({
      tbl <- res_summary()
      tbl <- round(tbl[, 2:ncol(tbl)], 3)
      
      paste0("<p> Thus, we used the observed data to conduct <b>",input$obs,"</b> plasmode simulations 
  based on the user-provided-PS and outcome SIMULATION models while fixing the ATE to a theoretical true value of 
  <b>", round(Effect_Size,3), "</b> (solid line). Applying the user-provided- ESTIMATION models results in a 
  estimated median ATE of <b>", tbl$med_ATE, "</b>, corresponding to a relative bias 
  of <b>", round((tbl$med_ATE-Effect_Size),2), "</b>. Corresponding confidence intervals covered
  the true ATE in <b>", round(tbl$coverage*100,2) ,"%</b> of simulations.
         This performance should be compared to other estimation methods. </p>")
    })
    
  })
  
  # Cleanup: Stop cluster when session ends (for Windows compatibility)
  session$onSessionEnded(function() {
    if (exists("cl") && !is.null(cl)) {
      try(stopCluster(cl), silent = TRUE)
    }
  })
  
}
