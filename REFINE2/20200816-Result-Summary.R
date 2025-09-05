

# Effect_Size <- 6.6
summarise.res <- function(boot1, wrt="all", Effect_Size=6.6){
  tmp <- NULL
  N_boot <- length(boot1)
  for(i in 1:N_boot){
    tmp <- rbind(tmp, boot1[[i]])
  }
  if ("t" %in% colnames(tmp)){
    t.vec <- as.double(as_tibble(tmp)$t)
  }else{
    t.vec <- NA
  }
  sim_corr1 <- as_tibble(tmp) %>% 
    mutate(ATE = as.double(ATE), SE = as.double(SE),
           lb = as.double(ATE) - 1.96*as.double(SE), 
           ub = as.double(ATE) + 1.96*as.double(SE),
           t = t.vec) %>% 
    add_column(iter = c(1:nrow(tmp))) %>% mutate(bias = as.double(ATE) - Effect_Size)
  
  #write_csv(plas_corr1, paste0(out_path, Sys.Date(), "-plasmode_sim_CV-TMLE_GLM_IPW_1000.csv"))
  
  # sim_res <- sim_corr1 %>% group_by(TYPE) %>% 
  #   summarize(b=length(ATE), med_t=median(t), mu_ATE = mean(ATE), med_ATE = median(ATE), 
  #             mu_SE = mean(SE), med_SE = median(SE), 
  #             mu_bias = mean(bias), med_bias = median(bias), 
  #             var = var(ATE), MSE = var + mu_bias^2,
  #             coverage = sum(lb <= Effect_Size & ub >= Effect_Size)/N_boot)
  sim_res <- sim_corr1 %>% group_by(TYPE) %>%
    summarize(
      b = length(ATE),
      mean_ATE = mean(ATE),
      mcse_mean_ATE = sd(ATE)/sqrt(n()),
      
      mean_SE = mean(SE),
      mcse_SE = sd(SE)/sqrt(n()),
      
      mean_bias = mean(bias),
      mcse_bias = sd(bias)/sqrt(n()),
      
      var = var(ATE),
      MSE = var + mean_bias^2,
      mcse_MSE = sd((ATE - Effect_Size)^2) / sqrt(n()),
      
      coverage = mean(lb <= Effect_Size & ub >= Effect_Size),
      mcse_coverage = sqrt(coverage * (1 - coverage) / n())
    )
  if (wrt !="all"){
    sim_res <- sim_res %>%  filter(TYPE %in% wrt)
  }

  # Visualize
  # sim_plot <- sim_corr1 %>% 
  #   ggplot(aes(x = iter, y = ATE, color = TYPE)) +
  #   geom_point() + geom_errorbar(aes(ymin = lb, ymax = ub)) + 
  #   geom_hline(aes(yintercept = Effect_Size)) + 
  #   # geom_hline(data = filter(sim_corr1, substring(TYPE,1,3) %in% c("AIP", "GLM")), aes(yintercept = mean(ATE)), linetype = "dashed") +
  #   # geom_hline(data = filter(sim_corr1, substring(TYPE,1,2) == "CV"), aes(yintercept = mean(ATE)), linetype = "dotted") +
  #   geom_hline(data = sim_corr1, aes(yintercept = mean(ATE)), linetype = "dotted") +
  # labs(x = "iteration", color = "Estimator") +
  #   facet_wrap(~TYPE) +
  #   theme(legend.position = "none")
  # sim_plot
  # return(list(sim_res,sim_plot))
  return(sim_res)
}

summarise.plot <- function(boot1, wrt="all", Effect_Size=6.6){
  tmp <- NULL
  N_boot <- length(boot1)
  for(i in 1:N_boot){
    tmp <- rbind(tmp, boot1[[i]])
  }
  if ("t" %in% colnames(tmp)){
    t.vec <- as.double(as_tibble(tmp)$t)
  }else{
    t.vec <- NA
  }
  sim_corr1 <- as_tibble(tmp) %>% 
    mutate(ATE = as.double(ATE), SE = as.double(SE),
           lb = as.double(ATE) - 1.96*as.double(SE), 
           ub = as.double(ATE) + 1.96*as.double(SE),
           t = t.vec) %>% 
    add_column(iter = c(1:nrow(tmp))) %>% mutate(bias = as.double(ATE) - Effect_Size)
  
  #write_csv(plas_corr1, paste0(out_path, Sys.Date(), "-plasmode_sim_CV-TMLE_GLM_IPW_1000.csv"))
  
  sim_res <- sim_corr1 %>% group_by(TYPE) %>% 
    summarize(b=length(ATE), med_t=median(t), mu_ATE = mean(ATE), med_ATE = median(ATE), 
              mu_SE = mean(SE), med_SE = median(SE), 
              mu_bias = mean(bias), med_bias = median(bias), 
              var = var(ATE), MSE = var + mu_bias^2,
              coverage = sum(lb <= Effect_Size & ub >= Effect_Size)/N_boot)
  if (wrt !="all"){
    sim_res <- sim_res %>%  filter(TYPE %in% wrt)
  }
  
  # Visualize
  sim_plot <- sim_corr1 %>% 
    ggplot(aes(x = iter, y = ATE, color = TYPE)) +
    geom_point() + geom_errorbar(aes(ymin = lb, ymax = ub)) + 
    geom_hline(aes(yintercept = Effect_Size)) + 
    # geom_hline(data = filter(sim_corr1, substring(TYPE,1,3) %in% c("AIP", "GLM")), aes(yintercept = mean(ATE)), linetype = "dashed") +
    # geom_hline(data = filter(sim_corr1, substring(TYPE,1,2) == "CV"), aes(yintercept = mean(ATE)), linetype = "dotted") +
    geom_hline(data = sim_corr1, aes(yintercept = mean(ATE)), linetype = "dotted") +
    labs(x = "iteration", color = "Estimator") +
    ggtitle(paste0("Plasmode Simulations (True ATE set to ", 
                   round(Effect_Size,3),")") )+
    facet_wrap(~TYPE) +
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5,face="bold"))
  # sim_plot
  # return(list(sim_res,sim_plot))
  return(sim_plot)
}

give.summary.res <- function(res.path, out.var,no_time){
  load(res.path)
  if (no_time){
    boot_colnames <- c("ATE", "SE", "TYPE")
  }else{
    boot_colnames <- c("ATE", "SE", "TYPE", "t", "no_cores")
  }
  for (i in 1:length(boot1)){
    colnames(boot1[[i]]) <- boot_colnames
  }
  # print(res.path)
  sum.obj <- summarise.res(boot1,Effect_Size=Effect_Size)[[1]]
  sum.var <- sum.obj[,out.var]
  if (out.var == "med_bias"){
    sum.var <- round(sum.var * 100,2)
  }else if (out.var == "med_SE"){
    sum.var <- round(sum.var* 100,2)
  }else{
    sum.var <- round(sum.var,2)
  }
  return(sum.var)
}

append.res <-function(out.table.str, res.path, out.var,no_time){
  # no_time <- T
  sum.var <- give.summary.res(res.path, out.var,no_time )
  out.table.str<- paste(out.table.str,"&", sum.var)
  return(out.table.str)
}

violin.df.by.situ <- function(path.lst, situ.lst, ci.axis.coeff=0.55,non.par=T,
                              folder ="RDataFiles")
{
  # Notice for DR estimators, we only plot non-parametric
  library(ggplot2)
  library(tidyverse)
  # setwd("~/Desktop/HuangGroup/cvtmle_plasmode/Code")
  # Effect_Size <- 6.6
  # Effect_Size <- 5.51788
  is.timing <- T; tim.str <- ifelse(is.timing,"-timing","")
  
  collapse.add.situ <- function(boot1,situ,par.type ="par"){
    tmp <- NULL
    N_boot <- length(boot1)
    for(i in 1:N_boot){
      colnames(boot1[[i]]) <- c("ATE", "SE", "TYPE", "t", "no_cores")
      tmp <- rbind(tmp, boot1[[i]])
    }
    sum.obj <- summarise.res(boot1,Effect_Size=Effect_Size)[[1]]
    sum.CI <- sum.obj[,"coverage"]$coverage
    tmp <- as_tibble(tmp) %>% mutate(situ = situ, coverage=rep(sum.CI, N_boot))
    if (tmp$TYPE[1] == "GLM + IPW"){
      tmp$TYPE <- rep("IPW", N_boot)
    }
    if (tmp$TYPE[1] == "ManuTMLE"){
      tmp$TYPE <- rep("TMLE", N_boot)
    }
    if (par.type=="non.par"){
      tmp$TYPE <-  paste0(tmp$TYPE, ".NS")
      par.type <- "NS"
    }else{
      if ((tmp$TYPE[1]  != "IPW") && (tmp$TYPE[1]  != "GComp")){
        tmp$TYPE <-  paste0(tmp$TYPE, ".S")
      }
      par.type <- "S"
    }
    tmp$par.type <- rep(par.type, N_boot)
    return(tmp)
  }
  violin.df <- data.frame()
  
  for (i in 1:length(path.lst)){
    #i<-1
    path <- path.lst[i]
    situ <- situ.lst[i]
    for (mtd in mtd.lst){
      # mtd <- "DCTMLE"
      res.path <- paste0(path,folder, "/result-",mtd,"-par",tim.str,".RData")
      load(res.path)
      violin.df <- rbind(violin.df, collapse.add.situ(boot1,situ))
      if (non.par){
        if ((mtd  != "IPW") && (mtd  != "GComp")){ # append non-par result
          res.path <- paste0(path,folder, "/result-",mtd,"-non-par",tim.str,".RData")
          load(res.path)
          violin.df <- rbind(violin.df, collapse.add.situ(boot1,situ,par.type = "non.par"))
        }
      }
    }
  }
  violin.df$TYPE <- factor(violin.df$TYPE, levels =factor.lst)
  # unique(violin.df$TYPE)s
  violin.df$situ <- factor(violin.df$situ, levels=situ.lst)
  violin.df$ATE <- as.double(violin.df$ATE)
  violin.df <- violin.df %>% mutate(bias = ATE - Effect_Size)
  return(violin.df)
}



violin.plot <- function(violin.df, ci.axis.coeff=0.55){
  library(cowplot)
  library(ggplot2)
  ggplot(violin.df, aes(x=TYPE)) +
    geom_violin(aes(y=bias, fill=situ)) +
    scale_y_continuous(
      # Features of the first axis
      name = "Bias",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~.*ci.axis.coeff, name="CI coverage")
    )+
    geom_hline(yintercept=0.95/ci.axis.coeff, linetype="dashed", color = "red")+
    geom_point(aes(y=coverage/ci.axis.coeff, fill=situ)) +facet_wrap(~ situ)+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    scale_fill_discrete(name = "Situation")
    
}

violin.plot.by.situ <- function(path.lst, situ.lst, ci.axis.coeff=0.55,folder ="RDataFiles")
{
  # Notice for DR estimators, we only plot non-parametric
  violin.df <- violin.df.by.situ(path.lst, situ.lst, ci.axis.coeff=0.55,folder =folder)
  violin.plot(violin.df, ci.axis.coeff)
}

violin.and.bar.plot <- function(violin.df, ci.axis.coeff=0.55){
  library(cowplot)
  library(ggplot2)
  # input: violin.df for one situation
  bias.plot <- ggplot(violin.df, aes(x=TYPE,y=bias, fill=par.type)) +
    geom_violin(width=1) +
    geom_boxplot(width=0.1, alpha=0.2, position = position_dodge(width=8.5)) +
    # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    scale_fill_brewer(palette="RdBu",name = "Estimator Class")+
    theme(axis.text=element_text(size=18, angle = 45),
          axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=22,face="bold"))
    # scale_fill_discrete()
    
  covg.df <- violin.df[seq(1,nrow(violin.df),by=100),]
  covg.plot <- ggplot(covg.df, aes(x=TYPE, y=coverage, fill=par.type)) +
    geom_bar(stat="identity",width = 0.5, position = position_dodge(width=8.5)) +
    geom_text(aes(label=coverage), vjust=-0.3, size=6.5)+
    geom_hline(yintercept=0.95, linetype="dashed", color = "red")+
    scale_fill_brewer(palette="RdBu",name = "Estimator Class")+
    theme(axis.text=element_text(size=18, angle = 45),
          axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=20),
          legend.title=element_text(size=22,face="bold"))
  
  # title <- ggdraw() + draw_label(paste0("Bias and Coverage Plot for ",unique(covg.df$situ)[1]), fontface='bold')
  # title <- ggdraw() + draw_label(paste0("Bias and Coverage Plot for
  # bottom_row <- plot_grid(nutrient_boxplot, tss_flow_plot, ncol = 2, labels = "AUTO")
  # plot_grid(title,bias.plot, covg.plot, nrow = 3, labels = c("", "", ""),
  #           rel_heights = c(0.2, 1, 1))
  plot_grid(bias.plot, covg.plot, nrow = 2, labels = c( "", ""),
            rel_heights = c(1, 1))
  
}

violin.and.bar.plot.by.situ <- function(path.lst, situ.lst, ci.axis.coeff=0.55,folder ="RDataFiles")
{
  # Notice for DR estimators, we only plot non-parametric
  violin.df <- violin.df.by.situ(path.lst, situ.lst, ci.axis.coeff=0.55,non.par = T,folder =folder)
  violin.and.bar.plot(violin.df, ci.axis.coeff)
}





