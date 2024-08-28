create_table <- function(model_object) {
  res_model <- model_object
  
  # turning into dataframe
  df_model  <- map_df(res_model, readRDS) 
  df_model <- df_model %>% filter(subset == "all")
  # creating a key in the same format as the one in correct_surv
  df_model <- df_model %>% 
    mutate(key = paste0(sim,"_",rr))
  
  # fixing survival for everyone who doesn't have a correct survival time
  df_model <- df_model %>% 
    mutate( Value = ifelse(Variable == "LY", 80-Value, Value))
  
  # no longer need key 
  df_model <- df_model %>%
    select(-key)
  
  df_output <- df_model %>% filter(discount == 1) %>% select(-rr)
  # df_output$intervention[df_output$intervention == "SoC"] <- "Targeted"
  df_output1 <- df_output %>% group_by(Variable, Type) %>% summarise(Mean = mean(Value))
  #write.csv(df_output1, "output.csv", row.names = F)
  
  df_delta <- df_model %>% 
    spread(intervention, Value) %>%   
    mutate(Delta = Targeted-SoC) %>% 
    select(-Targeted,-SoC)
  
  df_model <- df_model %>% 
    group_by(Variable, Type, subset, intervention, discount, rr) %>% 
    summarise( mean = mean(Value),
               q025 = quantile(Value, 0.025),
               q975 = quantile(Value, 0.975)) %>% ungroup()
  df_model <-df_model %>% 
    mutate(mean = case_when(Variable == "Cost"~ paste0(currency(mean,
                                                                symbol = "$" ,digits = 0)),
                            Variable != "Cost"~ paste0(round(mean, 2))),
           q025 = case_when(Variable == "Cost"~ paste0(currency(q025, symbol = "$" ,digits = 0)),
                            Variable != "Cost"~ paste0(round(q025, 2))),
           q975 = case_when(Variable == "Cost"~ paste0(currency(q975, symbol = "$" ,digits = 0)),
                            Variable != "Cost"~ paste0(round(q975, 2))                  
           )
    ) %>% 
    mutate( stat = paste0(mean, " (", q025, ";",q975,")"))
  
  df_model1 <- df_model %>% 
    select(-mean, - q025, -q975)
  
  
  df_delta1 <- df_delta %>% rename(Value = Delta) %>% 
    group_by(Variable, Type, subset, discount, rr) %>% 
    summarise( mean = mean(Value),
               q025 = quantile(Value, 0.025),
               q975 = quantile(Value, 0.975)) %>% ungroup() %>% 
    mutate(mean = case_when(Variable == "Cost"~ paste0(currency(mean, symbol = "$" ,digits = 0)),
                            Variable != "Cost"~ paste0(round(mean, 2))),
           q025 = case_when(Variable == "Cost"~ paste0(currency(q025, symbol = "$" ,digits = 0)),
                            Variable != "Cost"~ paste0(round(q025, 2))),
           q975 = case_when(Variable == "Cost"~ paste0(currency(q975, symbol = "$" ,digits = 0)),
                            Variable != "Cost"~ paste0(round(q975, 2))               
                            
           )
    ) %>% 
    mutate( stat = paste0(mean, " (", q025, ";",q975,")")) %>% 
    select(-mean, - q025, -q975)
  
  df_deltacea <- df_delta %>% 
    filter(Type == "Total" | Type == "QALY", subset == "all", discount == 1 ) %>% select(-Type)  %>% spread(Variable, Delta)
  
  df_model2 <- df_model1 %>% 
    bind_rows(mutate( df_delta1, intervention = "delta")) %>% 
    mutate(Variable = factor(Variable, levels = c("LY","QALY","Cost")),
           Type = factor(Type, levels = c("LY","QALY","Total","PLGG","AE","gen","Treatment","Radiation"))
    ) %>% 
    arrange(Variable, Type, discount, rr) %>% 
    spread(intervention, stat)
  
  df_model3 <- df_model2 %>% 
    arrange(subset, rr, Variable, Type, discount )
  
  
  
  tb2 <-df_model3  %>%
    filter(subset == "all") %>% 
    filter( !(Variable == "Cost" & discount == 0),
            !(Variable == "QALY" & discount == 0) ) %>%
    select(-Variable, - subset,-discount  ) %>% 
    rename(Variable = Type,
           Delta = delta,
           intervention = Targeted,
           control = SoC) 
  
  library(pander)
  library(forcats)
  
  # kable(tb2)
  tb2a <- tb2 %>% 
    mutate(Variable = fct_recode(Variable, `Life-years` = "LY", 
                                 QALY = "QALY",
                                 `Total Cost` = "Total")) %>%
    #select(Variable, intervention, control, Delta , rr)
    select(Variable, intervention, control, Delta )
  tb2a$Variable <- as.character(tb2a$Variable)
  tb2a$Variable[tb2a$Variable=="Treatment"] <- "Targeted costs"
  tb2a$Variable[tb2a$Variable=="PLGG"] <- "Non-Targeted PLGG related costs"
  tb2a$Variable[tb2a$Variable=="AE"] <- "Late effect costs"
  tb2a$Variable[tb2a$Variable=="gen"] <- "General population costs"
  tb2a$Variable[tb2a$Variable=="Radiation"] <- "Radiation costs"
  tb2a$Variable <- factor(tb2a$Variable, levels = c("Life-years","QALY","Total Cost","Non-Targeted PLGG related costs",
                                                    "Late effect costs" ,"General population costs", 
                                                    "Targeted costs","Radiation costs"))
  
  # CEAC
  df_deltacea1 <- df_deltacea %>%
    mutate(rr =factor(rr, levels = c(0,1), labels = c("No radiation benefit","Radiation benefit")),
           type = "all")
  df_means <- df_deltacea %>%
    group_by(rr) %>%
    summarise(Cost = mean(Cost),
              QALY = mean(QALY)) %>%
    mutate(type = "Mean Values") %>%
    mutate(rr1 = factor(rr, levels = c(0,1), labels = c("No radiation benefit","Radiation benefit") ))
  toplot_plane <-df_deltacea1 %>%
    ggplot() +
    geom_point(aes(x = QALY, y = Cost), alpha = 0.75 , size = 0.6) +
    # geom_point(aes(x = QALY, y = Cost), color = "black", data = df_means,shape = 2,show.legend = TRUE, size = 3) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    labs(color = "", shape = "Mean",  
         y="Discounted Total Costs",
         x="Discounted QALYs",
         title = "Cost-Effectiveness Acceptability Curve (CEAC)") +
    scale_y_continuous(labels = dollar_format(prefix="$")) +
    #facet_wrap(~ rr) +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  tablee <- kable(tb2a, format="latex", booktabs=TRUE) %>% 
    kable_styling(latex_options=c("scale_down","HOLD_position"))
  return(list(a=tablee,
         b=toplot_plane))
}

create_trace <- function(fileloc) {
  ### TRACE ###
  # Averaging cohort traces across simulation run
  l_prop_df <- list.files(paste0(fileloc, "/trace"),full.names = TRUE)[index]
  prop_df  <- map_df(l_prop_df[1], readRDS) 
  prop_df$cycle <- as.numeric(prop_df$cycle)
  prop_df$intervention <- ifelse(prop_df$intervention == "Targeted", 1, 0)
  prop_df <- apply(prop_df,2,as.numeric)
  prop_df <- as.data.frame(prop_df)
  prop_df <- prop_df %>% dplyr::select("cycle", "pre_prog", "auditory", "cardiovascular", "death_outside", 
                                       "death_plgg", "neurologic", "prog1", "SN", "stroke", "visual",      
                                       "death_cardio",
                                       "prog2", "SN_gen", "death_SN", "death_stroke", "cardiovascular_gen",
                                       "stroke_gen", "intervention")
  # OS: all-cause mortality and PLGG-related mortality
  p_surv <- p_surv_plgg <- as.data.frame(matrix(NA,nrow=nrow(prop_df), ncol=nn))
  p_surv$intervention  <- prop_df$intervention
  p_surv_plgg$intervention <- prop_df$intervention
  p_surv$cycle  <- prop_df$cycle
  p_surv_plgg$cycle <- prop_df$cycle
  p_surv[,1] <- 1 - (prop_df$death_outside + prop_df$death_cardio + prop_df$death_SN + prop_df$death_stroke +
                     prop_df$death_plgg)
  p_surv_plgg[,1] <- 1 - prop_df$death_plgg
  # PFS
  p_pfs <- as.data.frame(matrix(NA,nrow=nrow(prop_df), ncol=nn))
  p_pfs$intervention  <- prop_df$intervention
  p_pfs$cycle  <- prop_df$cycle
  p_pfs[,1] <- 1 -  (prop_df$prog1 + prop_df$prog2 + 
                     prop_df$death_outside + prop_df$death_cardio + prop_df$death_SN + prop_df$death_stroke +
                     prop_df$death_plgg)
  
  if (nn > 1) {
    for (i in 2:nn) {
      prop_df_i  <- map_df(l_prop_df[i], readRDS) 
      prop_df_i$cycle <- as.numeric(prop_df_i$cycle)
      prop_df_i$intervention <- ifelse(prop_df_i$intervention == "Targeted", 1, 0)
      prop_df_i <- prop_df_i %>% dplyr::select("cycle", "pre_prog", "auditory", "cardiovascular", "death_outside", 
                                               "death_plgg", "neurologic", "prog1", "SN", "stroke", "visual", "death_cardio",
                                               "prog2", "SN_gen", "death_SN", "death_stroke", "cardiovascular_gen",
                                               "stroke_gen", "intervention")
      prop_df_i <- apply(prop_df_i,2,as.numeric)
      prop_df_i <- as.data.frame(prop_df_i)
      prop_df <- prop_df + prop_df_i
      # calculate overall survival
      p_death_i <- prop_df_i$death_outside + prop_df_i$death_cardio + prop_df_i$death_SN + prop_df_i$death_stroke +
                   prop_df_i$death_plgg
      p_surv[,i] <- 1 - p_death_i
      # calculate plgg-related survival
      p_death_plgg_i <- prop_df_i$death_plgg
      p_surv_plgg[,i] <- 1 - p_death_plgg_i
      # calculate PFS
      p_outcome_i <- prop_df_i$prog1 + prop_df_i$prog2 + 
                     prop_df_i$death_outside + prop_df_i$death_cardio + prop_df_i$death_SN + prop_df_i$death_stroke +
                     prop_df_i$death_plgg
      p_pfs[,i] <- 1 - p_outcome_i
      
    }
    prop_df <- prop_df/nn
  }

  prop_df$death <- prop_df$death_outside + prop_df$death_cardio + prop_df$death_SN + prop_df$death_stroke
  prop_df$intervention <- ifelse(prop_df$intervention == 1, "Targeted", "SoC")
  
  prop_long <- prop_df %>% select(cycle, pre_prog, prog1, prog2, death_plgg, death, intervention) %>%
    tidyr::pivot_longer(-c(cycle,intervention), names_to = "state", values_to = "proportion") %>%
    mutate(state=case_when(state=="pre_prog"~"Pre-progression",
                           state=="prog1" ~ "1st progression",
                           state=="prog2" ~ "2nd+ progression",
                           state=="death_plgg" ~ "PLGG-related death",
                           state=="death" ~ "All-cause mortality"),
           state=factor(state, levels=c("Pre-progression", "1st progression", "2nd+ progression",
                                        "PLGG-related death", "All-cause mortality")),
           cycle=cycle/12)
  ggplot(data = prop_long, aes(x = cycle, y = proportion, color = state, group = state)) +
    geom_line() +
    labs(title = "Trace Plot of PLGG States Across Cycles",
         x = "Year",
         y = "Proportion",
         color = "State") +
    scale_color_manual(values = c("green","orange","red","navy", "black")) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10)) +
    #geom_hline(yintercept=0.4, linetype="dashed", color = "purple") +
    #geom_vline(xintercept=25, linetype="dashed", color = "purple") +
    theme_bw() +
    theme(legend.position = "bottom", legend.title=element_text(size=10), legend.text=element_text(size=7)) + 
    facet_wrap(~intervention) -> a
  
  ### AE ###
  l_AE_df <- list.files(paste0(fileloc, "/AE_track"),
                        full.names = TRUE)[index]
  AE_df <- map_df(l_AE_df[1], readRDS) 
  ## Targeted
  AE_df_Targeted <- AE_df %>% filter(intervention == "Targeted")
  AE_prop_Targeted <-apply(AE_df_Targeted %>% dplyr::select(-intervention), 2, function(x){prop.table(table(x))})
  df_AE_Targeted <- data.frame(cycle=0:960)
  df_AE_Targeted$cycle <- as.factor(df_AE_Targeted$cycle)
  for (i in 1:length(AE_prop_Targeted)){
    df_AE_prop <- as.data.frame(AE_prop_Targeted[[i]])
    colnames(df_AE_prop) <- c("cycle",names(AE_prop_Targeted)[i])
    df_AE_Targeted <- df_AE_Targeted %>% left_join(df_AE_prop, by = c("cycle"))
  }
  df_AE_Targeted[1,-1] <- 0
  df_AE_Targeted$cycle <- as.numeric(df_AE_Targeted$cycle)
  df_AE_Targeted[is.na(df_AE_Targeted)] <- 0
  if (nn == 1) {
    df_AE_Targeted$intervention <- "Targeted"
  }
  
  if (nn > 1) {
    for (i in 2:nn) {
      AE_df <- map_df(l_AE_df[i], readRDS) 
      AE_df_Targeted_i <- AE_df %>% filter(intervention == "Targeted")
      AE_prop_Targeted_i <-apply(AE_df_Targeted_i %>% dplyr::select(-intervention), 2, function(x){prop.table(table(x))})
      df_AE_Targeted_i <- data.frame(cycle=0:960)
      df_AE_Targeted_i$cycle <- as.factor(df_AE_Targeted_i$cycle)
      for (j in 1:length(AE_prop_Targeted_i)){
        df_AE_prop <- as.data.frame(AE_prop_Targeted_i[[j]])
        colnames(df_AE_prop) <- c("cycle",names(AE_prop_Targeted_i)[j])
        df_AE_Targeted_i <- df_AE_Targeted_i %>% left_join(df_AE_prop, by = c("cycle"))
      }
      df_AE_Targeted_i[1,-1] <- 0
      df_AE_Targeted_i$cycle <- as.numeric(df_AE_Targeted_i$cycle)
      df_AE_Targeted_i[is.na(df_AE_Targeted_i)] <- 0
      df_AE_Targeted <- df_AE_Targeted + df_AE_Targeted_i
    }
    df_AE_Targeted <- df_AE_Targeted/nn
    df_AE_Targeted$intervention <- "Targeted"
  }
  ## SoC
  AE_df_SoC <- AE_df %>% filter(intervention == "SoC")
  AE_prop_SoC <-apply(AE_df_SoC %>% dplyr::select(-intervention), 2, function(x){prop.table(table(x))})
  df_AE_SoC <- data.frame(cycle=0:960)
  df_AE_SoC$cycle <- as.factor(df_AE_SoC$cycle)
  for (i in 1:length(AE_prop_SoC)){
    df_AE_prop <- as.data.frame(AE_prop_SoC[[i]])
    colnames(df_AE_prop) <- c("cycle",names(AE_prop_SoC)[i])
    df_AE_SoC <- df_AE_SoC %>% left_join(df_AE_prop, by = c("cycle"))
  }
  df_AE_SoC[1,-1] <- 0
  df_AE_SoC$cycle <- as.numeric(df_AE_SoC$cycle)
  df_AE_SoC[is.na(df_AE_SoC)] <- 0
  if (nn == 1) {
    df_AE_SoC$intervention <- "SoC"
  }
  
  if (nn > 1) {
    for (i in 2:nn) {
      AE_df <- map_df(l_AE_df[i], readRDS) 
      AE_df_SoC_i <- AE_df %>% filter(intervention == "SoC")
      AE_prop_SoC_i <-apply(AE_df_SoC_i %>% dplyr::select(-intervention), 2, function(x){prop.table(table(x))})
      df_AE_SoC_i <- data.frame(cycle=0:960)
      df_AE_SoC_i$cycle <- as.factor(df_AE_SoC_i$cycle)
      for (j in 1:length(AE_prop_SoC_i)){
        df_AE_prop <- as.data.frame(AE_prop_SoC_i[[j]])
        colnames(df_AE_prop) <- c("cycle",names(AE_prop_SoC_i)[j])
        df_AE_SoC_i <- df_AE_SoC_i %>% left_join(df_AE_prop, by = c("cycle"))
      }
      df_AE_SoC_i[1,-1] <- 0
      df_AE_SoC_i$cycle <- as.numeric(df_AE_SoC_i$cycle)
      df_AE_SoC_i[is.na(df_AE_SoC_i)] <- 0
      df_AE_SoC <- df_AE_SoC + df_AE_SoC_i
    }
    df_AE_SoC <- df_AE_SoC/nn
    df_AE_SoC$intervention <- "SoC"
  }
  
  df_AE_Targeted0 <- as.data.frame(apply(df_AE_Targeted, 2, cumsum))
  df_AE_Targeted0$cycle <- df_AE_Targeted$cycle
  df_AE_Targeted0$intervention <- df_AE_Targeted$intervention
  df_AE_SoC0 <- as.data.frame(apply(df_AE_SoC, 2, cumsum))
  df_AE_SoC0$cycle <- df_AE_SoC$cycle
  df_AE_SoC0$intervention <- df_AE_SoC$intervention
  df_AE <- rbind(df_AE_Targeted0, df_AE_SoC0)
  
  prop_long <- df_AE %>% select(cycle,auditory,cardiovascular,neurologic,stroke,visual,intervention) %>%
    tidyr::pivot_longer(-c(cycle,intervention), names_to = "state", values_to = "proportion") %>%
    mutate(cycle=cycle/12)
  prop_long <- prop_long %>% mutate(intervention = ifelse(intervention == "SoC", "SoC = \nStandard Chemotherapy", "Targeted = \nDabrafenib-Trametinib"))
  ggplot(data = prop_long, aes(x = cycle, y = proportion, color = state, group = state)) +
    geom_line() +
    labs(title = "Cumulative Incidence of PLGG-related Late Effects Across Cycles",
         x = "Time from initiation of systemic therapy (Years)",
         y = "Proportion",
         color = "PLGG-related late effects") +
    theme_bw() +
    scale_y_continuous(limits=c(0,0.15)) + 
    theme(legend.position = "bottom", legend.text=element_text(size=6.5), legend.title=element_text(size=8.25),
          plot.title = element_text(size=11),axis.title.x=element_text(size=9.5), axis.title.y=element_text(size=9.5)) +
    facet_wrap(~intervention) -> b
  
  prop_long <- df_AE %>% select(cycle,SN,SN_gen,cardiovascular_gen,death_outside,death_SN,stroke_gen,death_cardio,death_stroke, 
                                intervention) %>%
    rename(`stroke general` = stroke_gen,
           `cardiovascular general` = cardiovascular_gen,
           `SN general` = SN_gen,
           `death outside` = death_outside,
           `death stroke` = death_stroke,
           `death SN` = death_SN,
           `death cardio` = death_cardio) %>% 
    tidyr::pivot_longer(-c(cycle,intervention), names_to = "state", values_to = "proportion") %>%
    mutate(cycle=cycle/12) 
  ggplot(data = prop_long, aes(x = cycle, y = proportion, color = state, group = state)) +
    geom_line() +
    labs(title = "Cumulative Incidence of Background Late Effects Across Cycles",
         x = "Time from initiation of systemic therapy (Years)",
         y = "Proportion",
         color = "Background late effects") +
    theme_bw() +
    scale_y_continuous(limits=c(0,0.3)) + 
    theme(legend.position = "right", legend.text=element_text(size=6.5), legend.title=element_text(size=8.25),
          plot.title = element_text(size=11), axis.title.x=element_text(size=9.5), axis.title.y=element_text(size=9.5)) +
    facet_wrap(~intervention) -> c
  
  ### SURVIVAL ###
  ## Overall
  p_surv$intervention <- ifelse(p_surv$intervention == 1, "Targeted = Dabrafenib-Trametinib", "SoC = Standard Chemotherapy")
  p_surv1 <- data.frame(avg = apply(p_surv %>% select(-cycle, -intervention), 1, mean),
                        lb = apply(p_surv %>% select(-cycle, -intervention), 1, function(x){quantile(x, 0.025)}),
                        ub = apply(p_surv %>% select(-cycle, -intervention), 1, function(x){quantile(x, 0.975)}))
  p_surv1$intervention <- p_surv$intervention
  p_surv1$cycle <- p_surv$cycle/12
  ggplot(p_surv1, aes(x = cycle)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = intervention), alpha = 0.15) +  
    geom_line(aes(y = avg, color = intervention, linetype = intervention)) +  
    labs(x = "Time from initiation of systemic therapy (Years)", y = "Survival Probability", 
         title = "Overall Survival with 95% CI") +
    theme_bw() + #facet_wrap(~intervention) + 
    theme(plot.title = element_text(size=9),
          axis.title.x=element_text(size=8), axis.title.y=element_text(size=8), 
          axis.text.y=element_text(size=6.2)) -> d
  
  ## PLGG-related
  p_surv_plgg$intervention <- ifelse(p_surv_plgg$intervention == 1, "Targeted = Dabrafenib-Trametinib", "SoC = Standard Chemotherapy")
  p_surv_plgg1 <- data.frame(avg = apply(p_surv_plgg %>% select(-cycle, -intervention), 1, mean),
                             lb = apply(p_surv_plgg %>% 
                                          select(-cycle, -intervention), 1, function(x){quantile(x, 0.025)}),
                             ub = apply(p_surv_plgg %>% 
                                          select(-cycle, -intervention), 1, function(x){quantile(x, 0.975)}))
  p_surv_plgg1$intervention <- p_surv_plgg$intervention
  p_surv_plgg1$cycle <- p_surv_plgg$cycle/12
  ggplot(p_surv_plgg1, aes(x = cycle)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = intervention), alpha = 0.15) +  
    geom_line(aes(y = avg, color = intervention, linetype = intervention)) +  
    labs(x = "Time from initiation of systemic therapy (Years)", y = "Survival Probability", 
         title = "PLGG-death-related Survival with 95% CI") +
    theme_bw() + #facet_wrap(~intervention) +  
    theme(plot.title = element_text(size=9), 
          axis.title.x=element_text(size=8), axis.title.y=element_text(size=8),
          axis.text.y=element_text(size=6.2))  -> e
  
  ## PFS
  p_pfs$intervention <- ifelse(p_pfs$intervention == 1, "Targeted = Dabrafenib-Trametinib", "SoC = Standard Chemotherapy")
  p_pfs1 <- data.frame(avg = apply(p_pfs %>% select(-cycle, -intervention), 1, mean),
                             lb = apply(p_pfs %>% 
                                          select(-cycle, -intervention), 1, function(x){quantile(x, 0.025)}),
                             ub = apply(p_pfs %>% 
                                          select(-cycle, -intervention), 1, function(x){quantile(x, 0.975)}))
  p_pfs1$intervention <- p_pfs$intervention
  p_pfs1$cycle <- p_pfs$cycle/12
  ggplot(p_pfs1, aes(x = cycle)) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = intervention), alpha = 0.15) +  
    geom_line(aes(y = avg, color = intervention, linetype = intervention)) +  
    labs(x = "Time from initiation of systemic therapy (Years)", y = "Survival Probability", 
         title = "Progression-Free Survival with 95% CI") +
    theme_bw() + #facet_wrap(~intervention) +  
    theme(plot.title = element_text(size=9), 
          axis.title.x=element_text(size=8), axis.title.y=element_text(size=8),
          axis.text.y=element_text(size=6.2))  -> f
  
  return(list(trace_plot=a, ae_plot1 = b, ae_plot2 = c, surv_plot1 = d, surv_plot2 = e, surv_plot3 = f))
}
