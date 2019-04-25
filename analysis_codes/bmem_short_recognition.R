# Required R package
library(lme4)
library(lmerTest)
library(ggplot2)

# Clear workspace
rm(list=ls())

# Define the local path where the data can be found
# you will need to set the correct path whre the RData file with the data is located
input_file_path="DEFINE PATH OF DATA HERE"

experiment_name = "bmem_short"

# load data
input_filename=paste(input_file_path, experiment_name, "_recognition.Rdata",sep="")  
load(file=input_filename)

# Statistics

# overall performance
# hit rate
hit_rate_all_participants = with(data=subset(recognition_data,isOld.==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))
hit_rate_mean=mean(hit_rate_all_participants)
hit_rate_sd=sd(hit_rate_all_participants)
print("Overall performance")
print(paste("Mean hit rate: ", round(hit_rate_mean*100,2), "% +- ", round(hit_rate_sd*100,2),"%", sep=""))
# correct rejection rate
cr_all_participants=with(data=subset(recognition_data,isOld.==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))
cr_rate_mean=mean(cr_all_participants)
cr_rate_sd=sd(cr_all_participants)
print(paste("Mean correct rejection rate: ", round(cr_rate_mean*100,2), "% +- ", round(cr_rate_sd*100,2),"%", sep=""))
# d prime
false_alarm_all_participants = 1-cr_all_participants
# d prime cannot be computed when hit rate / false alarm = 1 or 0. Adjust for this:
hit_rate_all_participants[hit_rate_all_participants==1] = 0.99
false_alarm_all_participants[false_alarm_all_participants==0] = 0.01
d_prime_all_participants=qnorm(hit_rate_all_participants) - qnorm(false_alarm_all_participants)
d_prime_mean=mean(d_prime_all_participants)
d_prime_sd=sd(d_prime_all_participants)
print(paste("Mean d prime: ", round(d_prime_mean,3), " +- ", round(d_prime_sd,3), sep=""))
# RT
# RT for hits
RT_hits_all_participants=with(data=subset(recognition_data,isOld.==1 & IsCorrectAnsOld==1), tapply(RT_isOld, subjectID, mean, na.rm=T))
RT_hits_mean=mean(RT_hits_all_participants)
RT_hits_sd=sd(RT_hits_all_participants)
print(paste("Mean RT for hits is: ", round(RT_hits_mean,3), " +- ", round(RT_hits_sd,3), sep=""))
# RT for misses
RT_misses_all_participants=with(data=subset(recognition_data,isOld.==1 & IsCorrectAnsOld==0), tapply(RT_isOld, subjectID, mean, na.rm=T))
RT_misses_mean=mean(RT_misses_all_participants, na.rm=T)
RT_misses_sd=sd(RT_misses_all_participants, na.rm=T)
print(paste("Mean RT for misses is: ", round(RT_misses_mean,3), " +- ", round(RT_misses_sd,3), sep=""))
# RT for correct rejections
RT_crs_all_participants=with(data=subset(recognition_data,isOld.==0 & IsCorrectAnsOld==1), tapply(RT_isOld, subjectID, mean, na.rm=T))
RT_crs_mean=mean(RT_crs_all_participants, na.rm=T)
RT_crs_sd=sd(RT_crs_all_participants, na.rm=T)
print(paste("Mean RT for correct rejections is: ", round(RT_crs_mean,3), " +- ", round(RT_crs_sd,3), sep=""))
# RT for false alarms
RT_fas_all_participants=with(data=subset(recognition_data,isOld.==0 & IsCorrectAnsOld==0), tapply(RT_isOld, subjectID, mean, na.rm=T))
RT_fas_mean=mean(RT_fas_all_participants, na.rm=T)
RT_fas_sd=sd(RT_fas_all_participants, na.rm=T)
print(paste("Mean RT for false alarms is: ", round(RT_fas_mean,3), " +- ", round(RT_fas_sd,3), sep=""))

# Compute means
means_isold_accuracy=c(mean(with(data=subset(recognition_probe_items,isGo.==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all go items
                        mean(with(data=subset(recognition_probe_items,isGo.==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all nogo items
                        mean(with(data=subset(recognition_probe_items,PairType=="HH"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all HH items  
                        mean(with(data=subset(recognition_probe_items,PairType=="HM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all HM items
                        mean(with(data=subset(recognition_probe_items,PairType=="LM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all LM items
                        mean(with(data=subset(recognition_probe_items,PairType=="LL"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all LL items
                        mean(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="HH"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HH Go items
                        mean(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="HM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HM Go items
                        mean(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="LM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LM Go items
                        mean(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="LL"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LL Go items
                        mean(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="HH"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HH NoGo items
                        mean(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="HM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HM NoGo items
                        mean(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="LM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LM NoGo items
                        mean(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="LL"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LL NoGo items
                        mean(with(data=recognition_probe_items, tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))) # all items
  )
  means_isold_accuracy=round(means_isold_accuracy*100,2)
  gonogo_category=c('Go','NoGo','All','All','All','All','Go','Go','Go','Go','NoGo','NoGo','NoGo','NoGo','All')
  value_category=c('All','All','HH','HM','LM','LL','HH','HM','LM','LL','HH','HM','LM','LL','All')
  accuracy_all_means=data.frame(row.names = c('HH','HM','LM','LL','All'))
  accuracy_all_means$Go=means_isold_accuracy[c(7:10,1)]
  accuracy_all_means$NoGo=means_isold_accuracy[c(11:14,2)]
  accuracy_all_means$All=means_isold_accuracy[c(3:6,15)]
  
  # compute SD values
  sd_isold_accuracy=c(sd(with(data=subset(recognition_probe_items,isGo.==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all go items
                         sd(with(data=subset(recognition_probe_items,isGo.==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all nogo items
                         sd(with(data=subset(recognition_probe_items,PairType=="HH"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all HH items  
                         sd(with(data=subset(recognition_probe_items,PairType=="HM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all HM items
                         sd(with(data=subset(recognition_probe_items,PairType=="LM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all LM items
                         sd(with(data=subset(recognition_probe_items,PairType=="LL"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all LL items
                         sd(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="HH"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HH Go items
                         sd(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="HM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HM Go items
                         sd(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="LM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LM Go items
                         sd(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="LL"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LL Go items
                         sd(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="HH"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HH NoGo items
                         sd(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="HM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HM NoGo items
                         sd(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="LM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LM NoGo items
                         sd(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="LL"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LL NoGo items
                         sd(with(data=recognition_probe_items, tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))) # all items
  )
  sd_isold_accuracy=round(sd_isold_accuracy*100,2)
  accuracy_all_sd=data.frame(row.names = c('HH','HM','LM','LL','All'))
  accuracy_all_sd$Go=sd_isold_accuracy[c(7:10,1)]
  accuracy_all_sd$NoGo=sd_isold_accuracy[c(11:14,2)]
  accuracy_all_sd$All=sd_isold_accuracy[c(3:6,15)]

  # Compute SE values
  se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }
  # compute SD values
  se_isold_accuracy=c(se(with(data=subset(recognition_probe_items,isGo.==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all go items
                      se(with(data=subset(recognition_probe_items,isGo.==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all nogo items
                      se(with(data=subset(recognition_probe_items,PairType=="HH"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all HH items  
                      se(with(data=subset(recognition_probe_items,PairType=="HM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all HM items
                      se(with(data=subset(recognition_probe_items,PairType=="LM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all LM items
                      se(with(data=subset(recognition_probe_items,PairType=="LL"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all LL items
                      se(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="HH"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HH Go items
                      se(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="HM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HM Go items
                      se(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="LM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LM Go items
                      se(with(data=subset(recognition_probe_items,isGo.==1 & PairType=="LL"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LL Go items
                      se(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="HH"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HH NoGo items
                      se(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="HM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HM NoGo items
                      se(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="LM"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LM NoGo items
                      se(with(data=subset(recognition_probe_items,isGo.==0 & PairType=="LL"), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LL NoGo items
                      se(with(data=recognition_probe_items, tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))) # all items
  )
  se_isold_accuracy=round(se_isold_accuracy*100,2)
  accuracy_all_se=data.frame(row.names = c('HH','HM','LM','LL','All'))
  accuracy_all_se$Go=se_isold_accuracy[c(7:10,1)]
  accuracy_all_se$NoGo=se_isold_accuracy[c(11:14,2)]
  accuracy_all_se$All=se_isold_accuracy[c(3:6,15)]
  
  # combine decriptive results
  descriptive_results_accuracy=accuracy_all_means
  descriptive_results_accuracy$Go=paste(accuracy_all_means$Go,"% (", accuracy_all_sd$Go, "%)",sep="")
  descriptive_results_accuracy$NoGo=paste(accuracy_all_means$NoGo,"% (", accuracy_all_sd$NoGo, "%)",sep="")
  descriptive_results_accuracy$All=paste(accuracy_all_means$All,"% (", accuracy_all_sd$All, "%)",sep="")
  
  # print descriptive results
  print("Accuracy old/new - mean percent of correct answers (standard deviation)")
  print("Only for items from the Go vs. NoGo probe comparisons")
  print(descriptive_results_accuracy)
  
  # Logistic regression analysis
  model_accuracy_by_item_type=glmer(IsCorrectAnsOld ~ 1 + isGo.+ (1|subjectID),data=recognition_probe_items,na.action=na.omit,family=binomial) 
  model_accuracy_by_item_type_and_value=glmer(IsCorrectAnsOld ~ 1 + isGo.+ PairType + (1|subjectID),data=recognition_probe_items,na.action=na.omit,family=binomial) 
  model_accuracy_by_item_type_int_value=glmer(IsCorrectAnsOld ~ 1 + isGo.* PairType + (1|subjectID),data=recognition_probe_items,na.action=na.omit,family=binomial) 
  accuracy_main_effect_value=anova(model_accuracy_by_item_type, model_accuracy_by_item_type_and_value) # check the main effect of value category
  accuracy_interaction_effect=anova(model_accuracy_by_item_type_and_value, model_accuracy_by_item_type_int_value) # check effect of interaction
  summary_model_accuracy_by_item_type_and_value=summary(model_accuracy_by_item_type_and_value)
  pval_accuracy_isgo=summary_model_accuracy_by_item_type_and_value$coefficients[2, "Pr(>|z|)"]/2
  odds_ratio_accuracy_isgo=exp(summary_model_accuracy_by_item_type_and_value$coefficients[2, "Estimate"])
  pval_accuracy_value_category=accuracy_main_effect_value$`Pr(>Chisq)`[2]/2
  pval_accuracy_interaction=accuracy_interaction_effect$`Pr(>Chisq)`[2]
  
  # print logistic regression analysis results
  print("results for accuracy")
  print("all Go vs. NoGo")
  print('one-sided p value:')
  print(pval_accuracy_isgo)
  print("value category (HH / HM / LM / LL)")
  print('one p value:')
  print(pval_accuracy_value_category)
  print("Interaction go/nogo and value category")
  print('two-sided p value:')
  print(pval_accuracy_interaction)
  
  # modelling HH vs. rest
  recognition_probe_items$is_hh=recognition_probe_items$PairType=="HH"
  # Logistic regression analysis
  accuracy_isold_results_all_items=summary(glmer(IsCorrectAnsOld ~ 1 + isGo. + is_hh + (1|subjectID),data=recognition_probe_items,na.action=na.omit,family=binomial)) 
  accuracy_isold_results_all_items_with_interaction=summary(glmer(IsCorrectAnsOld ~ 1 + isGo.*is_hh + (1|subjectID),data=recognition_probe_items,na.action=na.omit,family=binomial)) 
  accuracy_isold_results_HH_items=summary(glmer(IsCorrectAnsOld ~ 1 + isGo. + (1|subjectID),data=subset(recognition_probe_items,(recognition_probe_items$is_hh)),na.action=na.omit,family=binomial))
  accuracy_isold_results_rest_items=summary(glmer(IsCorrectAnsOld ~ 1 + isGo. + (1|subjectID),data=subset(recognition_probe_items,(!recognition_probe_items$is_hh)),na.action=na.omit,family=binomial))
  
  # organize statistics in table
  results_isold_accuracy=rbind(accuracy_isold_results_all_items$coefficients[2,],accuracy_isold_results_all_items$coefficients[3,],accuracy_isold_results_all_items_with_interaction$coefficients[4,],accuracy_isold_results_HH_items$coefficients[2,],accuracy_isold_results_rest_items$coefficients[2,])
  results_isold_accuracy=as.data.frame(results_isold_accuracy)
  colnames(results_isold_accuracy)[4] = "two-sided p"
  results_isold_accuracy$one_sided_p=results_isold_accuracy$`two-sided p`/2
  results_isold_accuracy$odds_ratio=round(exp(results_isold_accuracy$Estimate),3)
  CI_min=exp(results_isold_accuracy$Estimate-1.96*results_isold_accuracy$`Std. Error`)
  CI_max=exp(results_isold_accuracy$Estimate+1.96*results_isold_accuracy$`Std. Error`)
  results_isold_accuracy$CI=paste(round(CI_min,3)," - ", round(CI_max,3), sep="")
  results_isold_accuracy$category=c('all Go vs. NoGo', 'HH vs. the rest', 'all interaction Go/NoGo and HH/rest', 'HH Go vs. NoGo', 'rest Go vs. NoGo')
  
  # create dataframe for plot
  df_accuracy=data.frame(row.names = 1:8)
  df_accuracy$value_level=c("High","High","Medium-high","Medium-high","Medium-low","Medium-low", "Low","Low")
  df_accuracy$item_type=c("Go items","NoGo items","Go items","NoGo items","Go items","NoGo items","Go items","NoGo items")
  df_accuracy$means=c(accuracy_all_means$Go[1],accuracy_all_means$NoGo[1],accuracy_all_means$Go[2],accuracy_all_means$NoGo[2],accuracy_all_means$Go[3],accuracy_all_means$NoGo[3],accuracy_all_means$Go[4],accuracy_all_means$NoGo[4])
  df_accuracy$sd=c(accuracy_all_sd$Go[1],accuracy_all_sd$NoGo[1],accuracy_all_sd$Go[2],accuracy_all_sd$NoGo[2],accuracy_all_sd$Go[3],accuracy_all_sd$NoGo[3],accuracy_all_sd$Go[4],accuracy_all_sd$NoGo[4])
  df_accuracy$se=c(accuracy_all_se$Go[1],accuracy_all_se$NoGo[1],accuracy_all_se$Go[2],accuracy_all_se$NoGo[2],accuracy_all_se$Go[3],accuracy_all_se$NoGo[3],accuracy_all_se$Go[4],accuracy_all_se$NoGo[4])
  
  df_accuracy$value_level = factor(df_accuracy$value_level, levels = df_accuracy$value_level) # for the value levels to be shown on the graph in the correct order (and not sorted by ABC)
  
  # plot the results
  plot_accuracy= ggplot(data=df_accuracy, aes(x=item_type, y=means, fill=value_level)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8) +
    theme_bw() + # white background
    theme(legend.position="top",legend.title=element_blank()) + # position legend
    theme(axis.title.x=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank()) + # axis and background formating
    theme(aspect.ratio=1.2, text = element_text(size=24)) + # font size
    geom_errorbar(aes(ymin=means-se, ymax=means+se), width=1/6, position=position_dodge(0.8)) +
    scale_y_continuous("Mean hit rate of old / new responses (%)", expand = c(0,0), limit=c(0,100), breaks=seq(0,100,5)) + # define y axis properties
    #coord_cartesian(ylim=c(0,100)) +
    scale_fill_manual(values=c("#585858","#888888","#B8B8B8","#DCDCDC")) # color of bars
    
  # add significance to plot
    if (pval_accuracy_isgo < 0.001){
      go_nogo_accuracy_asteriks="***"  
    } else if (pval_accuracy_isgo < 0.01){
      go_nogo_accuracy_asteriks="**"
    } else if (pval_accuracy_isgo < 0.05){
      go_nogo_accuracy_asteriks="*"
    } else if (pval_accuracy_isgo < 0.07){
      go_nogo_accuracy_asteriks="+"
    } else {go_nogo_accuracy_asteriks=""}
    
  # show significance level on the plot, if difference between go and nogo is significant
    if (pval_accuracy_isgo<0.07) {
      Lines_hight=96
      i=1.5
      tmp_df=data.frame(x_val=c(i-1/2,i-1/2,i+1/2,i+1/2),y_val=c(Lines_hight-1,Lines_hight,Lines_hight,Lines_hight-1)) # define shape of an open rectangle above the bars
      plot_accuracy = plot_accuracy +
        geom_line(aes(x=c(tmp_df$x_val,tmp_df$x_val),y = c(tmp_df$y_val,tmp_df$y_val))) + # draw open rectangle
        geom_line(aes(x=tmp_df$x_val[c(2:3, 2:3, 2:3, 2:3)],y = tmp_df$y_val[c(2:3,2:3,2:3,2:3)])) + # draw open rectangle
        annotate("text", x = i, y = Lines_hight+1, label = (go_nogo_accuracy_asteriks) ,size=12) # differential effect significance asteriks
    }
  # show the plot
  print(plot_accuracy)
  
  # RT analysis
  print('- - - - - - - - -')
  print('RT old/new statistics - only correct responses!')
  print('- - - - - - - - -')
  # Subset relevant data
  recognition_probe_items_correct_isOld= subset(recognition_probe_items,recognition_probe_items$IsCorrectAnsOld==1)
  recognition_probe_items_correct_isGo= subset(recognition_probe_items,recognition_probe_items$IsCorrectAnsGo==1)
  
  # Compute means
  means_isold_RT=c(mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all go items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all nogo items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,PairType=="HH"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all HH items  
                         mean(with(data=subset(recognition_probe_items_correct_isOld,PairType=="HM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all HM items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,PairType=="LM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all LM items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,PairType=="LL"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all LL items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="HH"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # HH Go items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="HM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # HM Go items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="LM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # LM Go items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="LL"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # LL Go items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="HH"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # HH NoGo items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="HM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # HM NoGo items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="LM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # LM NoGo items
                         mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="LL"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # LL NoGo items
                         mean(with(data=recognition_probe_items_correct_isOld, tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T) # all items
  )
  means_isold_RT=round(means_isold_RT,3)
  gonogo_category=c('Go','NoGo','All','All','All','All','Go','Go','Go','Go','NoGo','NoGo','NoGo','NoGo','All')
  value_category=c('All','All','HH','HM','LM','LL','HH','HM','LM','LL','HH','HM','LM','LL','All')
  RT_all_means=data.frame(row.names = c('HH','HM','LM','LL','All'))
  RT_all_means$Go=means_isold_RT[c(7:10,1)]
  RT_all_means$NoGo=means_isold_RT[c(11:14,2)]
  RT_all_means$All=means_isold_RT[c(3:6,15)]
  
  # Compute standard deviation values
  sd_isold_RT=c(sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all go items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all nogo items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,PairType=="HH"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all HH items  
                   sd(with(data=subset(recognition_probe_items_correct_isOld,PairType=="HM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all HM items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,PairType=="LM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all LM items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,PairType=="LL"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # all LL items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="HH"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # HH Go items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="HM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # HM Go items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="LM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # LM Go items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="LL"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # LL Go items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="HH"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # HH NoGo items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="HM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # HM NoGo items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="LM"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # LM NoGo items
                   sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="LL"), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T), # LL NoGo items
                   sd(with(data=recognition_probe_items_correct_isOld, tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T) # all items
  )
  sd_isold_RT=round(sd_isold_RT,3)
  RT_all_sd=data.frame(row.names = c('HH','HM','LM','LL','All'))
  RT_all_sd=data.frame(row.names = c('HH','HM','LM','LL','All'))
  RT_all_sd$Go=sd_isold_RT[c(7:10,1)]
  RT_all_sd$NoGo=sd_isold_RT[c(11:14,2)]
  RT_all_sd$All=sd_isold_RT[c(3:6,15)]

  # Compute standard error values
  se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }
  se_isold_RT=c(se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all go items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all nogo items
                se(with(data=subset(recognition_probe_items_correct_isOld,PairType=="HH"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all HH items  
                se(with(data=subset(recognition_probe_items_correct_isOld,PairType=="HM"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all HM items
                se(with(data=subset(recognition_probe_items_correct_isOld,PairType=="LM"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all LM items
                se(with(data=subset(recognition_probe_items_correct_isOld,PairType=="LL"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all LL items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="HH"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # HH Go items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="HM"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # HM Go items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="LM"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # LM Go items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & PairType=="LL"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # LL Go items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="HH"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # HH NoGo items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="HM"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # HM NoGo items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="LM"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # LM NoGo items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & PairType=="LL"), tapply(RT_isOld, subjectID, mean, na.rm=T))), # LL NoGo items
                se(with(data=recognition_probe_items_correct_isOld, tapply(RT_isOld, subjectID, mean, na.rm=T))) # all items
  )
  se_isold_RT=round(se_isold_RT,3)
  RT_all_se=data.frame(row.names = c('HH','HM','LM','LL','All'))
  RT_all_se$Go=se_isold_RT[c(7:10,1)]
  RT_all_se$NoGo=se_isold_RT[c(11:14,2)]
  RT_all_se$All=se_isold_RT[c(3:6,15)]
  
  # combine descriptive results
  descriptive_results_RT=RT_all_means
  descriptive_results_RT$Go=paste(RT_all_means$Go," (", RT_all_sd$Go, ")",sep="")
  descriptive_results_RT$NoGo=paste(RT_all_means$NoGo," (", RT_all_sd$NoGo, ")",sep="")
  descriptive_results_RT$All=paste(RT_all_means$All," (", RT_all_sd$All, ")",sep="")
  
 # print descriptive results
  print("RT old/new - mean RT in seconds (standard deviation)")
  print("Only for items from the Go vs. NoGo probe comparisons")
  print("And only for correct responses")
  print(descriptive_results_RT)
  
  # linear regression analysis
  model_RT_by_item_type=lmer(RT_isOld ~ isGo. + (1|subjectID),data=recognition_probe_items_correct_isOld,na.action=na.omit)
  model_RT_by_item_type_and_value=lmer(RT_isOld ~ isGo. + PairType + (1|subjectID),data=recognition_probe_items_correct_isOld,na.action=na.omit) 
  model_RT_by_item_type_int_value=lmer(RT_isOld ~ isGo. * PairType + (1|subjectID),data=recognition_probe_items_correct_isOld,na.action=na.omit) 
  RT_main_effect_value=anova(model_RT_by_item_type, model_RT_by_item_type_and_value) # check the main effect of value category
  RT_interaction_effect=anova(model_RT_by_item_type_and_value, model_RT_by_item_type_int_value) # check effect of interaction
  summary_model_RT_by_item_type_and_value=summary(model_RT_by_item_type_and_value)
  pval_RT_isgo=summary_model_RT_by_item_type_and_value$coefficients[2,"Pr(>|t|)"]/2
  pval_RT_value_category=RT_main_effect_value$`Pr(>Chisq)`[2]/2
  pval_RT_interaction=RT_interaction_effect$`Pr(>Chisq)`[2]
  
  # print logistic regression analysis results
  print("results for RT")
  print("all Go vs. NoGo")
  print('one-sided p value:')
  print(pval_RT_isgo)
  print("value category (HH / HM / LM / LL)")
  print('one-sided p value:')
  print(pval_RT_value_category)
  print("Interaction go / nogo and value category:")
  print('two-sided p value:')
  print(pval_RT_interaction)

  # create dataframe for plot
  df_RT=data.frame(row.names = 1:8)
  df_RT$value_level=c("High","High","Medium-high","Medium-high","Medium-low","Medium-low", "Low","Low")
  df_RT$item_type=c("Go items","NoGo items","Go items","NoGo items","Go items","NoGo items","Go items","NoGo items")
  df_RT$means=c(RT_all_means$Go[1],RT_all_means$NoGo[1],RT_all_means$Go[2],RT_all_means$NoGo[2],RT_all_means$Go[3],RT_all_means$NoGo[3],RT_all_means$Go[4],RT_all_means$NoGo[4])
  df_RT$sd=c(RT_all_sd$Go[1],RT_all_sd$NoGo[1],RT_all_sd$Go[2],RT_all_sd$NoGo[2],RT_all_sd$Go[3],RT_all_sd$NoGo[3],RT_all_sd$Go[4],RT_all_sd$NoGo[4])
  df_RT$se=c(RT_all_se$Go[1],RT_all_se$NoGo[1],RT_all_se$Go[2],RT_all_se$NoGo[2],RT_all_se$Go[3],RT_all_se$NoGo[3],RT_all_se$Go[4],RT_all_se$NoGo[4])
  
  df_RT$value_level = factor(df_RT$value_level, levels = df_RT$value_level) # for the value levels to be shown on the graph in the correct order (and not sorted by ABC)
  
  # plot properties
  plot_RT= ggplot(data=df_RT, aes(x=item_type, y=means, fill=value_level)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8) +
    theme_bw() + # white background
    theme(legend.position="top",legend.title=element_blank()) + # position legend
    theme(axis.title.x=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank()) + # axis and background formating
    theme(aspect.ratio=1.2, text = element_text(size=24)) + # font size
    geom_errorbar(aes(ymin=means-se, ymax=means+se), width=1/6, position=position_dodge(0.8)) +
    scale_y_continuous("Mean RT of correct old / new responses (seconds)",limit=c(0,2), breaks=seq(0,2,0.02)) + # define y axis properties
    coord_cartesian(ylim=c(1.2,1.66)) +
    scale_fill_manual(values=c("#585858","#888888","#B8B8B8","#DCDCDC")) # color of bars
  
  # add significance to plot
  if (pval_RT_isgo < 0.001){
    go_nogo_RT_asteriks="***"  
  } else if (pval_RT_isgo < 0.01){
    go_nogo_RT_asteriks="**"
  } else if (pval_RT_isgo < 0.05){
    go_nogo_RT_asteriks="*"
  } else if (pval_RT_isgo < 0.07){
    go_nogo_RT_asteriks="+"
  } else {go_nogo_RT_asteriks=""}
  
  # show significance level on the plot, if difference between go and nogo is significant
  if (pval_RT_isgo<0.07) {
    Lines_hight=1.66
    i=1.5
    tmp_df=data.frame(x_val=c(i-1/2,i-1/2,i+1/2,i+1/2),y_val=c(Lines_hight-0.01,Lines_hight,Lines_hight,Lines_hight-0.01)) # define shape of an open rectangle above the bars
    plot_RT = plot_RT +
      geom_line(aes(x=c(tmp_df$x_val, tmp_df$x_val),y = c(tmp_df$y_val,tmp_df$y_val))) + # draw open rectangle
      geom_line(aes(x=tmp_df$x_val[c(2:3, 2:3, 2:3, 2:3)],y = tmp_df$y_val[c(2:3,2:3,2:3,2:3)])) + # draw open rectangle
      annotate("text", x = i, y = Lines_hight+0.01, label = (go_nogo_RT_asteriks) ,size=12) # differential effect significance asteriks
  }
  # plot the data
  print(plot_RT)
  
  # create a table for all the descriptive results in the current session - both accuracy and RT
  descriptive_results_All=data.frame(row.names=1:15)
  descriptive_results_All$item_type=c("Go","Go","Go","Go","Go","NoGo","NoGo","NoGo","NoGo","NoGo","All","All","All","All","All")
  descriptive_results_All$value_level=c("HH","HM","LM","LL","ALL")
  descriptive_results_All$accuracy=c(descriptive_results_accuracy$Go, descriptive_results_accuracy$NoGo, descriptive_results_accuracy$All)
  descriptive_results_All$RT=c(descriptive_results_RT$Go, descriptive_results_RT$NoGo, descriptive_results_RT$All)