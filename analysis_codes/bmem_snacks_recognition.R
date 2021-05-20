# Required R package
library(lme4)
library(lmerTest)
library(ggplot2)

# Clear workspace
rm(list=ls())

# Define the local path where the data can be found
# you will need to set the correct path where the RData file with the data is located
input_path="DEFINE PATH OF DATA HERE"

experiment_name = "bmem_snacks"

# define session (options: 2 or 3 or 2:3)
sessions = 3

# load data
input_filename=paste(input_path, experiment_name, "_recognition_all_sessions.Rdata",sep="")  
load(file=input_filename)

# run over sessions
for (session_num in sessions){
  if (session_num==2) {
    session_str = "Session 2"
  } else if (session_num==3) {
    session_str = "Follow-up"
  }
  # subset data for current session
  recognition_data = subset(recognition_data_all_sessions, recognition_data_all_sessions$session == session_str)
  recognition_data$subjectID = factor(recognition_data$subjectID)
  
  # make sure recognition data was read correctly
  length_recognition_across_participants=tapply(recognition_data$order,recognition_data$subjectID,length)
  unique_length_recognition=unique(length_recognition_across_participants)
  if (length(unique_length_recognition) > 1){
    print("Not all participants have the same number of trials!")
    print(length_recognition_across_participants)
  }
  
  
  # subset only probe items and old items
  recognition_probe_items = subset(recognition_data,recognition_data$isProbeItem)
  recognition_old_items = subset(recognition_data,recognition_data$isOld.==1)
  
  # check data is ok- all participants have the correct number of trials with probe items
  correct_num_probe_items=24
  num_probe_items_across_participants=tapply(recognition_probe_items$isGo.,recognition_probe_items$subjectID,length)
  unique_num_probe_items=unique(num_probe_items_across_participants)
  if (length(unique_num_probe_items) > 1 || unique_num_probe_items!=correct_num_probe_items){
    print(paste("Not all participants have the correct number of probe items!\nCorrect number is", correct_num_probe_items))
  }
  
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
                          mean(with(data=subset(recognition_probe_items,IsHighValue==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all HV items
                          mean(with(data=subset(recognition_probe_items,IsHighValue==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all LV items
                          mean(with(data=subset(recognition_probe_items,isGo.==1 & IsHighValue==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HV Go items
                          mean(with(data=subset(recognition_probe_items,isGo.==1 & IsHighValue==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LV Go items
                          mean(with(data=subset(recognition_probe_items,isGo.==0 & IsHighValue==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HV NoGo items
                          mean(with(data=subset(recognition_probe_items,isGo.==0 & IsHighValue==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LV NoGo items
                          mean(with(data=recognition_probe_items, tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))) # all items
  )
  means_isold_accuracy=round(means_isold_accuracy*100,2)
  gonogo_category=c('Go','NoGo','All','All','Go','Go','NoGo','NoGo')
  value_category=c('All','All','High-value','Low-value','High-value','Low-value','High-value','Low-value')
  accuracy_all_means=data.frame(row.names = c('High-value','Low-value','All'))
  accuracy_all_means$Go=means_isold_accuracy[c(5,6,1)]
  accuracy_all_means$NoGo=means_isold_accuracy[c(7,8,2)]
  accuracy_all_means$All=means_isold_accuracy[c(3:4,9)]
  
  # compute SD values
  sd_isold_accuracy=c(sd(with(data=subset(recognition_probe_items,isGo.==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all go items
                      sd(with(data=subset(recognition_probe_items,isGo.==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all nogo items
                      sd(with(data=subset(recognition_probe_items,IsHighValue==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all HV items
                      sd(with(data=subset(recognition_probe_items,IsHighValue==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all LV items
                      sd(with(data=subset(recognition_probe_items,isGo.==1 & IsHighValue==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HV Go items
                      sd(with(data=subset(recognition_probe_items,isGo.==1 & IsHighValue==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LV Go items
                      sd(with(data=subset(recognition_probe_items,isGo.==0 & IsHighValue==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HV NoGo items
                      sd(with(data=subset(recognition_probe_items,isGo.==0 & IsHighValue==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LV NoGo items
                      sd(with(data=recognition_probe_items, tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))) # all items
  )
  sd_isold_accuracy=round(sd_isold_accuracy*100,2) # from proportion to percents
  accuracy_all_sd=data.frame(row.names = c('High-value','Low-value','All'))
  accuracy_all_sd$Go=sd_isold_accuracy[c(5,6,1)]
  accuracy_all_sd$NoGo=sd_isold_accuracy[c(7,8,2)]
  accuracy_all_sd$All=sd_isold_accuracy[c(3:4,9)]
  
  # Compute SE values
  se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }
  se_isold_accuracy=c(se(with(data=subset(recognition_probe_items,isGo.==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all go items
                      se(with(data=subset(recognition_probe_items,isGo.==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all nogo items
                      se(with(data=subset(recognition_probe_items,IsHighValue==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all HV items
                      se(with(data=subset(recognition_probe_items,IsHighValue==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # all LV items
                      se(with(data=subset(recognition_probe_items,isGo.==1 & IsHighValue==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HV Go items
                      se(with(data=subset(recognition_probe_items,isGo.==1 & IsHighValue==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LV Go items
                      se(with(data=subset(recognition_probe_items,isGo.==0 & IsHighValue==1), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # HV NoGo items
                      se(with(data=subset(recognition_probe_items,isGo.==0 & IsHighValue==0), tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))), # LV NoGo items
                      se(with(data=recognition_probe_items, tapply(IsCorrectAnsOld, subjectID, mean, na.rm=T))) # all items
  )
  se_isold_accuracy=round(se_isold_accuracy*100,2) # from proportion to percents
  accuracy_all_se=data.frame(row.names = c('High-value','Low-value','All'))
  accuracy_all_se$Go=se_isold_accuracy[c(5,6,1)]
  accuracy_all_se$NoGo=se_isold_accuracy[c(7,8,2)]
  accuracy_all_se$All=se_isold_accuracy[c(3:4,9)]
  
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
  recognition_probe_items$go.ind = 1*recognition_probe_items$isGo.
  recognition_probe_items$high.ind = 1*recognition_probe_items$IsHighValue
  
  accuracy_isold_results_all_items_with_interaction=summary(glmer(IsCorrectAnsOld ~ 1 + isGo.*IsHighValue + (1+high.ind|subjectID),data=recognition_probe_items,na.action=na.omit,family=binomial)) 
  if (session_num==2) {
    accuracy_isold_results_all_items=summary(glmer(IsCorrectAnsOld ~ 1 + isGo.+IsHighValue + (1+high.ind|subjectID),data=recognition_probe_items,na.action=na.omit,family=binomial)) 
    accuracy_isold_results_HV_items=summary(glmer(IsCorrectAnsOld ~ 1 + isGo. + (1|subjectID),data=subset(recognition_probe_items,(recognition_probe_items$IsHighValue)),na.action=na.omit,family=binomial))
    accuracy_isold_results_LV_items=summary(glmer(IsCorrectAnsOld ~ 1 + isGo. + (1|subjectID),data=subset(recognition_probe_items,(!recognition_probe_items$IsHighValue)),na.action=na.omit,family=binomial))
  }
  if (session_num ==3) {
    accuracy_isold_results_all_items=summary(glmer(IsCorrectAnsOld ~ 1 + isGo.+IsHighValue + (1+high.ind+go.ind|subjectID),data=recognition_probe_items,na.action=na.omit,family=binomial)) 
    accuracy_isold_results_HV_items=summary(glmer(IsCorrectAnsOld ~ 1 + isGo. + (1+go.ind||subjectID),data=subset(recognition_probe_items,(recognition_probe_items$IsHighValue)),na.action=na.omit,family=binomial))
    accuracy_isold_results_LV_items=summary(glmer(IsCorrectAnsOld ~ 1 + isGo. + (1+go.ind||subjectID),data=subset(recognition_probe_items,(!recognition_probe_items$IsHighValue)),na.action=na.omit,family=binomial))
  }
  
  # organize statistics in table
  results_isold_accuracy=rbind(accuracy_isold_results_all_items$coefficients[2,],accuracy_isold_results_all_items$coefficients[3,],accuracy_isold_results_all_items_with_interaction$coefficients[4,],accuracy_isold_results_HV_items$coefficients[2,],accuracy_isold_results_LV_items$coefficients[2,])
  results_isold_accuracy=as.data.frame(results_isold_accuracy)
  colnames(results_isold_accuracy)[4] = "two-sided p"
  results_isold_accuracy$odds_ratio=round(exp(results_isold_accuracy$Estimate),3)
  CI_min=exp(results_isold_accuracy$Estimate-1.96*results_isold_accuracy$`Std. Error`)
  CI_max=exp(results_isold_accuracy$Estimate+1.96*results_isold_accuracy$`Std. Error`)
  results_isold_accuracy$CI=paste(round(CI_min,3)," - ", round(CI_max,3), sep="")
  results_isold_accuracy$category=c('all Go vs. NoGo', 'all HV vs. LV', 'all interaction Go/NoGo and HV/LV', 'HV Go vs. NoGo', 'LV Go vs. NoGo')
  results_isold_accuracy$one_sided_p=round(results_isold_accuracy$`two-sided p`/2,3)
  
  print(results_isold_accuracy)
  
  # create dataframe for plot
  df_accuracy=data.frame(row.names = 1:4)
  df_accuracy$value_level=c("High-value","High-value", "Low-value","Low-value")
  df_accuracy$item_type=c("Go items","NoGo items","Go items","NoGo items")
  df_accuracy$means=c(accuracy_all_means$Go[1],accuracy_all_means$NoGo[1],accuracy_all_means$Go[2],accuracy_all_means$NoGo[2])
  df_accuracy$sd=c(accuracy_all_sd$Go[1],accuracy_all_sd$NoGo[1],accuracy_all_sd$Go[2],accuracy_all_sd$NoGo[2])
  df_accuracy$se=c(accuracy_all_se$Go[1],accuracy_all_se$NoGo[1],accuracy_all_se$Go[2],accuracy_all_se$NoGo[2])
  
  # plot the results
  plot_accuracy= ggplot(data=df_accuracy, aes(x=item_type, y=means, fill=value_level)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8) +
    theme_bw() + # white background
    #theme(legend.position="top",legend.title=element_blank()) + # position legend
    theme(legend.position="none",legend.title=element_blank()) + # position legend
    theme(axis.title.x=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank()) + # axis and background formating
    theme(aspect.ratio=1.2, text = element_text(size=24)) + # font size
    geom_errorbar(aes(ymin=means-se, ymax=means+se), width=1/6, position=position_dodge(0.8)) +
    scale_y_continuous("Mean hit rate of old/new responses", expand = c(0,0), limit=c(0,100), breaks=seq(0,100,5)) + # define y axis properties
    #coord_cartesian(ylim=c(0,100)) +
    scale_fill_manual(values=c("#585858","#D8D8D8")) # color of bars
  
  # add significance to plot
    go_nogo_accuracy_one_sided_p=results_isold_accuracy$`two-sided p`[1]/2
    if (go_nogo_accuracy_one_sided_p < 0.001){
      go_nogo_accuracy_asteriks="***"  
    } else if (go_nogo_accuracy_one_sided_p < 0.01){
      go_nogo_accuracy_asteriks="**"
    } else if (go_nogo_accuracy_one_sided_p < 0.05){
      go_nogo_accuracy_asteriks="*"
    } else if (go_nogo_accuracy_one_sided_p < 0.07){
      go_nogo_accuracy_asteriks="+"
    } else {go_nogo_accuracy_asteriks=""}
    
  # show significance level on the plot, if difference between go and nogo is significant
    if (go_nogo_accuracy_one_sided_p<0.07) {
      Lines_hight=96
      i=1.5
      tmp_df=data.frame(x_val=c(i-1/2,i-1/2,i+1/2,i+1/2),y_val=c(Lines_hight-1,Lines_hight,Lines_hight,Lines_hight-1)) # define shape of an open rectangle above the bars
      plot_accuracy = plot_accuracy +
        geom_line(aes(x=tmp_df$x_val,y = tmp_df$y_val)) + # draw open rectangle
        geom_line(aes(x=tmp_df$x_val[c(2:3, 2:3)],y = tmp_df$y_val[c(2:3,2:3)])) + # draw open rectangle
        annotate("text", x = i, y = Lines_hight+1, label = (go_nogo_accuracy_asteriks) ,size=12) # differential effect significance asteriks
    }
  # show the plot
  print(plot_accuracy)
  
  # Perform t-test- because this is the test that was pre-registered (for bmem_snacks2)
  accuracy_by_sub=tapply(recognition_probe_items$IsCorrectAnsOld,list(isGo=recognition_probe_items$isGo.,subject=recognition_probe_items$subjectID),mean,na.rm=T)
  ttest_results_accuracy=t.test(accuracy_by_sub[2,],accuracy_by_sub[1,],paired=TRUE,alternative="greater")

  
  # RT analysis
  print('- - - - - - - - -')
  print('RT old/new statistics - only correct responses!')
  print('- - - - - - - - -')
  # Subset relevant data
  recognition_probe_items_correct_isOld= subset(recognition_probe_items,recognition_probe_items$IsCorrectAnsOld==1)
  recognition_probe_items_correct_isGo= subset(recognition_probe_items,recognition_probe_items$IsCorrectAnsGo==1)
  # Compute means
  means_isold_RT=       c(mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all go items
                          mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all nogo items
                          mean(with(data=subset(recognition_probe_items_correct_isOld,IsHighValue==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all HV items
                          mean(with(data=subset(recognition_probe_items_correct_isOld,IsHighValue==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all LV items
                          mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & IsHighValue==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # HV Go items
                          mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & IsHighValue==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # LV Go items
                          mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & IsHighValue==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # HV NoGo items
                          mean(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & IsHighValue==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # LV NoGo items
                          mean(with(data=recognition_probe_items_correct_isOld, tapply(RT_isOld, subjectID, mean, na.rm=T))) # all items
  )
  means_isold_RT=round(means_isold_RT,3)
  RT_all_means=data.frame(row.names = c('High-value','Low-value','All'))
  RT_all_means$Go=means_isold_RT[c(5,6,1)]
  RT_all_means$NoGo=means_isold_RT[c(7,8,2)]
  RT_all_means$All=means_isold_RT[c(3:4,9)]
  
  # Compute standard deviation values
  sd_isold_RT=c(sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all go items
                sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all nogo items
                sd(with(data=subset(recognition_probe_items_correct_isOld,IsHighValue==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all HV items
                sd(with(data=subset(recognition_probe_items_correct_isOld,IsHighValue==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all LV items
                sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & IsHighValue==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # HV Go items
                sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & IsHighValue==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # LV Go items
                sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & IsHighValue==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # HV NoGo items
                sd(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & IsHighValue==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # LV NoGo items
                sd(with(data=recognition_probe_items_correct_isOld, tapply(RT_isOld, subjectID, mean, na.rm=T))) # all items
  )
  sd_isold_RT=round(sd_isold_RT,2) # from proportion to percents
  RT_all_sd=data.frame(row.names = c('High-value','Low-value','All'))
  RT_all_sd$Go=sd_isold_RT[c(5,6,1)]
  RT_all_sd$NoGo=sd_isold_RT[c(7,8,2)]
  RT_all_sd$All=sd_isold_RT[c(3:4,9)]
  
  # Compute SE values
  se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }
  se_isold_RT=c(se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all go items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all nogo items
                se(with(data=subset(recognition_probe_items_correct_isOld,IsHighValue==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all HV items
                se(with(data=subset(recognition_probe_items_correct_isOld,IsHighValue==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # all LV items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & IsHighValue==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # HV Go items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==1 & IsHighValue==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # LV Go items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & IsHighValue==1), tapply(RT_isOld, subjectID, mean, na.rm=T))), # HV NoGo items
                se(with(data=subset(recognition_probe_items_correct_isOld,isGo.==0 & IsHighValue==0), tapply(RT_isOld, subjectID, mean, na.rm=T))), # LV NoGo items
                se(with(data=recognition_probe_items_correct_isOld, tapply(RT_isOld, subjectID, mean, na.rm=T))) # all items
  )
  se_isold_RT=round(se_isold_RT,2) # from proportion to percents
  RT_all_se=data.frame(row.names = c('High-value','Low-value','All'))
  RT_all_se$Go=se_isold_RT[c(5,6,1)]
  RT_all_se$NoGo=se_isold_RT[c(7,8,2)]
  RT_all_se$All=se_isold_RT[c(3:4,9)]
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
  
  # linear mixed model analysis
  RT_isold_results_all_items_with_interaction=summary(lmer(RT_isOld ~ isGo.*IsHighValue + (1+go.ind+high.ind||subjectID),data=recognition_probe_items_correct_isOld,na.action=na.omit))
  RT_isold_results_HV_items=summary(lmer(RT_isOld ~ isGo.+ (1+go.ind|subjectID),data=subset(recognition_probe_items_correct_isOld,(recognition_probe_items_correct_isOld$IsHighValue)),na.action=na.omit)) 
  if (session_num ==2) {
    RT_isold_results_all_items=summary(lmer(RT_isOld ~ isGo.+IsHighValue + (1+go.ind+high.ind|subjectID),data=recognition_probe_items_correct_isOld,na.action=na.omit))
    RT_isold_results_LV_items=summary(lmer(RT_isOld ~ isGo.+ (1|subjectID),data=subset(recognition_probe_items_correct_isOld,(!recognition_probe_items_correct_isOld$IsHighValue)),na.action=na.omit)) 
  }
  if (session_num ==3) {
    RT_isold_results_all_items=summary(lmer(RT_isOld ~ isGo.+IsHighValue + (1+go.ind+high.ind||subjectID),data=recognition_probe_items_correct_isOld,na.action=na.omit))
    RT_isold_results_LV_items=summary(lmer(RT_isOld ~ isGo.+ (1+go.ind|subjectID),data=subset(recognition_probe_items_correct_isOld,(!recognition_probe_items_correct_isOld$IsHighValue)),na.action=na.omit)) 
  }  
  # organize statistics in table
  results_isold_RT=rbind(RT_isold_results_all_items$coefficients[2,],RT_isold_results_all_items$coefficients[3,],RT_isold_results_all_items_with_interaction$coefficients[4,],RT_isold_results_HV_items$coefficients[2,],RT_isold_results_LV_items$coefficients[2,])
  results_isold_RT=as.data.frame(results_isold_RT)
  colnames(results_isold_RT)[5] = "two-sided p"
  results_isold_RT$category=c('all Go vs. NoGo', 'all HV vs. LV', 'all interaction Go/NoGo and HV/LV', 'HV Go vs. NoGo', 'LV Go vs. NoGo')
  results_isold_RT$one_sided_p=round(results_isold_RT$`two-sided p`/2,3)
  
  print(results_isold_RT)
  
  # create dataframe for plot
  df_RT=data.frame(row.names = 1:4)
  df_RT$value_level=c("High-value","High-value", "Low-value","Low-value")
  df_RT$item_type=c("Go items","NoGo items","Go items","NoGo items")
  df_RT$means=c(RT_all_means$Go[1],RT_all_means$NoGo[1],RT_all_means$Go[2],RT_all_means$NoGo[2])
  df_RT$sd=c(RT_all_sd$Go[1],RT_all_sd$NoGo[1],RT_all_sd$Go[2],RT_all_sd$NoGo[2])
  df_RT$se=c(RT_all_se$Go[1],RT_all_se$NoGo[1],RT_all_se$Go[2],RT_all_se$NoGo[2])
  # plot properties
  plot_RT= ggplot(data=df_RT, aes(x=item_type, y=means, fill=value_level)) +
    geom_bar(stat="identity", position=position_dodge(), width=0.8) +
    theme_bw() + # white background
    theme(legend.position="top",legend.title=element_blank()) + # position legend
    theme(axis.title.x=element_blank(),legend.position = "none",axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank()) + # axis and background formating
    theme(aspect.ratio=1.2, text = element_text(size=24)) + # font size
    geom_errorbar(aes(ymin=means-se, ymax=means+se), width=1/6, position=position_dodge(0.8)) +
    scale_y_continuous("Mean RT of correct old / new responses",limit=c(0,2), breaks=seq(0,2,0.02)) + # define y axis properties
    coord_cartesian(ylim=c(1.2,1.66)) +
    scale_fill_manual(values=c("#585858","#D8D8D8")) # color of bars
  
  # add significance to plot
  go_nogo_RT_one_sided_p=results_isold_RT$`two-sided p`[1]/2
  if (go_nogo_RT_one_sided_p < 0.001){
    go_nogo_RT_asteriks="***"  
  } else if (go_nogo_RT_one_sided_p < 0.01){
    go_nogo_RT_asteriks="**"
  } else if (go_nogo_RT_one_sided_p < 0.05){
    go_nogo_RT_asteriks="*"
  } else if (go_nogo_RT_one_sided_p < 0.07){
    go_nogo_RT_asteriks="+"
  } else {go_nogo_RT_asteriks=""}
  
  # show significance level on the plot, if difference between go and nogo is significant
  if (go_nogo_RT_one_sided_p<0.07) {
    Lines_hight=1.66
    i=1.5
    tmp_df=data.frame(x_val=c(i-1/2,i-1/2,i+1/2,i+1/2),y_val=c(Lines_hight-0.01,Lines_hight,Lines_hight,Lines_hight-0.01)) # define shape of an open rectangle above the bars
    plot_RT = plot_RT +
      geom_line(aes(x=tmp_df$x_val,y = tmp_df$y_val)) + # draw open rectangle
      geom_line(aes(x=tmp_df$x_val[c(2:3, 2:3)],y = tmp_df$y_val[c(2:3,2:3)])) + # draw open rectangle
      annotate("text", x = i, y = Lines_hight+0.01, label = (go_nogo_RT_asteriks) ,size=12) # differential effect significance asteriks
  }
  # plot the data
  print(plot_RT)
  
  # perform t-test because this is the test that we pre-registered (for bmem_snacks2)
  RT_by_sub=tapply(recognition_probe_items_correct_isOld$RT_isOld,list(isGo=recognition_probe_items_correct_isOld$isGo.,subject=recognition_probe_items_correct_isOld$subjectID),mean,na.rm=T)
  ttest_RT=t.test(RT_by_sub[2,],RT_by_sub[1,],paired=TRUE,alternative="less") 
  
  # create a table for all the descriptive results in the current session - both accuracy and RT
  descriptive_results_All=data.frame(row.names=1:9)
  descriptive_results_All$item_type=c("Go","Go","Go","NoGo","NoGo","NoGo","All","All","All")
  descriptive_results_All$value_level=c("HV","LV","ALL")
  descriptive_results_All$accuracy=c(descriptive_results_accuracy$Go, descriptive_results_accuracy$NoGo, descriptive_results_accuracy$All)
  descriptive_results_All$RT=c(descriptive_results_RT$Go, descriptive_results_RT$NoGo, descriptive_results_RT$All)
  
} # end of sessions loop
