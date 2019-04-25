# Required R package
library(lme4)
library(lmerTest)
library(ggplot2)
library(gridExtra)
library(reshape2)

# Clear workspace
rm(list=ls())

# choose sessions
sessionNum_probe = 2
sessionNum_recognition = 2

experiment_name="bmem_snacks"

# Define the local path where the data can be found
# you will need to set the correct path whre the RData file with the data is located
input_path="DEFINE PATH OF DATA HERE"

# load data
filename=paste(input_path,experiment_name,"_probe", sessionNum_probe, "_recognition", sessionNum_recognition, ".Rdata",sep="")
load(file=filename)

# organize data
probe_data$PairType3 = factor(probe_data$PairType2, levels = c("Low_Value","High_Value"))
probe_only_remembered_items$PairType3 = factor(probe_only_remembered_items$PairType2, levels = c("Low_Value","High_Value"))

# check that number of participants is ok
count_unique = function(x) {length(unique(x))}
num_participants_probe=count_unique(probe_data$subjectID)
print(paste("Number of participants in probe data: ", num_participants_probe, " Session: ", sessionNum_probe, SEP=""))
num_participants_recognition=count_unique(recognition_data$subjectID)
print(paste("Number of participants in recognition data: ", num_participants_recognition, " Session: ", sessionNum_recognition, SEP=""))

# analyze relationships between choices and memory, on a trial by trial basis

# probe analysis HV with regressors accuracy Go, 1-accuracy NoGo
HV_intercept_model=glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data,probe_data$PairType2=='High_Value' & !is.na(probe_data$Accuracy_old_goItem) & !is.na(probe_data$inv_Accuracy_old_nogoItem)),na.action=na.omit,family=binomial)
HV_accuracy_model=glmer(Outcome ~ 1 + Accuracy_old_goItem + inv_Accuracy_old_nogoItem + (1|subjectID),data=subset(probe_data,probe_data$PairType2=='High_Value' & !is.na(probe_data$Accuracy_old_goItem) & !is.na(probe_data$inv_Accuracy_old_nogoItem)),na.action=na.omit,family=binomial)
HV_accuracy=summary(HV_accuracy_model)
HV_accuracy_joint=anova(HV_intercept_model, HV_accuracy_model)

# probe analysis LV with regressors accuracy Go, 1-accuracy NoGo
LV_intercept_model=glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data,probe_data$PairType2=='Low_Value' & !is.na(probe_data$Accuracy_old_goItem) & !is.na(probe_data$inv_Accuracy_old_nogoItem)),na.action=na.omit,family=binomial)
LV_accuracy_model=glmer(Outcome ~ 1 + Accuracy_old_goItem + inv_Accuracy_old_nogoItem + (1|subjectID),data=subset(probe_data,probe_data$PairType2=='Low_Value' & !is.na(probe_data$Accuracy_old_goItem) & !is.na(probe_data$inv_Accuracy_old_nogoItem)),na.action=na.omit,family=binomial)
LV_accuracy=summary(LV_accuracy_model)
LV_accuracy_joint=anova(LV_intercept_model, LV_accuracy_model)

# interaction between memory regressors and the value level (high-value / low-value)
intercept_and_value_model=glmer(Outcome ~ 1 + Accuracy_old_goItem + inv_Accuracy_old_nogoItem + PairType3 + (1|subjectID),data=subset(probe_data,probe_data$PairType<=2),na.action=na.omit,family=binomial)
interaction_accuracy_model=glmer(Outcome ~ 1 + Accuracy_old_goItem + inv_Accuracy_old_nogoItem + PairType3 + Accuracy_old_goItem:PairType3 + inv_Accuracy_old_nogoItem:PairType3 + (1|subjectID),data=subset(probe_data,probe_data$PairType<=2),na.action=na.omit,family=binomial)
interaction_accuracy=summary(interaction_accuracy_model)
interaction_accuracy_joint=anova(intercept_and_value_model, interaction_accuracy_model)

# probe analysis HV with regressor RT difference NoGo-Go only for pairs where both items were remembered
HV_RT=summary(glmer(Outcome ~ 1 + RT_old_nogo_minus_go + (1|subjectID),data=subset(probe_only_remembered_items,(probe_only_remembered_items$PairType2=='High_Value')),na.action=na.omit,family=binomial)) 

# probe analysis LV with regressor RT difference NoGo-Go only for pairs where both items were remembered
LV_RT=summary(glmer(Outcome ~ 1 + RT_old_nogo_minus_go + (1|subjectID),data=subset(probe_only_remembered_items,(probe_only_remembered_items$PairType2=='Low_Value')),na.action=na.omit,family=binomial)) 

# probe analysis with regressors value (HV/LV), RT difference NoGo-Go and the interaction, only for pairs where both items were remembered
interaction_RT=summary(glmer(Outcome ~ 1 + PairType3 * RT_old_nogo_minus_go + (1|subjectID),data=subset(probe_only_remembered_items,(probe_only_remembered_items$PairType<=2)),na.action=na.omit,family=binomial)) 

# create a table with all results
df_results=data.frame(row.names=1:3)
#df_results$probe_session=sessionNum_probe
#df_results$recognition_session=sessionNum_recognition
df_results$value_level=c("HV", "LV", "Interaction with value level")
# accuracy Go
accuracy_Go_estimate=c(HV_accuracy$coefficients[2,1], LV_accuracy$coefficients[2,1],interaction_accuracy$coefficients[5,1])
accuracy_Go_odds_ratio=round(exp(accuracy_Go_estimate),3)
accuracy_Go_CI_min=round(exp(c(HV_accuracy$coefficients[2,1], LV_accuracy$coefficients[2,1],interaction_accuracy$coefficients[5,1])-(1.96*c(HV_accuracy$coefficients[2,2], LV_accuracy$coefficients[2,2],interaction_accuracy$coefficients[5,2]))),3)
accuracy_Go_CI_max=round(exp(c(HV_accuracy$coefficients[2,1], LV_accuracy$coefficients[2,1],interaction_accuracy$coefficients[5,1])+(1.96*c(HV_accuracy$coefficients[2,2], LV_accuracy$coefficients[2,2],interaction_accuracy$coefficients[5,2]))),3)
accuracy_Go_one_sided_p=round(c(HV_accuracy$coefficients[2,4], LV_accuracy$coefficients[2,4],interaction_accuracy$coefficients[5,4])/2,3)
df_results$accuracy_Go_P=accuracy_Go_one_sided_p
df_results$accuracy_Go_odds_ratio=accuracy_Go_odds_ratio
df_results$accuracy_Go_CI = paste(accuracy_Go_CI_min, '-',accuracy_Go_CI_max,sep="")
# inverse accuracy NoGo
inverse_accuracy_NoGo_estimate=c(HV_accuracy$coefficients[3,1], LV_accuracy$coefficients[3,1],interaction_accuracy$coefficients[6,1])
inverse_accuracy_NoGo_odds_ratio=round(exp(inverse_accuracy_NoGo_estimate),3)
inverse_accuracy_NoGo_CI_min=round(exp(c(HV_accuracy$coefficients[3,1], LV_accuracy$coefficients[3,1],interaction_accuracy$coefficients[6,1])-(1.96*c(HV_accuracy$coefficients[3,2], LV_accuracy$coefficients[3,2],interaction_accuracy$coefficients[6,2]))),3)
inverse_accuracy_NoGo_CI_max=round(exp(c(HV_accuracy$coefficients[3,1], LV_accuracy$coefficients[3,1],interaction_accuracy$coefficients[6,1])+(1.96*c(HV_accuracy$coefficients[3,2], LV_accuracy$coefficients[3,2],interaction_accuracy$coefficients[6,2]))),3)
inverse_accuracy_NoGo_one_sided_p=round(c(HV_accuracy$coefficients[3,4], LV_accuracy$coefficients[3,4],interaction_accuracy$coefficients[6,4])/2,3)
df_results$inverse_accuracy_NoGo_P=inverse_accuracy_NoGo_one_sided_p
df_results$inverse_accuracy_NoGo_odds_ratio=inverse_accuracy_NoGo_odds_ratio
df_results$inverse_accuracy_NoGo_CI = paste(inverse_accuracy_NoGo_CI_min, '-',inverse_accuracy_NoGo_CI_max,sep="")
# RT
RT_NoGo_minus_Go_estimate=c(HV_RT$coefficients[2,1], LV_RT$coefficients[2,1],interaction_RT$coefficients[4,1])
RT_NoGo_minus_Go_odds_ratio=round(exp(RT_NoGo_minus_Go_estimate),3)
RT_NoGo_minus_Go_CI_min=round(exp(c(HV_RT$coefficients[2,1], LV_RT$coefficients[2,1],interaction_RT$coefficients[4,1])-(1.96*c(HV_RT$coefficients[2,2], LV_RT$coefficients[2,2],interaction_RT$coefficients[4,2]))),3)
RT_NoGo_minus_Go_CI_max=round(exp(c(HV_RT$coefficients[2,1], LV_RT$coefficients[2,1],interaction_RT$coefficients[4,1])+(1.96*c(HV_RT$coefficients[2,2], LV_RT$coefficients[2,2],interaction_RT$coefficients[4,2]))),3)
RT_NoGo_minus_Go_one_sided_p=round(c(HV_RT$coefficients[2,4], LV_RT$coefficients[2,4],interaction_RT$coefficients[4,4])/2,3)
df_results$RT_NoGo_minus_Go_P=RT_NoGo_minus_Go_one_sided_p
df_results$RT_NoGo_minus_Go_odds_ratio=RT_NoGo_minus_Go_odds_ratio
df_results$RT_NoGo_minus_Go_CI = paste(RT_NoGo_minus_Go_CI_min, '-',RT_NoGo_minus_Go_CI_max,sep="")

print(df_results)

# visualization
# combinations between go and nogo accuracy - both together
# means
mean_props_combined=c(mean(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==probe_data$Accuracy_old_nogoItem & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HV both remembered or both forgotten
                      mean(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==1 & probe_data$Accuracy_old_nogoItem==0 & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HV Go remembered NoGo forgotten
                      mean(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==0 & probe_data$Accuracy_old_nogoItem==1 & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HV Go forgotten NoGo remembered
                      mean(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==probe_data$Accuracy_old_nogoItem & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # LV both remembered or both forgotten
                      mean(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==1 & probe_data$Accuracy_old_nogoItem==0 & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # LV Go remembered NoGo forgotten
                      mean(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==0 & probe_data$Accuracy_old_nogoItem==1 & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T)) # LV Go forgotten NoGo remembered

# SD values
SD_props_combined=c(sd(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==probe_data$Accuracy_old_nogoItem & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HV both remembered or both forgotten
                    sd(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==1 & probe_data$Accuracy_old_nogoItem==0 & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HV Go remembered NoGo forgotten
                    sd(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==0 & probe_data$Accuracy_old_nogoItem==1 & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HV Go forgotten NoGo remembered
                    sd(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==probe_data$Accuracy_old_nogoItem & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # LV both remembered or both forgotten
                    sd(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==1 & probe_data$Accuracy_old_nogoItem==0 & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # LV Go remembered NoGo forgotten
                    sd(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==0 & probe_data$Accuracy_old_nogoItem==1 & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T)) # LV Go forgotten NoGo remembered

# SE values
se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }
SE_props_combined=c(se(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==probe_data$Accuracy_old_nogoItem & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T))), # HV both remembered or both forgotten
                    se(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==1 & probe_data$Accuracy_old_nogoItem==0 & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T))), # HV Go remembered NoGo forgotten
                    se(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==0 & probe_data$Accuracy_old_nogoItem==1 & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T))), # HV Go forgotten NoGo remembered
                    se(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==probe_data$Accuracy_old_nogoItem & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T))), # LV both remembered or both forgotten
                    se(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==1 & probe_data$Accuracy_old_nogoItem==0 & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T))), # LV Go remembered NoGo forgotten
                    se(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==0 & probe_data$Accuracy_old_nogoItem==1 & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)))) # LV Go forgotten NoGo remembered

num_relevant_trials_combined = c(sum(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==probe_data$Accuracy_old_nogoItem & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T), # HV both remembered or both forgotten
                                 sum(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==1 & probe_data$Accuracy_old_nogoItem==0 & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T), # HV Go remembered NoGo forgotten
                                 sum(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==0 & probe_data$Accuracy_old_nogoItem==1 & probe_data$PairType==1 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T), # HV Go forgotten NoGo remembered
                                 sum(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==probe_data$Accuracy_old_nogoItem & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T), # LV both remembered or both forgotten
                                 sum(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==1 & probe_data$Accuracy_old_nogoItem==0 & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T), # LV Go remembered NoGo forgotten
                                 sum(with(data=subset(probe_data,probe_data$Accuracy_old_goItem==0 & probe_data$Accuracy_old_nogoItem==1 & probe_data$PairType==2 & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T)) # LV Go forgotten NoGo remembered

# put values in data frame for plot
df_accuracy_regressors_combined=data.frame(row.names = 1:6)
df_accuracy_regressors_combined$value_level[1:3] = "High value"
df_accuracy_regressors_combined$value_level[4:6] = "Low value"
df_accuracy_regressors_combined$accuracy_Go = df_accuracy_regressors_combined$accuracy_NoGo = factor(c("Both", "Remembered", "Forgotten"), levels = c("Remembered", "Forgotten", "Both"))
df_accuracy_regressors_combined$accuracy_NoGo = factor(c("Both", "Forgotten", "Remembered"), levels = c("Remembered", "Forgotten", "Both"))
df_accuracy_regressors_combined$means = mean_props_combined*100
df_accuracy_regressors_combined$SEs = SE_props_combined*100
df_accuracy_regressors_combined$SDs = SD_props_combined*100
df_accuracy_regressors_combined$category = factor(c("Both remembered / forgotten", "Go remembered, NoGo forgotten", "Go forgotten, NoGo remembered"), levels = c("Go remembered, NoGo forgotten", "Both remembered / forgotten", "Go forgotten, NoGo remembered"))
df_accuracy_regressors_combined$considered_trials = num_relevant_trials_combined

# plot accuracy Go & NoGo
plot_accuracy_regressors_combined= ggplot(data=df_accuracy_regressors_combined, aes(x=value_level, y=means, fill=category)) +
  geom_bar(stat="identity", position=position_dodge(0.7), width=0.7) +
  theme_bw() + # white background
  theme(axis.title.x=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank()) + # axis and background formating
  theme(aspect.ratio=1.2, text = element_text(size=24)) + # font size
  theme(legend.position="bottom", legend.direction = "vertical", legend.title=element_blank()) + # position legend
  geom_text(aes(label = paste(considered_trials), y=8), size=6, position = position_dodge(0.7), angle = 90) +
  geom_errorbar(aes(ymin=means-SEs, ymax=means+SEs), width=1/4, position=position_dodge(0.7)) +
  scale_y_continuous("Mean proportion of Go item choices", expand = c(0,0), limit=c(0,108), breaks=seq(0,100,5)) + # define y axis properties
  scale_fill_manual(values=c("Go remembered, NoGo forgotten"="#585858","Both remembered / forgotten"="#909090", "Go forgotten, NoGo remembered"="#DCDCDC")) # color of bars

# add significance asteriks
df_accuracy_significance = data.frame(row.names = 1:3)
df_accuracy_significance$pvals = c(HV_accuracy_joint$`Pr(>Chisq)`[2]/2, LV_accuracy_joint$`Pr(>Chisq)`[2]/2, interaction_accuracy_joint$`Pr(>Chisq)`[2]/2) # HV, LV, interaction
df_accuracy_significance$asteriks=""
df_accuracy_significance$asteriks[df_accuracy_significance$pvals<0.07]="+"
df_accuracy_significance$asteriks[df_accuracy_significance$pvals<0.05]="*"
df_accuracy_significance$asteriks[df_accuracy_significance$pvals<0.01]="**"
df_accuracy_significance$asteriks[df_accuracy_significance$pvals<0.001]="***"
Lines_hight=95

# add significance for HV and LV
for (i in 1:3){
  if (df_accuracy_significance$asteriks[i]!="") {
    if (i==3){ # interaction
      Lines_hight = Lines_hight + 7
      tmp_df=data.frame(x_val=c(1,1,2,2),y_val=c(Lines_hight,Lines_hight+1,Lines_hight+1,Lines_hight)) # define shape of an open rectangle above the bar
      tmp_df$category=factor(NA, levels = c("Go remembered, NoGo forgotten", "Both remembered / forgotten", "Go forgotten, NoGo remembered"))
      plot_accuracy_regressors_combined = plot_accuracy_regressors_combined +
        geom_line(data = tmp_df, aes(x=x_val,y = y_val)) + # draw open rectangle
        annotate("text", x = 1.5, y = Lines_hight+1, label = (df_accuracy_significance$asteriks[i]),size=14) # differential effect significance asteriks
    } else{
      tmp_df=data.frame(x_val=c(i-0.25,i-0.25,i+0.25,i+0.25),y_val=c(Lines_hight,Lines_hight+1,Lines_hight+1,Lines_hight)) # define shape of an open rectangle above the bar
      tmp_df$category=factor(NA,levels = c("Go remembered, NoGo forgotten", "Both remembered / forgotten", "Go forgotten, NoGo remembered"))
      plot_accuracy_regressors_combined = plot_accuracy_regressors_combined +
        geom_line(data = tmp_df, aes(x=x_val,y = y_val)) + # draw open rectangle
        annotate("text", x = i, y = Lines_hight+1, label = (df_accuracy_significance$asteriks[i]),size=14) # differential effect significance asteriks
    }
  }
}

# show the plot
print(plot_accuracy_regressors_combined)

# plot RT
RT_data_agg=as.data.frame(aggregate(probe_only_remembered_items,by=list(probe_only_remembered_items$subjectID_num, probe_only_remembered_items$PairType, probe_only_remembered_items$Outcome), mean, na.rm=TRUE))
Data_by_sub=melt(RT_data_agg,id = c("subjectID_num", "PairType", "Outcome") ,"RT_old_nogo_minus_go")
RT_means_for_plot=with(data=Data_by_sub, tapply(value, list(PairType, Outcome), mean))
RT_SEs_for_plot=with(data=Data_by_sub, tapply(value, list(PairType, Outcome), se))

df_RT_for_plot=data.frame(row.names = 1:4)
df_RT_for_plot$value_level[1:2]="High-value"
df_RT_for_plot$value_level[3:4]="Low-value"
df_RT_for_plot$item_type=c("Go chosen", "NoGo chosen", "Go chosen", "NoGo chosen")
df_RT_for_plot$means = c(RT_means_for_plot[1,2], RT_means_for_plot[1,1], RT_means_for_plot[2,2], RT_means_for_plot[2,1])
df_RT_for_plot$SEs = c(RT_SEs_for_plot[1,2], RT_SEs_for_plot[1,1], RT_SEs_for_plot[2,2], RT_SEs_for_plot[2,1])

# plot RT-choices
plot_RT= ggplot(data=df_RT_for_plot, aes(x=value_level, y=means, fill=item_type)) +
  geom_bar(stat="identity", position=position_dodge(0.7), width=0.7) +
  theme_bw() + # white background
  theme(legend.position="bottom",legend.title=element_blank()) + # position legend
  theme(axis.title.x=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank()) + # axis and background formating
  theme(aspect.ratio=1.2, text = element_text(size=24)) + # font size
  geom_errorbar(aes(ymin=means-SEs, ymax=means+SEs), width=1/4, position=position_dodge(0.7)) +
  scale_y_continuous(expression(paste(Delta, "RT (NoGo recognition RT - Go recognition RT)")), expand = c(0,0), limit=c(-0.24,0.4), breaks=seq(-0.24,0.4,0.04)) + # define y axis properties
  scale_fill_manual(values=c("#585858","#DCDCDC")) + # color of bars
  geom_abline(intercept = 0,slope=0,linetype =2, size = 1,aes()) # add line for y=0

# add significance asteriks
df_RT_significance = data.frame(row.names = 1:3)
df_RT_significance$pvals = c(HV_RT$coefficients[2, "Pr(>|z|)"]/2, LV_RT$coefficients[2, "Pr(>|z|)"]/2, interaction_RT$coefficients[4,"Pr(>|z|)"]/2) # HV, LV, interaction
df_RT_significance$estimates = c(HV_RT$coefficients[2, "Estimate"], LV_RT$coefficients[2, "Estimate"], interaction_RT$coefficients[4,"Estimate"]) # HV, LV, interaction)
df_RT_significance$asteriks=""
df_RT_significance$asteriks[df_RT_significance$pvals<0.07 & df_RT_significance$estimates>0]="+"
df_RT_significance$asteriks[df_RT_significance$pvals<0.05 & df_RT_significance$estimates>0]="*"
df_RT_significance$asteriks[df_RT_significance$pvals<0.01 & df_RT_significance$estimates>0]="**"
df_RT_significance$asteriks[df_RT_significance$pvals<0.001 & df_RT_significance$estimates>0]="***"
Lines_hight=0.33

# add significance for HV and LV
for (i in 1:3){
  if (df_RT_significance$asteriks[i]!="") {
    if (i==3){ # interaction
      Lines_hight = Lines_hight + 0.025
      tmp_df=data.frame(x_val=c(i/2-1/2,i/2-1/2,i/2+1/2,i/2+1/2),y_val=c(Lines_hight,Lines_hight+0.01,Lines_hight+0.01,Lines_hight)) # define shape of an open rectangle above the bar
      tmp_df$item_type=NA
      plot_RT = plot_RT +
        geom_line(data = tmp_df, aes(x=x_val,y = y_val)) + # draw open rectangle
        annotate("text", x = i/2, y = Lines_hight+0.01, label = (df_RT_significance$asteriks[i]),size=14) # differential effect significance asteriks
    } else{
      tmp_df=data.frame(x_val=c(i-0.7/4,i-0.7/4,i+0.7/4,i+0.7/4),y_val=c(Lines_hight,Lines_hight+0.01,Lines_hight+0.01,Lines_hight)) # define shape of an open rectangle above the bar
      tmp_df$item_type=NA
      plot_RT = plot_RT +
        geom_line(data = tmp_df, aes(x=x_val,y = y_val)) + # draw open rectangle
        annotate("text", x = i, y = Lines_hight+0.01, label = (df_RT_significance$asteriks[i]),size=14) # differential effect significance asteriks
      
    }
  }
}

# show the plot
print(plot_RT)