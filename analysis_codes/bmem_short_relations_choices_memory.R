# Required R package
library(lme4)
library(lmerTest)
library(ggplot2)
library(gridExtra)
library(reshape2)

# Clear workspace
rm(list=ls())

experiment_name="bmem_short"

# Define the local path where the data can be found
# you will need to set the correct path whre the RData file with the data is located
input_path="DEFINE PATH OF DATA HERE"

# load data
filename=paste(input_path,experiment_name,"_probe_recognition.Rdata",sep="")
load(file=filename)

# organize data
probe_data$PairType2 = factor(probe_data$PairType2, levels = c("Low_Value","High_Value"))
probe_only_remembered_items$PairType2 = factor(probe_only_remembered_items$PairType2, levels = c("Low_Value","High_Value"))
probe_data$PairType3 = factor(probe_data$PairType3, levels = c("HH","HM","LM","LL","Sanity1","Sanity2","Sanity3"))
probe_only_remembered_items$PairType3 = factor(probe_only_remembered_items$PairType3, levels = c("HH","HM","LM","LL","Sanity1","Sanity2","Sanity3"))

probe_data$isHH = probe_data$PairType3=="HH"
probe_only_remembered_items$isHH = probe_only_remembered_items$PairType3=="HH"

# check that number of participants is ok
count_unique = function(x) {length(unique(x))}
num_participants_probe=count_unique(probe_data$subjectID)
print(paste("Number of participants in probe data: ", num_participants_probe,sep = ""))
num_participants_recognition=count_unique(recognition_data$subjectID)
print(paste("Number of participants in recognition data: ", num_participants_recognition,sep = ""))

# analyze relations between choices and memory, on a trial by trial basis

# Accuracy
probe_data$accuracy_category = NA
probe_data$accuracy_category[(probe_data$Accuracy_old_goItem & probe_data$Accuracy_old_nogoItem) | (!probe_data$Accuracy_old_goItem & !probe_data$Accuracy_old_nogoItem)] = "Both remembered / forgotten"
probe_data$accuracy_category[probe_data$Accuracy_old_goItem & !probe_data$Accuracy_old_nogoItem] = "Go remembered, NoGo forgotten"
probe_data$accuracy_category[!probe_data$Accuracy_old_goItem & probe_data$Accuracy_old_nogoItem] = "Go forgotten, NoGo remembered"

# probe analysis HH with accuracy category
HH_accuracy_model=glmer(Outcome ~ 1 + accuracy_category + (1 + accuracy_category||subjectID),data=subset(probe_data,(probe_data$PairType3=='HH' & !is.na(probe_data$accuracy_category))),na.action=na.omit,family=binomial)
HH_accuracy=summary(HH_accuracy_model) 

# probe analysis HM with accuracy category
HM_accuracy_model=glmer(Outcome ~ 1 + accuracy_category + (1 + accuracy_category|subjectID),data=subset(probe_data,(probe_data$PairType3=='HM' & !is.na(probe_data$accuracy_category))),na.action=na.omit,family=binomial)
HM_accuracy=summary(HM_accuracy_model) 

# probe analysis LM with accuracy category
LM_accuracy_model=glmer(Outcome ~ 1 + accuracy_category + (1 + accuracy_category||subjectID),data=subset(probe_data,(probe_data$PairType3=='LM' & !is.na(probe_data$accuracy_category))),na.action=na.omit,family=binomial)
LM_accuracy=summary(LM_accuracy_model) 

# probe analysis LL with accuracy category
LL_accuracy_model=glmer(Outcome ~ 1 + accuracy_category + (1|subjectID),data=subset(probe_data,(probe_data$PairType3=='LL' & !is.na(probe_data$accuracy_category))),na.action=na.omit,family=binomial)
LL_accuracy=summary(LL_accuracy_model) 

# probe analysis Rest (all besides HH) with accuracy category
rest_accuracy_model=glmer(Outcome ~ 1 + accuracy_category + (1|subjectID),data=subset(probe_data,(probe_data$PairType3 %in% c('HM','LM','LL') & !is.na(probe_data$accuracy_category))),na.action=na.omit,family=binomial)
rest_accuracy=summary(rest_accuracy_model) 

# interaction between memory regressors and the value level (HH vs. the rest)
intercept_and_value_model=glmer(Outcome ~ 1 + accuracy_category + isHH + (1 + accuracy_category|subjectID),data=subset(probe_data,(probe_data$PairType<=2 | probe_data$PairType>=6) & !is.na(probe_data$accuracy_category)),na.action=na.omit,family=binomial)
interaction_accuracy_model=glmer(Outcome ~ 1 + accuracy_category + isHH + accuracy_category:isHH + (1 + accuracy_category|subjectID),data=subset(probe_data,(probe_data$PairType<=2 | probe_data$PairType>=6) & !is.na(probe_data$accuracy_category)),na.action=na.omit,family=binomial)
interaction_accuracy=summary(interaction_accuracy_model) 
interaction_accuracy_joint=anova(intercept_and_value_model, interaction_accuracy_model)

# RT
# probe analysis HH with regressor RT difference NoGo-Go only for pairs where both items were remembered
HH_RT=summary(glmer(Outcome ~ 1 + RT_old_nogo_minus_go + (1+RT_old_nogo_minus_go|subjectID),data=subset(probe_only_remembered_items,(probe_only_remembered_items$PairType3=='HH')),na.action=na.omit,family=binomial)) 

# probe analysis HM with regressor RT difference NoGo-Go only for pairs where both items were remembered
HM_RT=summary(glmer(Outcome ~ 1 + RT_old_nogo_minus_go + (1+RT_old_nogo_minus_go|subjectID),data=subset(probe_only_remembered_items,(probe_only_remembered_items$PairType3=='HM')),na.action=na.omit,family=binomial)) 

# probe analysis LM with regressor RT difference NoGo-Go only for pairs where both items were remembered
LM_RT=summary(glmer(Outcome ~ 1 + RT_old_nogo_minus_go + (1+RT_old_nogo_minus_go|subjectID),data=subset(probe_only_remembered_items,(probe_only_remembered_items$PairType3=='LM')),na.action=na.omit,family=binomial)) 

# probe analysis LL with regressor RT difference NoGo-Go only for pairs where both items were remembered
LL_RT=summary(glmer(Outcome ~ 1 + RT_old_nogo_minus_go + (1+RT_old_nogo_minus_go|subjectID),data=subset(probe_only_remembered_items,(probe_only_remembered_items$PairType3=='LL')),na.action=na.omit,family=binomial)) 

# probe analysis Rest (all besides HH) with regressor RT difference NoGo-Go only for pairs where both items were remembered
rest_RT=summary(glmer(Outcome ~ 1 + RT_old_nogo_minus_go + (1+RT_old_nogo_minus_go|subjectID),data=subset(probe_only_remembered_items,(probe_only_remembered_items$PairType3 %in% c('HM','LM','LL'))),na.action=na.omit,family=binomial)) 

# probe analysis with regressors value (HH vs. the rest), RT difference NoGo-Go and the interaction, only for pairs where both items were remembered
interaction_RT=summary(glmer(Outcome ~ 1 + isHH * RT_old_nogo_minus_go + (1+RT_old_nogo_minus_go+isHH|subjectID),data=subset(probe_only_remembered_items,(probe_only_remembered_items$PairType<=2 | probe_only_remembered_items$PairType>=6)),na.action=na.omit,family=binomial)) 

# create a table with all results - HH vs. rest
df_results=data.frame(row.names=1:3)
df_results$value_level=c("HH", "The rest", "Interaction with value level (HH vs. the rest")

# accuracy Go remembered NoGo forgotten
accuracy_Go_remembered_estimate=c(HH_accuracy$coefficients[3,1], rest_accuracy$coefficients[3,1], interaction_accuracy$coefficients[6,1])
accuracy_Go_remembered_odds_ratio=round(exp(accuracy_Go_remembered_estimate),3)
accuracy_Go_remembered_se=c(HH_accuracy$coefficients[3,2], rest_accuracy$coefficients[3,2], interaction_accuracy$coefficients[6,2])
accuracy_Go_remembered_CI_min=round(exp(accuracy_Go_remembered_estimate-(1.96*accuracy_Go_remembered_se)),3)
accuracy_Go_remembered_CI_max=round(exp(accuracy_Go_remembered_estimate+(1.96*accuracy_Go_remembered_se)),3)
accuracy_Go_remembered_one_sided_p=round(c(HH_accuracy$coefficients[3,4], rest_accuracy$coefficients[3,4], interaction_accuracy$coefficients[6,4])/2,3)
df_results$accuracy_Go_remembered_P=accuracy_Go_remembered_one_sided_p
df_results$accuracy_Go_remembered_odds_ratio=accuracy_Go_remembered_odds_ratio
df_results$accuracy_Go_remembered_CI = paste(accuracy_Go_remembered_CI_min, '-',accuracy_Go_remembered_CI_max,sep="")

# accuracy Go forgotten NoGo remembered
accuracy_Go_forgotten_estimate=c(HH_accuracy$coefficients[2,1], rest_accuracy$coefficients[2,1],interaction_accuracy$coefficients[5,1])
accuracy_Go_forgotten_odds_ratio=round(exp(accuracy_Go_forgotten_estimate),3)
accuracy_Go_forgotten_se=c(HH_accuracy$coefficients[2,2], rest_accuracy$coefficients[2,2],interaction_accuracy$coefficients[5,2])
accuracy_Go_forgotten_CI_min=round(exp(accuracy_Go_forgotten_estimate-(1.96*accuracy_Go_forgotten_se)),3)
accuracy_Go_forgotten_CI_max=round(exp(accuracy_Go_forgotten_estimate+(1.96*accuracy_Go_forgotten_se)),3)
accuracy_Go_forgotten_one_sided_p=round(c(HH_accuracy$coefficients[2,4], rest_accuracy$coefficients[2,4], interaction_accuracy$coefficients[5,4])/2,3)
df_results$accuracy_Go_forgotten_P=accuracy_Go_forgotten_one_sided_p
df_results$accuracy_Go_forgotten_odds_ratio=accuracy_Go_forgotten_odds_ratio
df_results$accuracy_Go_forgotten_CI = paste(accuracy_Go_forgotten_CI_min, '-',accuracy_Go_forgotten_CI_max,sep="")

# RT
RT_NoGo_minus_Go_estimate=c(HH_RT$coefficients[2,1], rest_RT$coefficients[2,1], interaction_RT$coefficients[4,1])
RT_NoGo_minus_Go_odds_ratio=round(exp(RT_NoGo_minus_Go_estimate),3)
RT_NoGo_minus_Go_se=c(HH_RT$coefficients[2,2], rest_RT$coefficients[2,2], interaction_RT$coefficients[4,2])
RT_NoGo_minus_Go_CI_min=round(exp(RT_NoGo_minus_Go_estimate-(1.96*RT_NoGo_minus_Go_se)),3)
RT_NoGo_minus_Go_CI_max=round(exp(RT_NoGo_minus_Go_estimate+(1.96*RT_NoGo_minus_Go_se)),3)
RT_NoGo_minus_Go_one_sided_p=round(c(HH_RT$coefficients[2,4], rest_RT$coefficients[2,4], interaction_RT$coefficients[4,4])/2,3)
df_results$RT_NoGo_minus_Go_P=RT_NoGo_minus_Go_one_sided_p
df_results$RT_NoGo_minus_Go_odds_ratio=RT_NoGo_minus_Go_odds_ratio
df_results$RT_NoGo_minus_Go_CI = paste(RT_NoGo_minus_Go_CI_min, '-',RT_NoGo_minus_Go_CI_max,sep="")

print(df_results)

# create a table with all results - HM, LM, LL
df_results_rest=data.frame(row.names=1:3)
df_results_rest$value_level=c("HM", "LM", "LL")

# accuracy Go remembered NoGo forgotten
accuracy_Go_remembered_estimate=c(HM_accuracy$coefficients[3,1], LM_accuracy$coefficients[3,1], LL_accuracy$coefficients[3,1])
accuracy_Go_remembered_odds_ratio=round(exp(accuracy_Go_remembered_estimate),3)
accuracy_Go_remembered_se=c(HM_accuracy$coefficients[3,2], LM_accuracy$coefficients[3,2], LL_accuracy$coefficients[3,2])
accuracy_Go_remembered_CI_min=round(exp(accuracy_Go_remembered_estimate-(1.96*accuracy_Go_remembered_se)),3)
accuracy_Go_remembered_CI_max=round(exp(accuracy_Go_remembered_estimate+(1.96*accuracy_Go_remembered_se)),3)
accuracy_Go_remembered_one_sided_p=round(c(HM_accuracy$coefficients[3,4], LM_accuracy$coefficients[3,4], LL_accuracy$coefficients[3,4])/2,3)
df_results_rest$accuracy_Go_remembered_P = accuracy_Go_remembered_one_sided_p
df_results_rest$accuracy_Go_remembered_odds_ratio = accuracy_Go_remembered_odds_ratio
df_results_rest$accuracy_Go_remembered_CI = paste(accuracy_Go_remembered_CI_min, '-',accuracy_Go_remembered_CI_max,sep="")

# accuracy Go forgotten NoGo remembered
accuracy_Go_forgotten_estimate=c(HM_accuracy$coefficients[2,1],LM_accuracy$coefficients[2,1],LL_accuracy$coefficients[2,1])
accuracy_Go_forgotten_odds_ratio=round(exp(accuracy_Go_forgotten_estimate),3)
accuracy_Go_forgotten_se=c(HM_accuracy$coefficients[2,2],LM_accuracy$coefficients[2,2],LL_accuracy$coefficients[2,2])
accuracy_Go_forgotten_CI_min=round(exp(accuracy_Go_forgotten_estimate-(1.96*accuracy_Go_forgotten_se)),3)
accuracy_Go_forgotten_CI_max=round(exp(accuracy_Go_forgotten_estimate+(1.96*accuracy_Go_forgotten_se)),3)
accuracy_Go_forgotten_one_sided_p=round(c(HM_accuracy$coefficients[2,4], LM_accuracy$coefficients[2,4], LL_accuracy$coefficients[2,4])/2,3)
df_results_rest$accuracy_Go_forgotten_P = accuracy_Go_forgotten_one_sided_p
df_results_rest$accuracy_Go_forgotten_odds_ratio = accuracy_Go_forgotten_odds_ratio
df_results_rest$accuracy_Go_forgotten_CI = paste(accuracy_Go_forgotten_CI_min, '-',accuracy_Go_forgotten_CI_max,sep="")

# RT
RT_NoGo_minus_Go_estimate=c(HM_RT$coefficients[2,1], LM_RT$coefficients[2,1], LL_RT$coefficients[2,1])
RT_NoGo_minus_Go_odds_ratio=round(exp(RT_NoGo_minus_Go_estimate),3)
RT_NoGo_minus_Go_se=c(HM_RT$coefficients[2,2], LM_RT$coefficients[2,2], LL_RT$coefficients[2,2])
RT_NoGo_minus_Go_CI_min=round(exp(RT_NoGo_minus_Go_estimate-(1.96*RT_NoGo_minus_Go_se)),3)
RT_NoGo_minus_Go_CI_max=round(exp(RT_NoGo_minus_Go_estimate+(1.96*RT_NoGo_minus_Go_se)),3)
RT_NoGo_minus_Go_one_sided_p=round(c(HM_RT$coefficients[2,4], LM_RT$coefficients[2,4], LL_RT$coefficients[2,4])/2,3)
df_results_rest$RT_NoGo_minus_Go_P = RT_NoGo_minus_Go_one_sided_p
df_results_rest$RT_NoGo_minus_Go_odds_ratio = RT_NoGo_minus_Go_odds_ratio
df_results_rest$RT_NoGo_minus_Go_CI = paste(RT_NoGo_minus_Go_CI_min, '-',RT_NoGo_minus_Go_CI_max,sep="")
  
print(df_results_rest)

# visualization
# means
mean_props_combined=c(mean(with(data=subset(probe_data,probe_data$accuracy_category=="Both remembered / forgotten" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HH both remembered or both forgotten
                      mean(with(data=subset(probe_data,probe_data$accuracy_category=="Go remembered, NoGo forgotten" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HH Go remembered NoGo forgotten
                      mean(with(data=subset(probe_data,probe_data$accuracy_category=="Go forgotten, NoGo remembered" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HH Go forgotten NoGo remembered
                      mean(with(data=subset(probe_data,probe_data$accuracy_category=="Both remembered / forgotten" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # rest both remembered or both forgotten
                      mean(with(data=subset(probe_data,probe_data$accuracy_category=="Go remembered, NoGo forgotten" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # rest Go remembered NoGo forgotten
                      mean(with(data=subset(probe_data,probe_data$accuracy_category=="Go forgotten, NoGo remembered" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T)) # rest Go forgotten NoGo remembered

# SD values
SD_props_combined=c(sd(with(data=subset(probe_data,probe_data$accuracy_category=="Both remembered / forgotten" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HH both remembered or both forgotten
                    sd(with(data=subset(probe_data,probe_data$accuracy_category=="Go remembered, NoGo forgotten" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HH Go remembered NoGo forgotten
                    sd(with(data=subset(probe_data,probe_data$accuracy_category=="Go forgotten, NoGo remembered" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # HH Go forgotten NoGo remembered
                    sd(with(data=subset(probe_data,probe_data$accuracy_category=="Both remembered / forgotten" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # rest both remembered or both forgotten
                    sd(with(data=subset(probe_data,probe_data$accuracy_category=="Go remembered, NoGo forgotten" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T), # rest Go remembered NoGo forgotten
                    sd(with(data=subset(probe_data,probe_data$accuracy_category=="Go forgotten, NoGo remembered" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)), na.rm=T)) # rest Go forgotten NoGo remembered

# SE values
se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }
SE_props_combined=c(se(with(data=subset(probe_data,probe_data$accuracy_category=="Both remembered / forgotten" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T))), # HH both remembered or both forgotten
                    se(with(data=subset(probe_data,probe_data$accuracy_category=="Go remembered, NoGo forgotten" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T))), # HH Go remembered NoGo forgotten
                    se(with(data=subset(probe_data,probe_data$accuracy_category=="Go forgotten, NoGo remembered" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T))), # HH Go forgotten NoGo remembered
                    se(with(data=subset(probe_data,probe_data$accuracy_category=="Both remembered / forgotten" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T))), # rest both remembered or both forgotten
                    se(with(data=subset(probe_data,probe_data$accuracy_category=="Go remembered, NoGo forgotten" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T))), # rest Go remembered NoGo forgotten
                    se(with(data=subset(probe_data,probe_data$accuracy_category=="Go forgotten, NoGo remembered" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, mean, na.rm=T)))) # rest Go forgotten NoGo remembered

num_relevant_trials_combined = c(sum(with(data=subset(probe_data,probe_data$accuracy_category=="Both remembered / forgotten" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T), # HH both remembered or both forgotten
                                 sum(with(data=subset(probe_data,probe_data$accuracy_category=="Go remembered, NoGo forgotten" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T), # HH Go remembered NoGo forgotten
                                 sum(with(data=subset(probe_data,probe_data$accuracy_category=="Go forgotten, NoGo remembered" & probe_data$PairType3=="HH" & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T), # HH Go forgotten NoGo remembered
                                 sum(with(data=subset(probe_data,probe_data$accuracy_category=="Both remembered / forgotten" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T), # rest both remembered or both forgotten
                                 sum(with(data=subset(probe_data,probe_data$accuracy_category=="Go remembered, NoGo forgotten" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T), # rest Go remembered NoGo forgotten
                                 sum(with(data=subset(probe_data,probe_data$accuracy_category=="Go forgotten, NoGo remembered" & probe_data$PairType3 %in% c("HM","LM","LL") & !is.na(probe_data$Outcome)), tapply(Outcome, subjectID, length)), na.rm=T)) # rest Go forgotten NoGo remembered

# put values in data frame for plot
df_accuracy=data.frame(row.names = 1:6)
df_accuracy$value_level[1:3] = "High-value"
df_accuracy$value_level[4:6] = "The rest"
df_accuracy$means = mean_props_combined*100
df_accuracy$SEs = SE_props_combined*100
df_accuracy$SDs = SD_props_combined*100
df_accuracy$category = factor(c("Both remembered / forgotten", "Go remembered, NoGo forgotten", "Go forgotten, NoGo remembered"), levels = c("Go remembered, NoGo forgotten", "Both remembered / forgotten", "Go forgotten, NoGo remembered"))
df_accuracy$considered_trials = num_relevant_trials_combined

df_accuracy$pval = NA
df_accuracy$pval[df_accuracy$category=="Go remembered, NoGo forgotten" & df_accuracy$value_level=="High-value"] = df_results$accuracy_Go_remembered_P[df_results$value_level=="HH"]
df_accuracy$pval[df_accuracy$category=="Go remembered, NoGo forgotten" & df_accuracy$value_level=="The rest"] = df_results$accuracy_Go_remembered_P[df_results$value_level=="The rest"]
df_accuracy$pval[df_accuracy$category=="Go forgotten, NoGo remembered" & df_accuracy$value_level=="High-value"] = df_results$accuracy_Go_forgotten_P[df_results$value_level=="HH"]
df_accuracy$pval[df_accuracy$category=="Go forgotten, NoGo remembered" & df_accuracy$value_level=="The rest"] = df_results$accuracy_Go_forgotten_P[df_results$value_level=="The rest"]
df_accuracy$asterisk=""
df_accuracy$asterisk[df_accuracy$pval<0.07]="+"
df_accuracy$asterisk[df_accuracy$pval<0.05]="*"
df_accuracy$asterisk[df_accuracy$pval<0.01]="**"
df_accuracy$asterisk[df_accuracy$pval<0.001]="***"

plot_accuracy= ggplot(data=df_accuracy, aes(x=value_level, y=means, fill=category)) +
  geom_bar(stat="identity", position=position_dodge(0.7), width=0.7) +
  theme_bw() + # white background
  theme(axis.title.x=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank()) + # axis and background formating
  theme(aspect.ratio=1.2, text = element_text(size=32)) + # font size
  theme(legend.position="bottom", legend.direction = "vertical", legend.title=element_blank()) + # position legend
  geom_text(aes(label = paste(considered_trials), y=8), size=12, position = position_dodge(0.7), angle = 90) +
  geom_text(aes(label = paste(asterisk), y = means + 8), size=12, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin=means-SEs, ymax=means+SEs), width=1/4, position=position_dodge(0.7)) +
  scale_y_continuous("Mean proportion of Go item choices", expand = c(0,0), limit=c(0,100), breaks=seq(0,100,5)) + # define y axis properties
  scale_fill_manual(limits = c("Go remembered, NoGo forgotten", "Both remembered / forgotten", "Go forgotten, NoGo remembered"), values=c("#585858","#909090","#DCDCDC")) # color of bars

# show the plot
print(plot_accuracy)

# plot RT
probe_only_remembered_items_go_nogo = subset(probe_only_remembered_items, probe_only_remembered_items$PairType <= 2 | probe_only_remembered_items$PairType >= 6)
RT_data_agg=as.data.frame(aggregate(probe_only_remembered_items_go_nogo,by=list(probe_only_remembered_items_go_nogo$subjectID_num, probe_only_remembered_items_go_nogo$isHH, probe_only_remembered_items_go_nogo$Outcome), mean, na.rm=TRUE))
Data_by_sub=melt(RT_data_agg,id = c("subjectID_num", "isHH", "Outcome") ,"RT_old_nogo_minus_go")
RT_means_for_plot=with(data=Data_by_sub, tapply(value, list(isHH, Outcome), mean))
RT_SEs_for_plot=with(data=Data_by_sub, tapply(value, list(isHH, Outcome), se))

df_RT_for_plot=data.frame(row.names = 1:4)
df_RT_for_plot$value_level = factor(c("High-value", "High-value", "The rest", "The rest"), levels=c("High-value", "The rest"))
df_RT_for_plot$item_type=c("Go chosen", "NoGo chosen", "Go chosen", "NoGo chosen")
df_RT_for_plot$means = c(RT_means_for_plot[2,2], RT_means_for_plot[2,1], RT_means_for_plot[1,2], RT_means_for_plot[1,1])
df_RT_for_plot$SEs = c(RT_SEs_for_plot[2,2], RT_SEs_for_plot[2,1], RT_SEs_for_plot[1,2], RT_SEs_for_plot[1,1])

df_RT_for_plot$pval = NA
df_RT_for_plot$pval[df_RT_for_plot$value_level=="High-value" & df_RT_for_plot$item_type == "Go chosen"] = df_results$RT_NoGo_minus_Go_P[df_results$value_level=="HH"]
df_RT_for_plot$pval[df_RT_for_plot$value_level=="The rest" & df_RT_for_plot$item_type == "Go chosen"] = df_results$RT_NoGo_minus_Go_P[df_results$value_level=="The rest"]
df_RT_for_plot$asterisk=""
df_RT_for_plot$asterisk[df_RT_for_plot$pval<0.07]="+"
df_RT_for_plot$asterisk[df_RT_for_plot$pval<0.05]="*"
df_RT_for_plot$asterisk[df_RT_for_plot$pval<0.01]="**"
df_RT_for_plot$asterisk[df_RT_for_plot$pval<0.001]="***"

# plot RT-choices
plot_RT= ggplot(data=df_RT_for_plot, aes(x=value_level, y=means, fill=item_type)) +
  geom_bar(stat="identity", position=position_dodge(0.7), width=0.7) +
  geom_text(aes(label = paste(asterisk), y = means + 0.08), size=12) +
  theme_bw() + # white background
  theme(legend.position="bottom",legend.title=element_blank()) + # position legend
  theme(axis.title.x=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank()) + # axis and background formating
  theme(aspect.ratio=1.2, text = element_text(size=32)) + # font size
  geom_errorbar(aes(ymin=means-SEs, ymax=means+SEs), width=1/4, position=position_dodge(0.7)) +
  scale_y_continuous(expression(paste(Delta, "RT (NoGo recognition RT - Go recognition RT)")), expand = c(0,0), limit=c(-0.1,0.28), breaks=seq(-0.1,0.28,0.02)) + # define y axis properties
  scale_fill_manual(values=c("#585858","#DCDCDC")) + # color of bars
  geom_abline(intercept = 0,slope=0,linetype =2, size = 1,aes()) # add line for y=0

# show the plot
print(plot_RT)