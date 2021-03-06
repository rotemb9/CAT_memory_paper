# Required R package
library(lme4)
library(lmerTest)
library(tidyverse)

# Clear workspace
rm(list=ls())

# load recognition data for all experiments
# Define the local path where the data can be found
# you will need to set the correct path whre the RData file with the data is located
input_path="DEFINE PATH OF DATA HERE"

recognition_data_all_experiments=c()
experiments_names = c("bmem_snacks", "bmem_snacks2","bmem_short")
for (exp_ind in 1:length(experiments_names)){
  experiment_name=experiments_names[exp_ind]
  if (experiment_name=="bmem_short"){
    input_filename=paste(input_path, experiment_name, "_recognition.Rdata",sep="")
    load(file=input_filename) 
    rm(recognition_probe_items, recognition_old_items)
    recognition_data$experiment=experiment_name
    recognition_data$session=NA
    recognition_data_all_experiments=rbind(recognition_data_all_experiments,recognition_data)
  } else {
    input_filename=paste(input_path, experiment_name, "_recognition_all_sessions.Rdata",sep="")
    load(file=input_filename) 
    recognition_data_all_sessions$experiment=experiment_name
    recognition_data_all_sessions$PairType=NA
    recognition_data_all_experiments=rbind(recognition_data_all_experiments,recognition_data_all_sessions)
  }

}

recognition_data=recognition_data_all_experiments
# subset data
recognition_probe_items=subset(recognition_data,recognition_data$isProbeItem)
rm(recognition_data_all_experiments, recognition_data_all_sessions)

num_exp_sessions = 2*length(experiments_names)
if ("bmem_short" %in% experiments_names){
  num_exp_sessions = num_exp_sessions - 1 # only one session for bmem_short
}

recognition_gonogo_descriptive_results = data.frame(row.names = 1:num_exp_sessions)
recognition_confidence_results = data.frame(row.names = 1:num_exp_sessions)
ind=0

# loop through experiments
for (exp_ind in 1:length(experiments_names)) {
  experiment_name=experiments_names[exp_ind]
  if (experiment_name!="bmem_short"){ # there are several sessions
    sessions = c('Session 2','Follow-up')
  } else {
    sessions = 1
  }
  for (curr_session in sessions) {
      ind = ind+1
      if (experiment_name!="bmem_short"){
        recognition_data_curr_session = subset(recognition_data,(recognition_data$experiment==experiment_name & recognition_data$session==curr_session))
        recognition_probe_items_curr_session = subset(recognition_probe_items, (recognition_probe_items$experiment==experiment_name & recognition_probe_items$session==curr_session))
      } else {
        recognition_data_curr_session = subset(recognition_data,(recognition_data$experiment==experiment_name))
        recognition_probe_items_curr_session = subset(recognition_probe_items, (recognition_probe_items$experiment==experiment_name))
      }
      recognition_gonogo_descriptive_results$experiment[ind]=experiment_name
      recognition_gonogo_descriptive_results$session[ind]=curr_session
      recognition_confidence_results$experiment[ind]=experiment_name
      recognition_confidence_results$session[ind]=curr_session

    # Go / NoGo recognition task
    # hit rate
    gonogo_hit_rate_all_participants=with(data=subset(recognition_data_curr_session,isGo.==1), tapply(IsCorrectAnsGo, subjectID, mean, na.rm=T))
    gonogo_hit_rate_mean=mean(gonogo_hit_rate_all_participants,na.rm=T)
    gonogo_hit_rate_sd=sd(gonogo_hit_rate_all_participants,na.rm=T)
    recognition_gonogo_descriptive_results$hit_rate[ind]=paste(round(gonogo_hit_rate_mean*100,2), " (", round(gonogo_hit_rate_sd*100,2), ")", sep="")
    # correct rejection
    gonogo_cr_all_participants=with(data=subset(recognition_data_curr_session,isGo.==0), tapply(IsCorrectAnsGo, subjectID, mean, na.rm=T))
    gonogo_cr_rate_mean=mean(gonogo_cr_all_participants,na.rm=T)
    gonogo_cr_rate_sd=sd(gonogo_cr_all_participants,na.rm=T)
    recognition_gonogo_descriptive_results$correct_rejection_rate[ind]=paste(round(gonogo_cr_rate_mean*100,2), " (", round(gonogo_cr_rate_sd*100,2),")", sep="")
    # d prime
    gonogo_false_alarm_all_participants = 1-gonogo_cr_all_participants
    # d prime cannot be computed when hit rate / false alarm = 1 or 0. Adjust for this:
    gonogo_hit_rate_all_participants[gonogo_hit_rate_all_participants==1] = 0.99
    gonogo_hit_rate_all_participants[gonogo_hit_rate_all_participants==0] = 0.01
    gonogo_false_alarm_all_participants[gonogo_false_alarm_all_participants==0] = 0.01
    gonogo_false_alarm_all_participants[gonogo_false_alarm_all_participants==1] = 0.99
    gonogo_d_prime_all_participants=qnorm(gonogo_hit_rate_all_participants) - qnorm(gonogo_false_alarm_all_participants)
    gonogo_d_prime_mean=mean(gonogo_d_prime_all_participants,na.rm=T)
    gonogo_d_prime_sd=sd(gonogo_d_prime_all_participants,na.rm=T)
    recognition_gonogo_descriptive_results$d_prime[ind]=paste(round(gonogo_d_prime_mean,3), " (", round(gonogo_d_prime_sd,3),")", sep="")
    
    # confidence levels
    recognition_probe_items_curr_session_correct_isOld= subset(recognition_probe_items_curr_session,recognition_probe_items_curr_session$IsCorrectAnsOld==1)
    recognition_probe_items_curr_session_correct_isOld$is_high_confidence_level=recognition_probe_items_curr_session_correct_isOld$subjectAnswerIsOld==1
    
    # descriptive: confidence level (True- high; false- low) for go vs nogo items (correct isold answers):
    #confidence_isgo=tapply(recognition_probe_items_curr_session_correct_isOld$subjectAnswerIsOld,recognition_probe_items_curr_session_correct_isOld$isGo.,mean,na.rm=T)
    confidence_go_all_participants=with(data=subset(recognition_probe_items_curr_session_correct_isOld, recognition_probe_items_curr_session_correct_isOld$isGo.), tapply(is_high_confidence_level, subjectID, mean, na.rm=T))
    confidence_nogo_all_participants=with(data=subset(recognition_probe_items_curr_session_correct_isOld, !recognition_probe_items_curr_session_correct_isOld$isGo.), tapply(is_high_confidence_level, subjectID, mean, na.rm=T))
    confidence_go_mean = mean(confidence_go_all_participants, na.rm=T)
    confidence_go_sd = sd(confidence_go_all_participants, na.rm=T)
    confidence_nogo_mean = mean(confidence_nogo_all_participants, na.rm=T)
    confidence_nogo_sd = sd(confidence_nogo_all_participants, na.rm=T)
    
    # confidence level go vs. nogo for correct is old answers (1=high confidence 0=low confidence)
    confidence_analysis=summary(glmer(is_high_confidence_level ~ 1 + isGo.+IsHighValue + (1|subjectID),data=recognition_probe_items_curr_session_correct_isOld,na.action=na.omit,family=binomial)) 
    recognition_probe_items_curr_session_correct_isOld$go.ind = 1*recognition_probe_items_curr_session_correct_isOld$isGo.
    recognition_probe_items_curr_session_correct_isOld$high.ind = 1*recognition_probe_items_curr_session_correct_isOld$IsHighValue
    if (experiment_name == "bmem_snacks" & curr_session == "Session 2") {
      confidence_analysis=summary(glmer(is_high_confidence_level ~ 1 + isGo.+IsHighValue + (1 + go.ind|subjectID),data=recognition_probe_items_curr_session_correct_isOld,na.action=na.omit,family=binomial)) 
    }
    if (experiment_name == "bmem_snacks" & curr_session == "Follow-up") {
      confidence_analysis=summary(glmer(is_high_confidence_level ~ 1 + isGo.+IsHighValue + (1 + high.ind|subjectID),data=recognition_probe_items_curr_session_correct_isOld,na.action=na.omit,family=binomial)) 
    }
    
    if (experiment_name == "bmem_snacks2" & curr_session == "Session 2") {
      confidence_analysis=summary(glmer(is_high_confidence_level ~ 1 + isGo.+IsHighValue + (1 + go.ind |subjectID),data=recognition_probe_items_curr_session_correct_isOld,na.action=na.omit,family=binomial)) 
    }  
    
    if (experiment_name == "bmem_snacks2" & curr_session == "Follow-up") {
      confidence_analysis=summary(glmer(is_high_confidence_level ~ 1 + isGo.+IsHighValue + (1 + go.ind + high.ind||subjectID),data=recognition_probe_items_curr_session_correct_isOld,na.action=na.omit,family=binomial)) 
    }     
    if (experiment_name == "bmem_short" & curr_session == 1) {
      confidence_analysis=summary(glmer(is_high_confidence_level ~ 1 + isGo.+IsHighValue + (1 + go.ind + high.ind||subjectID),data=recognition_probe_items_curr_session_correct_isOld,na.action=na.omit,family=binomial)) 
    }      
    recognition_confidence_results$confidence_Go[ind]=paste(round(confidence_go_mean,3)," (", round(confidence_go_sd,3),")",sep="")
    recognition_confidence_results$confidence_NoGo[ind]=paste(round(confidence_nogo_mean,3)," (", round(confidence_nogo_sd,3),")",sep="")
    confidence_CI_min=round(exp(confidence_analysis$coefficients[2,1]-1.96*confidence_analysis$coefficients[2,2]),3)
    confidence_CI_max=round(exp(confidence_analysis$coefficients[2,1]+1.96*confidence_analysis$coefficients[2,2]),3)
    recognition_confidence_results$P[ind]=round(confidence_analysis$coefficients[2,4],3)
    recognition_confidence_results$odds_ratio[ind]=round(exp(confidence_analysis$coefficients[2,1]),3)
    recognition_confidence_results$CI[ind]=paste(confidence_CI_min, "-", confidence_CI_max, sep="")
    }
}
print("Go / NoGo recognition task")
print(recognition_gonogo_descriptive_results)
print("Confidence level analysis")
print(recognition_confidence_results)

# Analysis of RT - correct vs. incorrect responses, Go vs. NoGo
# Only for Experiment 1 (bmem_snacks2), where the recognition RT effect was found
recognition_data_exp1=subset(recognition_data,recognition_data$experiment %in% c("bmem_snacks2"))
recognition_probe_items_exp1 = subset(recognition_probe_items,recognition_probe_items$experiment %in% c("bmem_snacks2"))

# RT old/new Go vs. NoGo for hits vs. misses
# means
rt_hits_go=mean(with(data=subset(recognition_probe_items_exp1,isGo.==1 & IsCorrectAnsOld==1), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T)
rt_hits_nogo=mean(with(data=subset(recognition_probe_items_exp1,isGo.==0 & IsCorrectAnsOld==1), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T)
rt_misses_go=mean(with(data=subset(recognition_probe_items_exp1,isGo.==1 & IsCorrectAnsOld==0), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T)
rt_misses_nogo=mean(with(data=subset(recognition_probe_items_exp1,isGo.==0 & IsCorrectAnsOld==0), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T)
# SDs
rt_hits_go_sd=sd(with(data=subset(recognition_probe_items_exp1,isGo.==1 & IsCorrectAnsOld==1), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T)
rt_hits_nogo_sd=sd(with(data=subset(recognition_probe_items_exp1,isGo.==0 & IsCorrectAnsOld==1), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T)
rt_misses_go_sd=sd(with(data=subset(recognition_probe_items_exp1,isGo.==1 & IsCorrectAnsOld==0), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T)
rt_misses_nogo_sd=sd(with(data=subset(recognition_probe_items_exp1,isGo.==0 & IsCorrectAnsOld==0), tapply(RT_isOld, subjectID, mean, na.rm=T)),na.rm=T)
# linear regression
RT_diff_hits=summary(lmer(RT_isOld ~ isGo. + IsHighValue + (1+isGo.|subjectID),data=subset(recognition_probe_items_exp1,IsCorrectAnsOld==1),na.action=na.omit))
RT_diff_misses=summary(lmer(RT_isOld ~ isGo. + IsHighValue + (1+isGo.|subjectID),data=subset(recognition_probe_items_exp1,IsCorrectAnsOld==0),na.action=na.omit))
interaction_RTs=summary(lmer(RT_isOld ~ isGo. * IsCorrectAnsOld + IsHighValue + (1 + isGo.|subjectID),data=recognition_probe_items_exp1,na.action=na.omit))
# num trials
num_hits_go_items=length(recognition_probe_items_exp1$RT_isOld[recognition_probe_items_exp1$IsCorrectAnsOld==1 & recognition_probe_items_exp1$isGo.==1])
num_hits_nogo_items=length(recognition_probe_items_exp1$RT_isOld[recognition_probe_items_exp1$IsCorrectAnsOld==1 & recognition_probe_items_exp1$isGo.==0])
num_misses_go_items=length(recognition_probe_items_exp1$RT_isOld[recognition_probe_items_exp1$IsCorrectAnsOld==0 & recognition_probe_items_exp1$isGo.==1])
num_misses_nogo_items=length(recognition_probe_items_exp1$RT_isOld[recognition_probe_items_exp1$IsCorrectAnsOld==0 & recognition_probe_items_exp1$isGo.==0])

# print results
print("RT Go vs. NoGo for hits:")
print(paste("RT Go hits = ",round(rt_hits_go,3), " +- ", round(rt_hits_go_sd,3), "   RT NoGo hits = ", round(rt_hits_nogo,3), " +- ", round(rt_hits_nogo_sd,3), sep=""))
print(paste("one-sided p = ",round(RT_diff_hits$coefficients[2,5]/2,3),sep=""))
print("RT Go vs. NoGo for misses:")
print(paste("RT Go misses = ",round(rt_misses_go,3), " +- ", round(rt_misses_go_sd,3), "   RT NoGo misses = ", round(rt_misses_nogo,3), " +- ", round(rt_misses_nogo_sd,3), sep=""))
print(paste("one-sided p = ",round(RT_diff_misses$coefficients[2,5]/2,3),sep=""))
print("Only for probe items")
print(paste("Number of hits go = ", num_hits_go_items, sep=""))
print(paste("Number of hits nogo = ", num_hits_nogo_items, sep=""))
print(paste("Number of misses go = ", num_misses_go_items, sep=""))
print(paste("Number of misses nogo = ", num_misses_nogo_items, sep=""))
print("Interaction Go/NoGo and hits/misses:")
print(paste("one-sided p = ",round(interaction_RTs$coefficients[5,5]/2,3)))

# RT go/nogo Go vs. NoGo for hits vs. misses
print("Go / NoGo recognition task!")
# means
rt_gonogo_correct_go=mean(with(data=subset(recognition_probe_items_exp1,isGo.==1 & IsCorrectAnsGo==1), tapply(RT_isGo, subjectID, mean, na.rm=T)),na.rm=T)
rt_gonogo_correct_nogo=mean(with(data=subset(recognition_probe_items_exp1,isGo.==0 & IsCorrectAnsGo==1), tapply(RT_isGo, subjectID, mean, na.rm=T)),na.rm=T)
rt_gonogo_incorrect_go=mean(with(data=subset(recognition_probe_items_exp1,isGo.==1 & IsCorrectAnsGo==0), tapply(RT_isGo, subjectID, mean, na.rm=T)),na.rm=T)
rt_gonogo_incorrect_nogo=mean(with(data=subset(recognition_probe_items_exp1,isGo.==0 & IsCorrectAnsGo==0), tapply(RT_isGo, subjectID, mean, na.rm=T)),na.rm=T)
# SDs
rt_gonogo_correct_go_sd=sd(with(data=subset(recognition_probe_items_exp1,isGo.==1 & IsCorrectAnsGo==1), tapply(RT_isGo, subjectID, mean, na.rm=T)),na.rm=T)
rt_gonogo_correct_nogo_sd=sd(with(data=subset(recognition_probe_items_exp1,isGo.==0 & IsCorrectAnsGo==1), tapply(RT_isGo, subjectID, mean, na.rm=T)),na.rm=T)
rt_gonogo_incorrect_go_sd=sd(with(data=subset(recognition_probe_items_exp1,isGo.==1 & IsCorrectAnsGo==0), tapply(RT_isGo, subjectID, mean, na.rm=T)),na.rm=T)
rt_gonogo_incorrect_nogo_sd=sd(with(data=subset(recognition_probe_items_exp1,isGo.==0 & IsCorrectAnsGo==0), tapply(RT_isGo, subjectID, mean, na.rm=T)),na.rm=T)
# logistic regression
RT_gonogo_diff_correct=summary(lmer(RT_isGo ~ isGo. + IsHighValue + (1+isGo.|subjectID),data=subset(recognition_probe_items_exp1,IsCorrectAnsGo==1),na.action=na.omit))
RT_gonogo_diff_incorrect=summary(lmer(RT_isGo ~ isGo. + IsHighValue + (1+IsHighValue|subjectID),data=subset(recognition_probe_items_exp1,IsCorrectAnsGo==0),na.action=na.omit))
interaction_RTs_gonogo=summary(lmer(RT_isGo ~ isGo. * IsCorrectAnsGo + IsHighValue + (1+IsHighValue|subjectID),data=recognition_probe_items_exp1,na.action=na.omit))
# num trials
num_correct_go_items_gonogo=length(recognition_probe_items_exp1$RT_isGo[recognition_probe_items_exp1$IsCorrectAnsGo==1 & recognition_probe_items_exp1$isGo.==1])
num_correct_nogo_items_gonogo=length(recognition_probe_items_exp1$RT_isGo[recognition_probe_items_exp1$IsCorrectAnsGo==1 & recognition_probe_items_exp1$isGo.==0])
num_incorrect_go_items_gonogo=length(recognition_probe_items_exp1$RT_isGo[recognition_probe_items_exp1$IsCorrectAnsGo==0 & recognition_probe_items_exp1$isGo.==1])
num_incorrect_nogo_items_gonogo=length(recognition_probe_items_exp1$RT_isGo[recognition_probe_items_exp1$IsCorrectAnsGo==0 & recognition_probe_items_exp1$isGo.==0])

# print results
print("RT Go vs. NoGo for correct responses:")
print(paste("RT Go correct = ",round(rt_gonogo_correct_go,3), " +- ", round(rt_gonogo_correct_go_sd,3), "   RT NoGo correct = ", round(rt_gonogo_correct_nogo,3), " +- ", round(rt_gonogo_correct_nogo_sd,3), sep=""))
print(paste("two-sided p = ",round(RT_gonogo_diff_correct$coefficients[2,5],3),sep=""))
print("RT Go vs. NoGo for incorrect responses:")
print(paste("RT Go incorrect = ",round(rt_gonogo_incorrect_go,3), " +- ", round(rt_gonogo_incorrect_go_sd,3), "   RT NoGo incorrect = ", round(rt_gonogo_incorrect_nogo,3), " +- ", round(rt_gonogo_incorrect_nogo_sd,3), sep=""))
print(paste("two-sided p = ",round(RT_gonogo_diff_incorrect$coefficients[2,5],3),sep=""))
print("Only for probe items")
print(paste("Number of correct responses gonogo go = ", num_correct_go_items_gonogo, sep=""))
print(paste("Number of correct responses gonogo nogo = ", num_correct_nogo_items_gonogo, sep=""))
print(paste("Number of incorrect responses gonogo go = ", num_incorrect_go_items_gonogo, sep=""))
print(paste("Number of incorrect responses gonogo nogo = ", num_incorrect_nogo_items_gonogo, sep=""))
print("Interaction Go/NoGo and correct/incorrect responses:")
print(paste("two-sided p = ",round(interaction_RTs_gonogo$coefficients[5,5],3)))


### Probe RT analysis
# read probe data from all experiments and sessions and combine them

probe_data_all_experiments=c()
data_files = c("bmem_snacks2_probe2_recognition2.Rdata",
               "bmem_snacks2_probe3_recognition3.Rdata",
               "bmem_short_probe_recognition.Rdata",
               "bmem_snacks_probe2_recognition2.Rdata",
               "bmem_snacks_probe3_recognition3.Rdata")
experiment_names = c("bmem_snacks2",
                     "bmem_snacks2",
                     "bmem_short",
                     "bmem_snacks",
                     "bmem_snacks")
sessions = c(2,3,1,2,3)

for (data_ind in 1:length(data_files)){
  full_filename = file.path(input_path,data_files[data_ind])
  load(file=full_filename)
  probe_data$experiment_name = experiment_names[data_ind]
  probe_data$session = sessions[data_ind]
  if (experiment_names[data_ind] != "bmem_short"){
    probe_data$PairType3 = factor(probe_data$PairType2, levels = c("Low_Value","High_Value"))
  }
  rm(recognition_data,recognition_probe_items, recognition_probe_items_correct_isGo, recognition_probe_items_correct_isOld,probe_only_remembered_items)
  probe_data_all_experiments=rbind(probe_data_all_experiments,probe_data)
}

rm(probe_data)

probe_data_all_experiments$Outcome = as.logical(probe_data_all_experiments$Outcome)
probe_data_all_experiments$Accuracy_old_chosen_item[!is.na(probe_data_all_experiments$Outcome) & probe_data_all_experiments$Outcome & probe_data_all_experiments$PairType2 %in% c("High_Value", "Low_Value")] = probe_data_all_experiments$Accuracy_old_goItem[!is.na(probe_data_all_experiments$Outcome) & probe_data_all_experiments$Outcome & probe_data_all_experiments$PairType2 %in% c("High_Value", "Low_Value")]
probe_data_all_experiments$Accuracy_old_chosen_item[!is.na(probe_data_all_experiments$Outcome) & !probe_data_all_experiments$Outcome & probe_data_all_experiments$PairType2 %in% c("High_Value", "Low_Value")] = probe_data_all_experiments$Accuracy_old_nogoItem[!is.na(probe_data_all_experiments$Outcome) & !probe_data_all_experiments$Outcome & probe_data_all_experiments$PairType2 %in% c("High_Value", "Low_Value")]
probe_data_all_experiments$Accuracy_old_unchosen_item[!is.na(probe_data_all_experiments$Outcome) & probe_data_all_experiments$Outcome & probe_data_all_experiments$PairType2 %in% c("High_Value", "Low_Value")] = probe_data_all_experiments$Accuracy_old_nogoItem[!is.na(probe_data_all_experiments$Outcome) & probe_data_all_experiments$Outcome & probe_data_all_experiments$PairType2 %in% c("High_Value", "Low_Value")]
probe_data_all_experiments$Accuracy_old_unchosen_item[!is.na(probe_data_all_experiments$Outcome) & !probe_data_all_experiments$Outcome & probe_data_all_experiments$PairType2 %in% c("High_Value", "Low_Value")] = probe_data_all_experiments$Accuracy_old_goItem[!is.na(probe_data_all_experiments$Outcome) & !probe_data_all_experiments$Outcome & probe_data_all_experiments$PairType2 %in% c("High_Value", "Low_Value")]

# model
probe_RT_HV = summary(lmer(RT ~ 1 + Outcome + Accuracy_old_chosen_item + (1 + Outcome + Accuracy_old_chosen_item|subjectID),data=subset(probe_data_all_experiments,probe_data_all_experiments$PairType2 == "High_Value"),na.action=na.omit)) 
probe_RT_LV = summary(lmer(RT ~ 1 + Outcome + Accuracy_old_chosen_item + (1 + Outcome + Accuracy_old_chosen_item|subjectID),data=subset(probe_data_all_experiments,probe_data_all_experiments$PairType2 == "Low_Value"),na.action=na.omit)) 

# Get the numbers
relevant_data = probe_data_all_experiments[probe_data_all_experiments$PairType2 %in% c("High_Value","Low_Value"),]

mean_RT_HV_chosen_rem = round(mean(with(data=subset(probe_data_all_experiments,PairType2=="High_Value" & Accuracy_old_chosen_item), tapply(RT, subjectID, mean, na.rm=T)),na.rm=T),3)
mean_RT_HV_chosen_for = round(mean(with(data=subset(probe_data_all_experiments,PairType2=="High_Value" & !Accuracy_old_chosen_item), tapply(RT, subjectID, mean, na.rm=T)),na.rm=T),3)
mean_RT_LV_chosen_rem = round(mean(with(data=subset(probe_data_all_experiments,PairType2=="Low_Value" & Accuracy_old_chosen_item), tapply(RT, subjectID, mean, na.rm=T)),na.rm=T),3)
mean_RT_LV_chosen_for = round(mean(with(data=subset(probe_data_all_experiments,PairType2=="Low_Value" & !Accuracy_old_chosen_item), tapply(RT, subjectID, mean, na.rm=T)),na.rm=T),3)

mean_RT_HV_chosen_go = round(mean(with(data=subset(probe_data_all_experiments,PairType2=="High_Value" & Outcome), tapply(RT, subjectID, mean, na.rm=T)),na.rm=T),3)
mean_RT_HV_chosen_nogo = round(mean(with(data=subset(probe_data_all_experiments,PairType2=="High_Value" & !Outcome), tapply(RT, subjectID, mean, na.rm=T)),na.rm=T),3)
mean_RT_LV_chosen_go = round(mean(with(data=subset(probe_data_all_experiments,PairType2=="Low_Value" & Outcome), tapply(RT, subjectID, mean, na.rm=T)),na.rm=T),3)
mean_RT_LV_chosen_nogo = round(mean(with(data=subset(probe_data_all_experiments,PairType2=="Low_Value" & !Outcome), tapply(RT, subjectID, mean, na.rm=T)),na.rm=T),3)

### probe: Go remembered, NoGo forgotten vs. Go forgotten, NoGo remembered
probe_data_all_experiments$accuracy_category = NA
probe_data_all_experiments$accuracy_category[(probe_data_all_experiments$Accuracy_old_goItem & probe_data_all_experiments$Accuracy_old_nogoItem) | (!probe_data_all_experiments$Accuracy_old_goItem & !probe_data_all_experiments$Accuracy_old_nogoItem)] = "Both remembered / forgotten"
probe_data_all_experiments$accuracy_category[probe_data_all_experiments$Accuracy_old_goItem & !probe_data_all_experiments$Accuracy_old_nogoItem] = "Go remembered, NoGo forgotten"
probe_data_all_experiments$accuracy_category[!probe_data_all_experiments$Accuracy_old_goItem & probe_data_all_experiments$Accuracy_old_nogoItem] = "Go forgotten, NoGo remembered"
probe_data_all_experiments$accuracy_category = relevel(factor(probe_data_all_experiments$accuracy_category), ref = "Both remembered / forgotten")
probe_data_all_experiments$accuracy_category2[probe_data_all_experiments$Accuracy_old_chosen_item == probe_data_all_experiments$Accuracy_old_unchosen_item] = "Both remembered / forgotten"
probe_data_all_experiments$accuracy_category2[probe_data_all_experiments$Accuracy_old_chosen_item & !probe_data_all_experiments$Accuracy_old_unchosen_item] = "Chosen remembered, unchosen forgotten"
probe_data_all_experiments$accuracy_category2[!probe_data_all_experiments$Accuracy_old_chosen_item & probe_data_all_experiments$Accuracy_old_unchosen_item] = "Chosen forgotten, unchosen remembered"
probe_data_all_experiments$accuracy_category2 = relevel(factor(probe_data_all_experiments$accuracy_category2), ref = "Both remembered / forgotten")
probe_data_all_experiments$accuracy_category_without_both = probe_data_all_experiments$accuracy_category
probe_data_all_experiments$accuracy_category_without_both[probe_data_all_experiments$accuracy_category_without_both == "Both remembered / forgotten"] = NA
probe_data_all_experiments$accuracy_category_without_both = droplevels(probe_data_all_experiments$accuracy_category_without_both)

probe_choices_HV_accuracy_category = summary(glmer(Outcome ~ 1 + accuracy_category_without_both + (1 + accuracy_category_without_both|subjectID) + (1|experiment_name),data=subset(probe_data_all_experiments,(probe_data_all_experiments$PairType2=="High_Value")),na.action=na.omit,family=binomial))
probe_choices_LV_accuracy_category = summary(glmer(Outcome ~ 1 + accuracy_category_without_both + (1 + accuracy_category_without_both|subjectID) + (1|experiment_name),data=subset(probe_data_all_experiments,(probe_data_all_experiments$PairType2=="Low_Value")),na.action=na.omit,family=binomial))

## get the statistics for the manuscript
# means
mean_HV_go_rem_nogo_for = mean(with(data=subset(probe_data_all_experiments,PairType2=="High_Value" & accuracy_category_without_both=="Go remembered, NoGo forgotten"), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T)
mean_HV_go_for_nogo_rem = mean(with(data=subset(probe_data_all_experiments,PairType2=="High_Value" & accuracy_category_without_both=="Go forgotten, NoGo remembered"), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T)
mean_LV_go_rem_nogo_for = mean(with(data=subset(probe_data_all_experiments,PairType2=="Low_Value" & accuracy_category_without_both=="Go remembered, NoGo forgotten"), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T)
mean_LV_go_for_nogo_rem = mean(with(data=subset(probe_data_all_experiments,PairType2=="Low_Value" & accuracy_category_without_both=="Go forgotten, NoGo remembered"), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T)

# OR
probe_choices_HV_accuracy_category_coef = as.data.frame(probe_choices_HV_accuracy_category$coefficients)
HV_OR = round(exp(probe_choices_HV_accuracy_category_coef$Estimate[2]),3)
probe_choices_LV_accuracy_category_coef = as.data.frame(probe_choices_LV_accuracy_category$coefficients)
LV_OR = round(exp(probe_choices_LV_accuracy_category_coef$Estimate[2]),3)

# CI
HV_CI_min=round(exp(probe_choices_HV_accuracy_category_coef$Estimate[2]-1.96*probe_choices_HV_accuracy_category_coef$`Std. Error`[2]),3)
HV_CI_max=round(exp(probe_choices_HV_accuracy_category_coef$Estimate[2]+1.96*probe_choices_HV_accuracy_category_coef$`Std. Error`[2]),3)
LV_CI_min=round(exp(probe_choices_LV_accuracy_category_coef$Estimate[2]-1.96*probe_choices_LV_accuracy_category_coef$`Std. Error`[2]),3)
LV_CI_max=round(exp(probe_choices_LV_accuracy_category_coef$Estimate[2]+1.96*probe_choices_LV_accuracy_category_coef$`Std. Error`[2]),3)
