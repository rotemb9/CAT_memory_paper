# Required R package
library(lme4)
library(ggplot2)

# Clear workspace
rm(list=ls())

# Define the local path where the data can be found
# you will need to set the correct path whre the RData file with the data is located
input_path="DEFINE PATH OF DATA HERE"

experiment_name = "bmem_short"

# load data
input_filename=paste(input_path, experiment_name, "_probe.Rdata",sep="")  
load(file=input_filename)

# calculate means
means=c(mean(with(data=subset(probe_data,PairType==1), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T),#HH
              mean(with(data=subset(probe_data,PairType==6), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T),#HM
              mean(with(data=subset(probe_data,PairType==7), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T),#LM
              mean(with(data=subset(probe_data,PairType==2), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T),#LL
              mean(with(data=subset(probe_data,PairType==3), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T), # sanity HV/LV
              mean(with(data=subset(probe_data,PairType==4 | PairType==5), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T)) # sanity same value

# Standard error function
se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }

# calculate standard error
SEM=c(se(with(data=subset(probe_data,PairType==1), tapply(Outcome, subjectID, mean, na.rm=T))), # HH
      se(with(data=subset(probe_data,PairType==6), tapply(Outcome, subjectID, mean, na.rm=T))), # HM
      se(with(data=subset(probe_data,PairType==7), tapply(Outcome, subjectID, mean, na.rm=T))), # LM
      se(with(data=subset(probe_data,PairType==2), tapply(Outcome, subjectID, mean, na.rm=T))), # LL
      se(with(data=subset(probe_data,PairType==3), tapply(Outcome, subjectID, mean, na.rm=T))), # sanity HV/LV
      se(with(data=subset(probe_data,PairType==4 | PairType==5), tapply(Outcome, subjectID, mean, na.rm=T)))) # sanity same value

# Logistic regression analysis
HH_results=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data,(probe_data$PairType3=='HH')),na.action=na.omit,family=binomial)) 
HM_results=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data,(probe_data$PairType3=='HM')),na.action=na.omit,family=binomial))
LM_results=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data,(probe_data$PairType3=='LM')),na.action=na.omit,family=binomial))
LL_results=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data,(probe_data$PairType3=='LL')),na.action=na.omit,family=binomial))
#sanity_high_low_results=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data,(probe_data$PairType==3)),na.action=na.omit,family=binomial)) 
#sanity_same_value_results=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data,(probe_data$PairType==4 | probe_data$PairType==5)),na.action=na.omit,family=binomial)) 

# Summary dataframe for plots
# ---------------------------
Results_df=as.data.frame(rbind(HH_results$coefficients[1,],HM_results$coefficients[1,],LM_results$coefficients[1,],LL_results$coefficients[1,]))
# p-value is one sided
Results_df$`Pr(>|z|)` = Results_df$`Pr(>|z|)`/2

colnames(Results_df)=colnames(HH_results$coefficients)
# add indicator for comparison type
Results_df$pairtype = c("High", "Medium-\nhigh", "Medium-\nlow", "Low")

# Add mean proportion of trials participants chose Go, and standard error of the means (SEM)
Results_df$means=NA
Results_df$means=as.vector(means[1:4])
Results_df$SEM=as.vector(SEM[1:4])
# Significance indicator
Results_df$asteriks=""
Results_df$asteriks[Results_df$`Pr(>|z|)`<0.07]="+"
Results_df$asteriks[Results_df$`Pr(>|z|)`<0.05]="*"
Results_df$asteriks[Results_df$`Pr(>|z|)`<0.01]="**"
Results_df$asteriks[Results_df$`Pr(>|z|)`<0.001]="***"
Results_df$asteriks[Results_df$Estimate<0]=""

Results_df_go_nogo=Results_df[1:4,]
value_categories = c("High", "Medium-\nhigh", "Medium-\nlow", "Low")
Results_df_go_nogo$category = value_categories

# asteriks locations
asteriks_baseline_loc=0.6
Results_df_go_nogo$asteriks_loc = asteriks_baseline_loc
Results_df_go_nogo$asteriks_loc[Results_df_go_nogo$asteriks=="+"] = asteriks_baseline_loc + 0.02

# Bar Plot - Proportions
plot_proportions=ggplot(data=Results_df_go_nogo, aes(x=category, y=means, fill=category)) +
  geom_bar(width=0.7,colour="black",position=position_dodge(0.7), stat="identity") + # Bar plot
  theme_bw() + # white background
  theme(legend.position="none") + # remove legend
  theme(axis.title.x=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank()) + # axis and background formating
  theme(aspect.ratio=2,text = element_text(size=18)) + # aspect ratio and font size
  geom_errorbar(position=position_dodge(1), width=1/4, aes(ymin=Results_df_go_nogo$means-Results_df_go_nogo$SEM, ymax=Results_df_go_nogo$means+Results_df_go_nogo$SEM))  + # add error bar of SEM
  scale_y_continuous("Proportion of trials Go items were chosen",limit=c(0,1),breaks=seq(0, 1, 0.1),expand = c(0,0)) + # define y axis properties
  scale_fill_manual(values=c("#585858","#DCDCDC", "#888888","#B8B8B8")) + # color of bars
  geom_abline(intercept = (0.5),slope=0,linetype =2, size = 1,show.legend=TRUE,aes()) + # chance level 50% reference line
  geom_text(position=position_dodge(1),aes(y=asteriks_loc,label=(asteriks)),size=12) + # significance asteriks
  scale_x_discrete(limits = value_categories)

print(plot_proportions)

# create statistics table
print("Probe statistics table")
probe_results=data.frame(row.names=c("HH","HM","LM","LL"))
probe_results$value_category = c("High", "Medium-high", "Medium-low", "Low")
probe_results$Mean=paste(round(Results_df_go_nogo$means*100,2), "%", sep="")
probe_results$SEM=paste(round(Results_df_go_nogo$SEM*100,1), "%", sep="")
odds_ratio=exp(Results_df_go_nogo$Estimate)
probe_results$one_sided_p=round(Results_df_go_nogo$`Pr(>|z|)`,3)
probe_results$odds_ratio=round(odds_ratio,3)
CI_min=exp(Results_df_go_nogo$Estimate-1.96*Results_df_go_nogo$SEM)
CI_max=exp(Results_df_go_nogo$Estimate+1.96*Results_df_go_nogo$SEM)
probe_results$CI=paste(round(CI_min,3), " - ", round(CI_max,3),sep="")
print(probe_results)
