# Required R package
library(lme4)
library(ggplot2)

# Clear workspace
rm(list=ls())

# Define the local path where the data can be found
# you will need to set the correct path whre the RData file with the data is located
input_file_path="DEFINE PATH OF DATA HERE"

experiment_name = "bmem_snacks2"

sessions = 2:3

# load data
if (length(sessions==2)){
  input_filename=paste(input_file_path, experiment_name, "_probe_all_sessions.Rdata",sep="")  
} else {
  input_filename=paste(input_file_path, experiment_name, "_probe_session_", as.character(sessions),".Rdata",sep="")
}
load(file=input_filename)

# loop through sessions
for (session_num in sessions){
  probe_data_curr_session = subset(probe_data, probe_data$session==session_num)

  # calculate means
  means=c(mean(with(data=subset(probe_data_curr_session,PairType==1), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T),
                mean(with(data=subset(probe_data_curr_session,PairType==2), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T),
                mean(with(data=subset(probe_data_curr_session,PairType==3), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T),
                mean(with(data=subset(probe_data_curr_session,PairType==4 | PairType==5), tapply(Outcome, subjectID, mean, na.rm=T)),na.rm=T))
  
  # Standard error function
  se = function(x) { out=sqrt(var(x, na.rm = TRUE)/length(which(!is.na(x)))) }
  
  # calculate standard error
  SEM=c(se(with(data=subset(probe_data_curr_session,PairType==1), tapply(Outcome, subjectID, mean, na.rm=T))),
        se(with(data=subset(probe_data_curr_session,PairType==2), tapply(Outcome, subjectID, mean, na.rm=T))),
        se(with(data=subset(probe_data_curr_session,PairType==3), tapply(Outcome, subjectID, mean, na.rm=T))),
        se(with(data=subset(probe_data_curr_session,PairType==4 | PairType==5), tapply(Outcome, subjectID, mean, na.rm=T))))
  
  # Logistic regression analysis
  high_value_results=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data_curr_session,(probe_data_curr_session$PairType2=='High_Value')),na.action=na.omit,family=binomial)) 
  low_value_results=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data_curr_session,(probe_data_curr_session$PairType2=='Low_Value')),na.action=na.omit,family=binomial))
  high_low_difference_results=summary(glmer(Outcome ~ 1 + PairType + (1 + PairType|subjectID),data=subset(probe_data_curr_session,probe_data_curr_session$PairType %in% c(1,2)),na.action=na.omit,family=binomial)) #effect of Go choice for HV vs LV
  sanity_high_low_results=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data_curr_session,(probe_data_curr_session$PairType==3)),na.action=na.omit,family=binomial)) 
  sanity_same_value_results=summary(glmer(Outcome ~ 1 + (1|subjectID),data=subset(probe_data_curr_session,(probe_data_curr_session$PairType==4 | probe_data_curr_session$PairType==5)),na.action=na.omit,family=binomial)) 
  
  # Summary dataframe for plots
  # ---------------------------
  Results_df=as.data.frame(rbind(high_value_results$coefficients[1,],low_value_results$coefficients[1,],high_low_difference_results$coefficients[2,],sanity_high_low_results$coefficients[1,],sanity_same_value_results$coefficients[1,]))
  # p-value is one sided
  Results_df$`Pr(>|z|)` = Results_df$`Pr(>|z|)`/2
  
  colnames(Results_df)=colnames(high_value_results$coefficients)
  # add indicator for comparison type
  Results_df$pairtype = c("high-value", "low-value", "high-low difference", "sanity high-low", "sanity same value")
  
  # Add mean proportion of trials participants chose Go, and standard error of the means (SEM)
  Results_df$means=NA
  Results_df$means[!Results_df$pairtype=="high-low difference"]=as.vector(means)
  Results_df$SEM[!Results_df$pairtype=="high-low difference"]=as.vector(SEM)
  # Significance indicator
  Results_df$asteriks=""
  Results_df$asteriks[Results_df$`Pr(>|z|)`<0.07]="+"
  Results_df$asteriks[Results_df$`Pr(>|z|)`<0.05]="*"
  Results_df$asteriks[Results_df$`Pr(>|z|)`<0.01]="**"
  Results_df$asteriks[Results_df$`Pr(>|z|)`<0.001]="***"
  
  Results_df$session=as.factor(session_num) 
  
  if (session_num == sessions[1]){
    Results_df_all_sessions=Results_df
  } else {
    Results_df_all_sessions=rbind(Results_df_all_sessions, Results_df)
  }

} # up to here the loop on sessions


Results_df_go_nogo=Results_df_all_sessions[Results_df_all_sessions$pairtype=='high-value' | Results_df_all_sessions$pairtype=='low-value',]
Results_df_go_nogo$label=paste("Session",as.character(Results_df_go_nogo$session),sep=" ")
Results_df_go_nogo$label[Results_df_go_nogo$label=="Session 3"] = "Follow-up Session"
Results_df_go_nogo$label=factor(Results_df_go_nogo$label, levels = c("Session 1","Session 2", "Follow-up Session"))

# asteriks locations
asteriks_baseline_loc=0.7
Results_df_go_nogo$asteriks_loc = asteriks_baseline_loc
Results_df_go_nogo$asteriks_loc[Results_df_go_nogo$asteriks=="+"] = asteriks_baseline_loc + 0.02

# Bar Plot - Proportions
plot_proportions=ggplot(data=Results_df_go_nogo, aes(x=label, y=means, fill=pairtype)) +
  geom_bar(width=0.7,colour="black",position=position_dodge(0.7), stat="identity") + # Bar plot
  theme_bw() + # white background
  theme(legend.position="none") + # remove legend
  theme(axis.title.x=element_blank(),axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank()) + # axis and background formating
  theme(aspect.ratio=2/length(sessions),text = element_text(size=24)) + # aspect ratio and font size
  geom_errorbar(position=position_dodge(0.7), width=1/4, aes(ymin=Results_df_go_nogo$means-Results_df_go_nogo$SEM, ymax=Results_df_go_nogo$means+Results_df_go_nogo$SEM))  + # add error bar of SEM
  scale_y_continuous("Proportion of trials Go items were chosen",limit=c(0,1),breaks=seq(0, 1, 0.1),expand = c(0,0)) + # define y axis properties
  scale_fill_manual(values=c("#585858","#D8D8D8")) + # color of bars
  geom_abline(intercept = (0.5),slope=0,linetype =2, size = 1,show.legend=TRUE,aes()) + # chace level 50% reference line
  geom_text(position=position_dodge(0.7),aes(y=asteriks_loc,label=(asteriks)),size=12) # significance asteriks
# add high-value - low-value differential effect significance asteriks
for (i in sessions) {
  row_loc=(i-2)*5+3
  if (Results_df_all_sessions$`Pr(>|z|)`[row_loc]<0.05 & Results_df_all_sessions$Estimate[row_loc]<0) {
    session_ind=c(i,i+length(sessions))
    Lines_hight=asteriks_baseline_loc + 0.1
    i=i-1
    tmp_df=data.frame(x_val=c(i-1/4,i-1/4,i+1/4,i+1/4),y_val=c(Lines_hight-0.01,Lines_hight,Lines_hight,Lines_hight-0.01)) # define shape of an open rectangle above the bar
    tmp_df$pairtype=NA
    plot_proportions = plot_proportions +
      geom_line(data = tmp_df, aes(x=x_val,y = y_val)) + # draw open rectangle
      annotate("text", x = i, y = Lines_hight+0.02, label = (Results_df_all_sessions$asteriks[row_loc]),size=12) # differential effect significance asteriks
  }
}
dev.new(width=1.5*length(sessions), height=8)
print(plot_proportions)

# create statistics table
print("Probe statistics table")
probe_results=data.frame(row.names=1:(2*length(sessions)))
probe_results$Session=c('Session 2', 'Session 2', 'Follow-up Session', 'Follow-up Session')
probe_results$value_category=c("High", "Low")
probe_results$Mean=paste(round(Results_df_go_nogo$means*100,2), "%", sep="")
probe_results$SEM=paste(round(Results_df_go_nogo$SEM*100,1),"%", sep="")
probe_results$P=round(Results_df_go_nogo$`Pr(>|z|)`,3)
probe_results$odds_ratio=round(exp(Results_df_go_nogo$Estimate),3)
CI_min=exp(Results_df_go_nogo$Estimate-1.96*Results_df_go_nogo$SEM)
CI_max=exp(Results_df_go_nogo$Estimate+1.96*Results_df_go_nogo$SEM)
probe_results$CI=paste(round(CI_min,3), " - ", round(CI_max, 3), sep="")
print(probe_results)
