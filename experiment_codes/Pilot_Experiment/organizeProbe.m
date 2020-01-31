function [trialsPerRun] = organizeProbe(subjectID, order, mainPath, block, numRunsPerBlock)

% function [trialsPerRun] = organizeProbe_Israel(subjectID, order, mainPath, block, numRunsPerBlock)
%
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ================== by Rotem Botvinik Nezer March 2016 ===================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function organizes the matrices for each block of the probe session of the boost
% (cue-approach) task, divided to number of runs as requested (1 or 2 or 4 would
% work). Pay attention that number of runs should be a divisor of number of
% comparisons.

% This function is for the version where only 40 items are being trained,
% and the sanity checks are only on the NOGO items- 2*2 for HV vs. LV; 2*2
% for HV vs. HV (similar ranking); 2*2 LV vs. LV (similar ranking)


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_allstim_order*.txt'' --> created by sortBDM_Israel


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stimuliForProbe_order%d_block_%d_run%d.txt'


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID =  'bmem_snacks_999';
% subjectID =  'bmem_snacks_998'; % it is important to test both orders
% isMRI = 0;
% mainPath = pwd;
% numRunsPerBlock = 1;
% block = 1;


tic

rng shuffle

%==============================================
%% 'GLOBAL VARIABLES'
%==============================================

outputPath = [mainPath '/Output'];

%==============================================
%% 'Read in data'
%==============================================

%   'read in sorted file'
% - - - - - - - - - - - - - - - - -

file = dir([mainPath '/Output/' subjectID '_stopGoList_allstim_order*']);
fid = fopen([mainPath '/Output/' sprintf(file(length(file)).name)]);
data = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sortbdm
stimName = data{1};
% bidIndex = data{3};
% bidValue = data{4};
fclose(fid);

%==============================================
%%   'DATA ORGANIZATION'
%==============================================

% determine stimuli to use based on order number
%-----------------------------------------------------------------
switch order
    case 1
        %   comparisons of interest
        % - - - - - - - - - - - - - - -
        HH_beep =   [7 10 12 13 15 18]; % HH_beep
        HH_nobeep = [8 9 11 14 16 17]; % HH_nobeep
        
        HM_beep = [24 25 27 30 32 33]; % HM beep
        HM_nobeep = [23 26 28 29 31 34]; % HM nobeep
        
        LM_beep =   [47 50 52 53 55 58]; % LM_beep
        LM_nobeep = [48 49 51 54 56 57]; % LM_nobeep
        
        LL_beep = [64 65 67 70 72 73]; % LL_beep
        LL_nobeep = [63 66 68 69 71 74]; % LL_nobeep
        
        
        %   sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        % HV vs. LV
        sanityHV_nobeep = [5 6]; % HV_nobeep
        sanityLV_nobeep = [75 76]; % LV_nobeep
        
        % HV vs. HV
        sanityHV_nobeep1 = [19 22]; % HV_nobeep
        sanityHV_nobeep2 = [20 21]; % HV_nobeep
        
        % LV vs. LV
        sanityLV_nobeep1 = [59 62]; % LV_nobeep
        sanityLV_nobeep2 = [60 61]; % LV_nobeep
        
    case 2
        
        %   comparisons of interest
        % - - - - - - - - - - - - - - -
        HH_beep =   [8 9 11 14 16 17]; % HH_beep
        HH_nobeep = [7 10 12 13 15 18]; % HH_nobeep
        
        HM_beep = [23 26 28 29 31 34]; % HM beep 
        HM_nobeep = [24 25 27 30 32 33]; % HM nobeep     
        
        LM_beep = [48 49 51 54 56 57]; % LM_beep
        LM_nobeep =   [47 50 52 53 55 58]; % LM_nobeep
        
        LL_beep = [63 66 68 69 71 74]; % LL_beep
        LL_nobeep = [64 65 67 70 72 73]; % LL_nobeep

        %   sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        % HV vs. LV
        sanityHV_nobeep = [5 6]; % HV_nobeep
        sanityLV_nobeep = [75 76]; % LV_nobeep
        
        % HV vs. HV
        sanityHV_nobeep1 = [20 21]; % HV_nobeep
        sanityHV_nobeep2 = [19 22]; % HV_nobeep
        
        % LV vs. LV
        sanityLV_nobeep1 = [60 61]; % LV_nobeep
        sanityLV_nobeep2 = [59 62]; % LV_nobeep
        
end % end switch order


%   add multiple iterations of each item presentation
%-----------------------------------------------------
times_each_sanity_check = 2;

%   TRIAL TYPE 1 - HH: HighValue Go vs. HighValue NoGo
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numHHbeepItems = length(HH_beep);
numHHnobeepItems = length(HH_nobeep);

HH_beep_new = repmat(HH_beep,numHHbeepItems,1);
HH_beep_new = HH_beep_new(:)';
HH_nobeep_new = repmat(HH_nobeep,1,numHHnobeepItems);
HH_nobeep_new = HH_nobeep_new(:)';

[shuffle_HH_beep_new,shuff_HH_beep_new_ind] = Shuffle(HH_beep_new);
shuffle_HH_nobeep_new = HH_nobeep_new(shuff_HH_beep_new_ind);


%   TRIAL TYPE 6 - HM: HighValue Go vs. HighValue NoGo
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numHMbeepItems = length(HM_beep);
numHMnobeepItems = length(HM_nobeep);

HM_beep_new = repmat(HM_beep,numHMbeepItems,1);
HM_beep_new = HM_beep_new(:)';
HM_nobeep_new = repmat(HM_nobeep,1,numHMnobeepItems);
HM_nobeep_new = HM_nobeep_new(:)';

[shuffle_HM_beep_new,shuff_HM_beep_new_ind] = Shuffle(HM_beep_new);
shuffle_HM_nobeep_new = HM_nobeep_new(shuff_HM_beep_new_ind);


%   TRIAL TYPE 2 - LL: LowValue Go vs. LowValue NoGo
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numLLbeepItems = length(LL_beep);
numLLnobeepItems = length(LL_nobeep);

LL_beep_new = repmat(LL_beep,numLLbeepItems,1);
LL_beep_new = LL_beep_new(:)';
LL_nobeep_new = repmat(LL_nobeep,1,numLLnobeepItems);
LL_nobeep_new = LL_nobeep_new(:)';

[shuffle_LL_beep_new,shuff_LL_beep_new_ind] = Shuffle(LL_beep_new);
shuffle_LL_nobeep_new = LL_nobeep_new(shuff_LL_beep_new_ind);

%   TRIAL TYPE 7 - LM: LowValue Go vs. LowValue NoGo
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numLMbeepItems = length(LM_beep);
numLMnobeepItems = length(LM_nobeep);

LM_beep_new = repmat(LM_beep,numLMbeepItems,1);
LM_beep_new = LM_beep_new(:)';
LM_nobeep_new = repmat(LM_nobeep,1,numLMnobeepItems);
LM_nobeep_new = LM_nobeep_new(:)';

[shuffle_LM_beep_new,shuff_LM_beep_new_ind] = Shuffle(LM_beep_new);
shuffle_LM_nobeep_new = LM_nobeep_new(shuff_LM_beep_new_ind);

%   TRIAL TYPE 3: HighValue NoGo(Stop) vs. LowValue NoGo(Stop)
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numSanityHVnobeepItems = length(sanityHV_nobeep);
numSanityLVnobeepItems = length(sanityLV_nobeep);

sanityHV_nobeep_new = repmat(sanityHV_nobeep,numSanityHVnobeepItems,1);
sanityHV_nobeep_new = sanityHV_nobeep_new(:)';
sanityLV_nobeep_new = repmat(sanityLV_nobeep,1,numSanityLVnobeepItems);
sanityLV_nobeep_new = sanityLV_nobeep_new(:)';

[shuffle_sanityHV_nobeep_new,shuff_sanityHV_nobeep_new_ind] = Shuffle(sanityHV_nobeep_new);
shuffle_sanityLV_nobeep_new = sanityLV_nobeep_new(shuff_sanityHV_nobeep_new_ind);

%   TRIAL TYPE 4: HighValue NoGo vs. HighValue NoGo
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numSanityHHnobeep1Items = length(sanityHV_nobeep1);
numSanityHHnobeep2Items = length(sanityHV_nobeep2);

sanityHH_nobeep1_new = repmat(sanityHV_nobeep1,numSanityHHnobeep1Items,1);
sanityHH_nobeep1_new = sanityHH_nobeep1_new(:)';
sanityHH_nobeep2_new = repmat(sanityHV_nobeep2,numSanityHHnobeep2Items,1);
sanityHH_nobeep2_new = sanityHH_nobeep2_new(:)';

[shuffle_sanityHH_nobeep1_new,shuff_sanityHH_nobeep1_new_ind] = Shuffle(sanityHH_nobeep1_new);
shuffle_sanityHH_nobeep2_new = sanityHH_nobeep2_new(shuff_sanityHH_nobeep1_new_ind);


%   TRIAL TYPE 5: LowValue NoGo vs. LowValue NoGo
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numSanityLLnobeep1Items = length(sanityLV_nobeep1);
numSanityLLnobeep2Items = length(sanityLV_nobeep2);

sanityLL_nobeep1_new = repmat(sanityLV_nobeep1,numSanityLLnobeep1Items,1);
sanityLL_nobeep1_new = sanityLL_nobeep1_new(:)';
sanityLL_nobeep2_new = repmat(sanityLV_nobeep2,numSanityLLnobeep2Items,1);
sanityLL_nobeep2_new = sanityLL_nobeep2_new(:)';

[shuffle_sanityLL_nobeep1_new,shuff_sanityLL_nobeep1_new_ind] = Shuffle(sanityLL_nobeep1_new);
shuffle_sanityLL_nobeep2_new = sanityLL_nobeep2_new(shuff_sanityLL_nobeep1_new_ind);


%   randomize all possible comparisons for all trial types
%-----------------------------------------------------------------
numComparisonsHH = numHHbeepItems*numHHnobeepItems;
numComparisonsLL = numLLbeepItems*numLLnobeepItems;
numComparisonsHM = numHMbeepItems*numHMnobeepItems;
numComparisonsLM = numLMbeepItems*numLMnobeepItems;

numComparisons_gonogo = numComparisonsHH + numComparisonsLL + numComparisonsHM + numComparisonsLM;
numSanity = (numSanityHVnobeepItems^2 + numSanityHHnobeep1Items^2 + numSanityLLnobeep1Items^2)*times_each_sanity_check; % all 3 types of sanity checks
total_num_trials = numComparisons_gonogo + numSanity;
trialsPerRun = total_num_trials/numRunsPerBlock;

stimnum1 = zeros(numRunsPerBlock,trialsPerRun);
stimnum2 = zeros(numRunsPerBlock,trialsPerRun);
leftname = cell(numRunsPerBlock,trialsPerRun);
rightname = cell(numRunsPerBlock,trialsPerRun);
%pairType = zeros(numRunsPerBlock,trialsPerRun);


numComparisonsPerRun = numComparisons_gonogo/numRunsPerBlock;
numSanityPerRun = numSanity/numRunsPerBlock;
pairType1 = ones(numRunsPerBlock,numComparisonsHH/numRunsPerBlock);
pairType2 = 2*ones(numRunsPerBlock,numComparisonsLL/numRunsPerBlock);
pairType6 = 6*ones(numRunsPerBlock,numComparisonsHM/numRunsPerBlock);
pairType7 = 7*ones(numRunsPerBlock,numComparisonsLM/numRunsPerBlock);
pairType3 = 3*ones(numRunsPerBlock,numSanityHVnobeepItems*numSanityLVnobeepItems*times_each_sanity_check/numRunsPerBlock);
pairType4 = 4*ones(numRunsPerBlock,numSanityHHnobeep1Items*numSanityHHnobeep2Items*times_each_sanity_check/numRunsPerBlock);
pairType5 = 5*ones(numRunsPerBlock,numSanityLLnobeep1Items*numSanityLLnobeep2Items*times_each_sanity_check/numRunsPerBlock);
pairType = [pairType1 pairType2 pairType3 pairType4 pairType5 pairType6 pairType7];
% pairType(1:numRunsPerBlock,1:numComparisonsHH/numRunsPerBlock) = 1;
% pairType(1:numRunsPerBlock,1:numComparisonsPerRun/2) = 1;
% pairType(1:numRunsPerBlock,numComparisonsPerRun/2+1:numComparisonsPerRun) = 2;
% pairType(1:numRunsPerBlock,numComparisonsPerRun+1:numComparisonsPerRun+numSanityPerRun/3) = 3;
% pairType(1:numRunsPerBlock,numComparisonsPerRun+numSanityPerRun/3+1:numComparisonsPerRun+numSanityPerRun*2/3) = 4;
% pairType(1:numRunsPerBlock,numComparisonsPerRun+numSanityPerRun*2/3+1:numComparisonsPerRun+numSanityPerRun) = 5;

leftGo = ones(numRunsPerBlock,total_num_trials./numRunsPerBlock);
leftGo(:,[1:numComparisonsPerRun/4 numComparisonsPerRun/2+1:numComparisonsPerRun*3/4 1+numComparisonsPerRun:numComparisonsPerRun+numSanityPerRun/2]) = 0;

for numRun = 1:numRunsPerBlock
    pairType(numRun,:) = Shuffle(pairType(numRun,:));
    leftGo(numRun,:) = Shuffle(leftGo(numRun,:));
end % end for numRun = 1:numRunsPerBlock

HH_beep = shuffle_HH_beep_new;
HH_nobeep = shuffle_HH_nobeep_new;
HM_beep = shuffle_HM_beep_new;
HM_nobeep = shuffle_HM_nobeep_new;
LM_beep = shuffle_LM_beep_new;
LM_nobeep = shuffle_LM_nobeep_new;
LL_beep = shuffle_LL_beep_new;
LL_nobeep = shuffle_LL_nobeep_new;

sanityHV_nobeep = shuffle_sanityHV_nobeep_new;
sanityLV_nobeep = shuffle_sanityLV_nobeep_new;
sanityHV_nobeep1 = shuffle_sanityHH_nobeep1_new;
sanityHV_nobeep2 = shuffle_sanityHH_nobeep2_new;
sanityLV_nobeep1 = shuffle_sanityLL_nobeep1_new;
sanityLV_nobeep2 = shuffle_sanityLL_nobeep2_new;

sanityHV_nobeep = repmat(sanityHV_nobeep,1,times_each_sanity_check);
sanityLV_nobeep = repmat(sanityLV_nobeep,1,times_each_sanity_check);
sanityHV_nobeep1 = repmat(sanityHV_nobeep1,1,times_each_sanity_check);
sanityHV_nobeep2 = repmat(sanityHV_nobeep2,1,times_each_sanity_check);
sanityLV_nobeep1 = repmat(sanityLV_nobeep1,1,times_each_sanity_check);
sanityLV_nobeep2 = repmat(sanityLV_nobeep2,1,times_each_sanity_check);

HH = 1;
HM = 1;
LM = 1;
LL = 1;
HL_G = 1;
sanity_HH1 = 1;
sanity_LL1 = 1;

for numRun = 1:numRunsPerBlock
    
    % Create stimuliForProbe.txt for this run
    fid1 = fopen([outputPath '/' sprintf('%s_stimuliForProbe_order%d_block_%d_run%d.txt',subjectID,order,block,numRun)], 'w');
    
    for trial = 1:trialsPerRun % trial num within block
        switch pairType(numRun,trial)
            case 1
                
                % HH: HighValue Go vs. HighValue NoGo
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = HH_beep(HH);
                stimnum2(numRun,trial) = HH_nobeep(HH);
                HH = HH+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 6
                
                % HM: HighValue Go vs. HighValue NoGo
                % - - - - - - - - - - - - - - - - - - -
                stimnum1(numRun,trial) = HM_beep(HM);
                stimnum2(numRun,trial) = HM_nobeep(HM);
                HM = HM+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 2
                
                % LL: LowValue Go vs. LowValue NoGo
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = LL_beep(LL);
                stimnum2(numRun,trial) = LL_nobeep(LL);
                LL = LL+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 7
                
                % LM: LowValue Go vs. LowValue NoGo
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = LM_beep(LM);
                stimnum2(numRun,trial) = LM_nobeep(LM);
                LM = LM+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 3
                % HighValue NoGo vs. LowValue NoGo
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = sanityHV_nobeep(HL_G);
                stimnum2(numRun,trial) = sanityLV_nobeep(HL_G);
                HL_G = HL_G+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 4
                % HighValue NoGo vs. HighValue NoGo
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = sanityHV_nobeep1(sanity_HH1);
                stimnum2(numRun,trial) = sanityHV_nobeep2(sanity_HH1);
                sanity_HH1 = sanity_HH1+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 5
                % HighValue NoGo vs. HighValue NoGo
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = sanityLV_nobeep1(sanity_LL1);
                stimnum2(numRun,trial) = sanityLV_nobeep2(sanity_LL1);
                sanity_LL1 = sanity_LL1+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
        end % end switch pairtype
        
        fprintf(fid1, '%d\t %d\t %d\t %d\t %s\t %s\t \n', stimnum1(numRun,trial),stimnum2(numRun,trial),leftGo(numRun,trial),pairType(numRun,trial),leftname{numRun,trial},rightname{numRun,trial});
    end % end for trial = 1:total_num_trials
    
    fprintf(fid1, '\n');
    fclose(fid1);
end % end for numRun = 1:numRunsPerBlocks


%---------------------------------------------------------------------
% create a data structure with info about the run and all the matrices
%---------------------------------------------------------------------
outfile = strcat(outputPath,'/', sprintf('%s_stimuliForProbe_order%d_block_%d_%d_trials_%d_runs_%s.mat',subjectID,order,block,total_num_trials,numRunsPerBlock,date));

% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;
run_info.script_name = mfilename;

save(outfile);


end % end function

