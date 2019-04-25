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
        HV_beep =   [7 10 12 13 15 18]; % HV_beep
        HV_nobeep = [8 9 11 14 16 17]; % HV_nobeep
        
        LV_beep =   [44 45 47 50 52 53]; % LV_beep
        LV_nobeep = [43 46 48 49 51 54]; % LV_nobeep
        
        
        %   sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        % HV vs. LV
        sanityHV_nobeep = [5 6]; % HV_nobeep
        sanityLV_nobeep = [55 56]; % LV_nobeep
        
        % HV vs. HV
        sanityHH_nobeep1 = [19 22]; % HV_nobeep
        sanityHH_nobeep2 = [20 21]; % HV_nobeep
        
        % LV vs. LV
        sanityLL_nobeep1 = [39 42]; % LV_nobeep
        sanityLL_nobeep2 = [40 41]; % LV_nobeep
        
    case 2
        
        %   comparisons of interest
        % - - - - - - - - - - - - - - -
        HV_beep =   [8 9 11 14 16 17]; % HV_beep
        HV_nobeep = [7 10 12 13 15 18]; % HV_nobeep
        
        
        LV_beep =   [43 46 48 49 51 54]; % LV_beep
        LV_nobeep = [44 45 47 50 52 53]; % LV_nobeep
        
        
        %   sanity check comparisons - just NOGO
        % - - - - - - - - - - - - - -
        % HV vs. LV
        sanityHV_nobeep = [5 6]; % HV_nobeep
        sanityLV_nobeep = [55 56]; % LV_nobeep
        
        % HV vs. HV
        sanityHH_nobeep1 = [20 21]; % HV_nobeep
        sanityHH_nobeep2 = [19 22]; % HV_nobeep
        
        % LV vs. LV
        sanityLL_nobeep1 = [40 41]; % LV_nobeep
        sanityLL_nobeep2 = [39 42]; % LV_nobeep
        
end % end switch order


%   add multiple iterations of each item presentation
%-----------------------------------------------------


%   TRIAL TYPE 1: HighValue Go vs. HighValue NoGo(Stop)
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numHVbeepItems = length(HV_beep);
numHVnobeepItems = length(HV_nobeep);

HV_beep_new = repmat(HV_beep,numHVbeepItems,1);
HV_beep_new = HV_beep_new(:)';
HV_nobeep_new = repmat(HV_nobeep,1,numHVnobeepItems);
HV_nobeep_new = HV_nobeep_new(:)';

[shuffle_HV_beep_new,shuff_HV_beep_new_ind] = Shuffle(HV_beep_new);
shuffle_HV_nobeep_new = HV_nobeep_new(shuff_HV_beep_new_ind);



%   TRIAL TYPE 2: LowValue Go vs. LowValue NoGo(Stop)
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numLVbeepItems = length(LV_beep);
numLVnobeepItems = length(LV_nobeep);

LV_beep_new = repmat(LV_beep,numLVbeepItems,1);
LV_beep_new = LV_beep_new(:)';
LV_nobeep_new = repmat(LV_nobeep,1,numLVnobeepItems);
LV_nobeep_new = LV_nobeep_new(:)';

[shuffle_LV_beep_new,shuff_LV_beep_new_ind] = Shuffle(LV_beep_new);
shuffle_LV_nobeep_new = LV_nobeep_new(shuff_LV_beep_new_ind);


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
numSanityHHnobeep1Items = length(sanityHH_nobeep1);
numSanityHHnobeep2Items = length(sanityHH_nobeep2);

sanityHH_nobeep1_new = repmat(sanityHH_nobeep1,numSanityHHnobeep1Items,1);
sanityHH_nobeep1_new = sanityHH_nobeep1_new(:)';
sanityHH_nobeep2_new = repmat(sanityHH_nobeep2,numSanityHHnobeep2Items,1);
sanityHH_nobeep2_new = sanityHH_nobeep2_new(:)';

[shuffle_sanityHH_nobeep1_new,shuff_sanityHH_nobeep1_new_ind] = Shuffle(sanityHH_nobeep1_new);
shuffle_sanityHH_nobeep2_new = sanityHH_nobeep2_new(shuff_sanityHH_nobeep1_new_ind);


%   TRIAL TYPE 5: LowValue NoGo vs. LowValue NoGo
% - - - - - - - - - - - - - - - - - - - - - - - - - - -
numSanityLLnobeep1Items = length(sanityLL_nobeep1);
numSanityLLnobeep2Items = length(sanityLL_nobeep2);

sanityLL_nobeep1_new = repmat(sanityLL_nobeep1,numSanityLLnobeep1Items,1);
sanityLL_nobeep1_new = sanityLL_nobeep1_new(:)';
sanityLL_nobeep2_new = repmat(sanityLL_nobeep2,numSanityLLnobeep2Items,1);
sanityLL_nobeep2_new = sanityLL_nobeep2_new(:)';

[shuffle_sanityLL_nobeep1_new,shuff_sanityLL_nobeep1_new_ind] = Shuffle(sanityLL_nobeep1_new);
shuffle_sanityLL_nobeep2_new = sanityLL_nobeep2_new(shuff_sanityLL_nobeep1_new_ind);


%   randomize all possible comparisons for all trial types
%-----------------------------------------------------------------
numComparisonsHV = numHVbeepItems^2;
numComparisonsLV = numLVbeepItems^2;
numComparisons = numComparisonsHV + numComparisonsLV;
numSanity = numSanityHVnobeepItems^2 + numSanityHHnobeep1Items^2 + numSanityLLnobeep1Items^2; % all 3 types of sanity checks
total_num_trials = numComparisons + numSanity;
trialsPerRun = total_num_trials/numRunsPerBlock;

stimnum1 = zeros(numRunsPerBlock,trialsPerRun);
stimnum2 = zeros(numRunsPerBlock,trialsPerRun);
leftname = cell(numRunsPerBlock,trialsPerRun);
rightname = cell(numRunsPerBlock,trialsPerRun);
pairType = zeros(numRunsPerBlock,trialsPerRun);


numComparisonsPerRun = numComparisons/numRunsPerBlock;
numSanityPerRun = numSanity/numRunsPerBlock;
pairType(1:numRunsPerBlock,1:numComparisonsPerRun/2) = 1;
pairType(1:numRunsPerBlock,numComparisonsPerRun/2+1:numComparisonsPerRun) = 2;
pairType(1:numRunsPerBlock,numComparisonsPerRun+1:numComparisonsPerRun+numSanityPerRun/3) = 3;
pairType(1:numRunsPerBlock,numComparisonsPerRun+numSanityPerRun/3+1:numComparisonsPerRun+numSanityPerRun*2/3) = 4;
pairType(1:numRunsPerBlock,numComparisonsPerRun+numSanityPerRun*2/3+1:numComparisonsPerRun+numSanityPerRun) = 5;

leftGo = ones(numRunsPerBlock,total_num_trials./numRunsPerBlock);
leftGo(:,[1:numComparisonsPerRun/4 numComparisonsPerRun/2+1:numComparisonsPerRun*3/4 1+numComparisonsPerRun:numComparisonsPerRun+numSanityPerRun/2]) = 0;

for numRun = 1:numRunsPerBlock
    pairType(numRun,:) = Shuffle(pairType(numRun,:));
    leftGo(numRun,:) = Shuffle(leftGo(numRun,:));
end % end for numRun = 1:numRunsPerBlock

HV_beep = shuffle_HV_beep_new;
HV_nobeep = shuffle_HV_nobeep_new;
LV_beep = shuffle_LV_beep_new;
LV_nobeep = shuffle_LV_nobeep_new;

sanityHV_nobeep = shuffle_sanityHV_nobeep_new;
sanityLV_nobeep = shuffle_sanityLV_nobeep_new;
sanityHH_nobeep1 = shuffle_sanityHH_nobeep1_new;
sanityHH_nobeep2 = shuffle_sanityHH_nobeep2_new;
sanityLL_nobeep1 = shuffle_sanityLL_nobeep1_new;
sanityLL_nobeep2 = shuffle_sanityLL_nobeep2_new;

HH = 1;
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
                
                % HighValue Go vs. HighValue NoGo(Stop)
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = HV_beep(HH);
                stimnum2(numRun,trial) = HV_nobeep(HH);
                HH = HH+1;
                if leftGo(numRun,trial) == 1
                    leftname(numRun,trial) = stimName(stimnum1(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum2(numRun,trial));
                else
                    leftname(numRun,trial) = stimName(stimnum2(numRun,trial));
                    rightname(numRun,trial) = stimName(stimnum1(numRun,trial));
                end
                
            case 2
                
                % LowValue Go vs. LowValue NoGo(Stop)
                % - - - - - - - - - - - - - - - - - - -
                
                stimnum1(numRun,trial) = LV_beep(LL);
                stimnum2(numRun,trial) = LV_nobeep(LL);
                LL = LL+1;
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
                
                stimnum1(numRun,trial) = sanityHH_nobeep1(sanity_HH1);
                stimnum2(numRun,trial) = sanityHH_nobeep2(sanity_HH1);
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
                
                stimnum1(numRun,trial) = sanityLL_nobeep1(sanity_LL1);
                stimnum2(numRun,trial) = sanityLL_nobeep2(sanity_LL1);
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

