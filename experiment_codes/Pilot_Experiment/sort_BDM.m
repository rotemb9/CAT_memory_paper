function [] = sort_BDM(subjectID,order,outputPath)

% function [] = sort_BDM_Israel(subjectID,order,outputPath)

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Rotem Botvinik May 2015 ====================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function sorts the stimuli according to the BDM results.
% This function is a version in which only the 40 of the items are included
% in the training


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   [mainPath '\Output\' subjectID '_BDM1_*.txt']


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_allstim_order%d.txt', order


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% order = 1;
% outputPath = [pwd '/Output'];
% subjectNum = '999';

tic

rng shuffle

num_stimuli = 80; % how many items were ranked

%=========================================================================
%%  read in info from BDM1.txt
%=========================================================================

BDM1_file = dir([outputPath '/' subjectID '_BDM1*.txt']);
fid = fopen([outputPath '/' BDM1_file(end).name]);
BDM1_data = textscan(fid, '%s\t%d\t%f\t%s\t%f\t%f\t%f' , 'HeaderLines', 1); %read in data as new matrix   
fclose(fid);


%=========================================================================
%%  Create matrix sorted by descending bid value
%========================================================================

[bids_sort,trialnum_sort_bybid] = sort(BDM1_data{5},'descend');

bid_sortedM(:,1) = trialnum_sort_bybid; % trialnums organized by descending bid amt
bid_sortedM(:,2) = bids_sort; % bids sorted large to small
bid_sortedM(:,3) = 1:num_stimuli; % stimrank

stimnames_sorted_by_bid = BDM1_data{4}(trialnum_sort_bybid);


%=========================================================================
%%   The ranking of the stimuli determine the stimtype
%=========================================================================

if order == 1

    bid_sortedM([       7 10 12 13 15 18      ], 4) = 11; % HH_beep
    bid_sortedM([       8  9 11 14 16 17      ], 4) = 12; % HH_nobeep
    bid_sortedM([      24 25 27 30 32 33      ], 4) = 11; % HM_beep
    bid_sortedM([      23 26 28 29 31 34      ], 4) = 12; % HM_nobeep
    bid_sortedM([      47 50 52 53 55 58      ], 4) = 22; % LM_beep 
    bid_sortedM([      48 49 51 54 56 57      ], 4) = 24; % LM_nobeep    
    bid_sortedM([      64 65 67 70 72 73      ], 4) = 22; % LL_beep 
    bid_sortedM([      63 66 68 69 71 74      ], 4) = 24; % LL_nobeep
    bid_sortedM([         1:4    35:40        ], 4) = 12; % HV_nobeep not in probe
    bid_sortedM([        41:46   77:80        ], 4) = 24; % LV_nobeep not in probe
    bid_sortedM([         5:6    19:22        ], 4) = 12; % HV_nobeep sanity
    bid_sortedM([        59:62   75:76        ], 4) = 24; % LV_nobeep sanity   

    else

    bid_sortedM([       8  9 11 14 16 17      ], 4) = 11; % HH_beep
    bid_sortedM([       7 10 12 13 15 18      ], 4) = 12; % HH_nobeep
    bid_sortedM([      23 26 28 29 31 34      ], 4) = 11; % HM_beep
    bid_sortedM([      24 25 27 30 32 33      ], 4) = 12; % HM_nobeep
    bid_sortedM([      48 49 51 54 56 57      ], 4) = 22; % LM_beep 
    bid_sortedM([      47 50 52 53 55 58      ], 4) = 24; % LM_nobeep    
    bid_sortedM([      63 66 68 69 71 74      ], 4) = 22; % LL_beep 
    bid_sortedM([      64 65 67 70 72 73      ], 4) = 24; % LL_nobeep
    bid_sortedM([         1:4    35:40        ], 4) = 12; % HV_nobeep not in probe
    bid_sortedM([        41:46   77:80        ], 4) = 24; % LV_nobeep not in probe
    bid_sortedM([         5:6    19:22        ], 4) = 12; % HV_nobeep sanity
    bid_sortedM([        59:62   75:76        ], 4) = 24; % LV_nobeep sanity 

end % end if order == 1

%=========================================================================
%%  create stopGoList_allstim.txt
%   this file is used during probe
%=========================================================================

fid2 = fopen([outputPath '/' subjectID sprintf('_stopGoList_allstim_order%d.txt', order)], 'w');    

for i = 1:length(bid_sortedM)
    fprintf(fid2, '%s\t%d\t%d\t%d\t%d\t\n', stimnames_sorted_by_bid{i,1},bid_sortedM(i,4),bid_sortedM(i,3),bid_sortedM(i,2),bid_sortedM(i,1)); 
end
fprintf(fid2, '\n');
fclose(fid2);  

end % end function