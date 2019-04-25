function run_bmem_snacks2(sessionNum)
% function run_bmem_snacks2()
% this function runs the different sessions of the experiment bmem_snacks2
% input: sessionNum- the number of session- 1/2/3
% if no input entered, the experimenter will be asked to enter session
% number.
% session1: BDM, BDM sorting, training_demo,training, personal details, BDM resolve.
% session2: recognition_demo,recognition (with confidence levels), probe_demo,probe1, BDM2, BDM resolve,
% probe resolve.
% follow_up: recognition_demo,recognition (with confidence levels), probe_demo,probe2, BDM3, BDM resolve,
% probe resolve.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FUNCTIONS REQUIRED TO RUN PROPERLY: ----------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   --- Cue-approach codes: ---
% % %   'run_bmem_snacks_session1'
% % %   'run_bmem_snacks_session2'
% % %   'run_bmem_snacks_follow_up'
% % %   'sort_BDM'
% % %   'trainingDemo'
% % %   'training'
% % %   'personal_details
% % %   'organizeProbe'
% % %   'probeDemo'
% % %   'probe'
% % %   'disp_probeResolve'
% % %   'disp_resolveBDM'
% % %   'recognition_confidence_demo_new'
% % %   'recognition_confidence_new'
% % %   --- Other codes: ---
% % %  'CenterText'
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FOLDERS REQUIRED TO RUN PROPERLY: ------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   'Onset_files': a folder with the onset files for the probe.
% % %   'Output': the folder for the output files- results.
% % %   'Stim': with the bmp files of all the stimuli for the cue-approach
% % %    task (old stimuli).
% % %   'Instructions': a folder with the jpg images of the instructions
% % %   for each part
% % %   'Stim\recognitionNew': with the bmp files of the new stimuli
% % %   (stimuli that are not included in the cue-approach tasks, only in the
% % %   recognition task, as new stimuli).

if nargin <1 || ~ismember(sessionNum,1:3)
    sessionNum = input('Please enter session number- 1/2/3: ');
    while ~exist('sessionNum','var') || ~ismember(sessionNum,1:3)
       sessionNum = input('Please enter session number- 1/2/3: ');
    end % end while
end % end if nargin

switch sessionNum
    case 1
        run_bmem_snacks2_session1()
    case 2
        run_bmem_snacks2_session2()
    case 3
        run_bmem_snacks2_follow_up()
end % end switch sessionNum

end % end function