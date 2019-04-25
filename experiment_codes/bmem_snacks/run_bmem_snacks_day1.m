function run_bmem_snacks_day1()

% function run_bmem_snacks_day1()

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ================= by Rotem Botvinik Nezer July 2016 ===============
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function runs all the parts of the bmem_snacks experiment
% DAY1: BDM, BDM sorting, training, personal details, BDM fractals, probe, BDM resolve,
% probe resolvse.
% DAY2: recognition (with confidence levels), probe, BDM2, BDM resolve,
% probe resolve.
% DAY30: recognition (with confidence levels), probe, BDM3, BDM resolve,
% probe resolve.

% The try-catch is because the mac caused an error with opening the screen
% from time to time, so we want to prevent it from failing.

% This version is for running only 40 items in training!
% 20 runs in the training session

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FUNCTIONS REQUIRED TO RUN PROPERLY: ----------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   --- Cue-approach codes: ---
% % %   'sort_BDM'
% % %   'trainingDemo'
% % %   'training'
% % %   'personal_details
% % %   'organizeProbe'
% % %   'probeDemo'
% % %   'probe'
% % %   'disp_probeResolve'
% % %   'disp_resolveBDM'
% % %   'recognition_confidence'

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

tic

rng shuffle

% =========================================================================
% Get input args and check if input is ok
% =========================================================================

% %---dummy info for debugging purposes --------
% subjectNum =  '999'; % for debugging. Real subjects should start from 101
% and on
% subjectNum = '998'; % to test both order 1 and 2
% isMRI = 0; % 0- not MRI; 1- MRI experiment (and then test_comp = 1)
% sessionNum = 1; % if not follow-up
% mainPath = pwd;

% - - - - - - - - - - - - - - - - - - - - - - - -
% fixed variables for this behavioral experiment
% - - - - - - - - - - - - - - - - - - - - - - - -
experimentName = 'bmem_snacks'; % hard-coded since true for all subjects
sessionNum = 1; % this is the code for day1 - the first session
isMRI = 0; % this is a behavioral study. change if switch to MRI...
mainPath = pwd;
outputPath = [mainPath '/Output'];

% get time and date
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];

subjectID_ok = 0;

while subjectID_ok == 0
    subjectNum = input('Subject number (only the digits): ');
    while isempty(subjectNum)
        disp('ERROR: no value entered. Please try again.');
        subjectNum = input('Subject number (only the digits):');
    end
    % order number (same across all tasks\runs for a single subject. Should be
    % counterbalanced between 1,2 across subjects)
    % define order according to subjectNum being even/odd
    if mod(subjectNum, 2) == 0
        order = 2;
    else
        order = 1;
    end
    
    subjectID = [experimentName '_' num2str(subjectNum)];
    disp(['subjectID is: ' num2str(subjectNum)]);
    disp(['order is: ' num2str(order)]);
    disp(['session is: ' num2str(sessionNum)]);
    
    % read files from the output folder to make sure the subject number is
    % correct and there aren't any files for this subject (if session1) or
    % there are files from previous sessions and not for the current one (for
    % follow up sessions)
    subject_files = dir([outputPath '/' subjectID '*']);
    if ~isempty(subject_files)
        warning_msg = ['There are already ' num2str(length(subject_files)) ' files for this subject- ' subjectID '. Please make sure you entered the right number'];
        set(groot,'defaultUicontrolFontSize', 16);
        warning_answer = questdlg(warning_msg,'Warning!','it is OK', 'change subject number','it is OK');
    else
        subjectID_ok = 1;
    end
    
    if exist('warning_msg', 'var') && strcmp(warning_answer, 'it is OK')
        subjectID_ok = 1;
    end
    
end

% open a txt file for crashing logs
fid_crash = fopen([outputPath '/' subjectID '_crashingLogs' num2str(sessionNum) '_' timestamp '.txt'], 'a');

% =========================================================================
% PART 1: BDM (python)
% =========================================================================
system(['/usr/local/bin/python2.7 ' mainPath '/BDM_demo.py ' subjectID]);

system(['/usr/local/bin/python2.7 ' mainPath '/BDM.py ' subjectID]);

% =========================================================================
% Sort stimuli according to the BDM ranking (BDM = PART 1)
% =========================================================================
sort_BDM(subjectID,order,outputPath);


% =========================================================================
% Training (including demo) - (Training = PART 2)
% =========================================================================

% Set number of runs
% -------------------------------------------------
total_num_runs_training = 20;


% for debugging:
% total_num_runs_training = 4;

runInd = 1; % The first run to run with the training function (should be changes if for example the first runs are with eye-tracking)

crashedDemoTraining = 0;
keepTrying = 1;
while keepTrying < 10
    try
        trainingDemo(subjectID,order,mainPath,isMRI);
        keepTrying = 10;
    catch
        sca;
        crashedDemoTraining = crashedDemoTraining + 1;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - TRAINING DEMO!');
    end
end
fprintf(fid_crash,'training demo crashed:\t %d\n', crashedDemoTraining);

set(groot,'defaultUicontrolFontSize', 16)
demo_again = questdlg('Do you want more practice?','No','Yes','No','Yes');
if strcmp(demo_again, 'No')
    demo_again = 0;
else
    demo_again = 1;
end

if demo_again == 1
    trainingDemo(subjectID,order,mainPath,isMRI);
end

Ladder1IN = 750;
Ladder2IN = 750;

crashedTraining = 0;
keepTrying = 1;
while keepTrying < 10
    try
        training(subjectID,order,mainPath,isMRI,runInd,total_num_runs_training,Ladder1IN,Ladder2IN);
        keepTrying = 10;
    catch
        sca;
        crashedTraining = crashedTraining + 1;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - TRAINING!');
    end
end
fprintf(fid_crash,'training crashed:\t %d\n', crashedTraining);

%==========================================================
%%   'Part 3.1 - Personal Details'
%==========================================================
% PART 3
% Getting personal details from the subject (gender, dominant hand, height,
% weight, occupation
personal_details(subjectID, order, outputPath, sessionNum);

%==========================================================
%%   'Part 3.2 - Fractals - liking'
%==========================================================
% This part is used to give a break between the training and probe, to
% allow consolidation (at least partly).
% This part is used to get data regarding subjects' liking of other stimuli
% for the other experiments, such as faces / fractals.

system(['/usr/local/bin/python2.7 ' mainPath '/BDM_fractals.py ' subjectID]);

questdlg('Thank you!!! Please call the experimenter','This part is over','Continue','Continue');


%=================================
%%   'Part 4 - probe_demo & probe'
%=================================

numBlocks = 2; % Define how many blocks for the probe session. Each block includes all comparisons, one time each.
numRunsPerBlock = 1; % Define the required number of runs per block
% This is the version with only 84 comparisons (36 HV GO-NOGO; 36 LV
% GO-NOGO; 4 sanity NOGO HV-LV, 4 sanity NOGO HV, 4 sanity NOGO LV), so the block should not be divided into runs

crashedDemoProbe = 0;
keepTrying = 1;
while keepTrying < 10
    try
        probeDemo(subjectID, order, mainPath, isMRI, sessionNum);
        keepTrying = 10;
    catch
        sca;
        crashedDemoProbe = crashedDemoProbe + 1;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - PROBE DEMO!');
    end
end
fprintf(fid_crash,'probe demo crashed:\t %d\n', crashedDemoProbe);

set(groot,'defaultUicontrolFontSize', 16)
demo_again = questdlg('Do you want more practice?','No','Yes','No','Yes');
if strcmp(demo_again, 'No')
    demo_again = 0;
else
    demo_again = 1;
end

if demo_again==1
    probeDemo(subjectID, order, mainPath, isMRI, sessionNum);
end

% Run blocks. Before each block, stimuli need to be organized in a list and
% divided to the required number of runs
for ind = 1:numBlocks
    block = (sessionNum-1)*numBlocks + ind;
    % Organize the stimuli for the probe
    % ===================================
    [trialsPerRun] = organizeProbe(subjectID, order, mainPath, block, numRunsPerBlock);
    crashedProbe = 0;
    for numRun = 1:numRunsPerBlock
        keepTrying = 1;
        while keepTrying < 10
            try
                probe(subjectID, order, mainPath, isMRI, sessionNum, block, numRun, numRunsPerBlock, trialsPerRun);
                keepTrying=10;
            catch
                sca;
                keepTrying = keepTrying+1;
                crashedProbe = crashedProbe + 1;
                disp('CODE HAD CRASHED - PROBE!');
            end
        end
        fprintf(fid_crash,'probe crashed:\t %d\n', crashedProbe);
    end % end for numRun = 1:numRunsPerBlock
end % end for block = 1:numBlocks


%==========================================================
%%   'post-task'
%==========================================================

% The experimenter needs to press '123' for the results to be shown

msgbox('Thank you!!! Please call the experimenter to finish the experiment.');

code_to_continue = 0;
while code_to_continue ~= 123
    code_to_continue = inputdlg('Enter the code when you are ready for the probe lottery');
    code_to_continue = cell2mat(code_to_continue);
    code_to_continue = str2double(code_to_continue);
end

disp_resolveBDM(subjectID, mainPath, sessionNum);

code_to_continue = 0;
while code_to_continue ~= 123
    code_to_continue = inputdlg('Enter the code when you are ready for the probe lottery');
    code_to_continue = cell2mat(code_to_continue);
    code_to_continue = str2double(code_to_continue);
end

disp_probeResolve(subjectID, sessionNum, outputPath, numBlocks);


fclose(fid_crash);

end % end function
