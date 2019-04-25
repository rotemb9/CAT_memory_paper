function run_bmem_short(debugging_mode)

% function run_bmem_short(debugging_mode)
%
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ================= by Rotem Botvinik Nezer July 2017 ===============
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%
% This function runs all the parts of the bmem_short experiment
% BDM, BDM sorting, training (one run), personal details,
% filler tasks (ranking fractals and familiar and unfamiliar faces)
% recognition,probe, familiarity, BDM resolve, probe resolvse.
%
% The try-catch is because the mac caused an error with opening the screen
% from time to time, so we want to prevent it from failing.
%
% This version is for running only 80 items in BDM and training
% only 1 in the training session
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FUNCTIONS REQUIRED TO RUN PROPERLY: ----------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   --- Cue-approach codes: ---
% % %   'BDM'
% % %   'sort_BDM'
% % %   'trainingDemo'
% % %   'training'
% % %   'personal_details'
% % %   'ScaleRanking'
% % %   'ScaleRanking_familiarity'
% % %   'recognition_confidence'
% % %   'recognition_confidence_demo'
% % %   'organizeProbe'
% % %   'probeDemo'
% % %   'probe'
% % %   'familiarity'
% % %   'disp_probeResolve'
% % %   'disp_resolveBDM'
%
% % %   --- Other codes: ---
% % %  'CenterText'
% % %  'Eventflag'
% % %  'GenFlag'
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ---------------- FOLDERS REQUIRED TO RUN PROPERLY: ------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % %   'Onset_files': a folder with the onset files for the training and the probe.
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

if nargin < 1
    debugging_mode = 0;
end

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
experimentName = 'bmem_short'; % hard-coded since true for all subjects
sessionNum = 1; % this is the code for day1 - the first session
isMRI = 0; % this is a behavioral study. change if switch to MRI...
mainPath = pwd;
outputPath = [mainPath '/Output'];
use_eyetracker_general = 1; % in this version we use eye-tracking

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

% =========================================================================
% PART 1: BDM
% =========================================================================

% Ask if you want to use Eye Tracker (DIALOG BOX)
% if use_eyetracker_general==1
%     use_eyetracker = 0;
%     ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'BDM demo', 'Yes', 'No', 'No');
%     if strcmp(ask_if_want_eyetracker, 'Yes')
%         use_eyetracker = 1;
%     end
% end

use_eyetracker = 0;

% demo
keepTrying = 1;
while keepTrying < 10
    try
        BDM(subjectID, sessionNum, 'demo', use_eyetracker)
        keepTrying = 10;
    catch
        sca;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - BDM DEMO!');
    end
end

% Ask if you want to use Eye Tracker (DIALOG BOX)
if use_eyetracker_general==1
    use_eyetracker = 0;
    ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'Yes', 'Yes', 'No', 'No');
    if strcmp(ask_if_want_eyetracker, 'Yes')
        use_eyetracker = 1;
    end
end

% full
keepTrying = 1;
while keepTrying < 10
    try
        BDM(subjectID, sessionNum, 'BDM', use_eyetracker, debugging_mode)
        keepTrying = 10;
    catch
        sca;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - BDM!');
    end
end

% =========================================================================
% Sort stimuli according to the BDM ranking
% =========================================================================
sort_BDM(subjectID,order,outputPath);


% =========================================================================
% Training (including demo) - (Training = PART 2)
% =========================================================================

% Set number of runs
% -------------------------------------------------
total_num_runs_training = 1;
num_parts_of_training = 1; % we use eyetracking and want to calibrate and validat it once at the middle of training
parts = linspace(1,total_num_runs_training,num_parts_of_training+1);
runInd = ceil(parts(1:end-1));
last_run = [ceil(parts(2:end-1)-1) parts(end)];

%% Ask if you want to use Eye Tracker (DIALOG BOX)
% ================================================
% if use_eyetracker_general==1
%     use_eyetracker = 0;
%     ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'training demo', 'Yes', 'No', 'No');
%     if strcmp(ask_if_want_eyetracker, 'Yes')
%         use_eyetracker = 1;
%     end
% end

use_eyetracker = 0;

% demo
keepTrying = 1;
while keepTrying < 10
    try
        trainingDemo(subjectID,order,mainPath,isMRI, use_eyetracker);
        keepTrying = 10;
    catch
        sca;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - TRAINING DEMO!');
    end
end   

set(groot,'defaultUicontrolFontSize', 16)
demo_again = questdlg('Do you want more practice?','training demo','Yes','No','Yes');
if strcmp(demo_again, 'No')
    demo_again = 0;
else
    demo_again = 1;
end

if demo_again == 1
    keepTrying = 1;
    while keepTrying < 10
        try
            trainingDemo(subjectID,order,mainPath,isMRI, use_eyetracker);
            keepTrying = 10;
        catch
            sca;
            keepTrying = keepTrying + 1;
            disp('CODE HAD CRASHED - TRAINING DEMO!');
        end
    end
end

Ladder1IN = 750;
Ladder2IN = 750;

if use_eyetracker_general==1
    use_eyetracker = 0;
    ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'training', 'Yes', 'No', 'No');
    if strcmp(ask_if_want_eyetracker, 'Yes')
        use_eyetracker = 1;
    end
end

for part_num = 1:num_parts_of_training
    keepTrying = 1;
    while keepTrying < 10
        try
            [Ladder1IN,Ladder2IN] = training(subjectID,order,mainPath,isMRI,runInd(part_num),total_num_runs_training,Ladder1IN,Ladder2IN,use_eyetracker,last_run(part_num), debugging_mode);
            keepTrying = 10;
        catch
            sca;
            keepTrying = keepTrying + 1;
            disp('CODE HAD CRASHED - TRAINING!');
        end
    end
    if part_num ~= num_parts_of_training
        questdlg('Please call the experimenter.','','OK','OK');
    end
end

%==========================================================
%%   'Part 3 - Personal Details'
%==========================================================
% PART 3
% Getting personal details from the subject (gender, dominant hand, height,
% weight, occupation
personal_details(subjectID, order, outputPath, sessionNum);

%==========================================================
%%   'Part 4 & 5 - filler tasks - scale rankings and familiarity'
%==========================================================
% This part is used to give a break between the training and recognition and probe.
% This part is used to get data regarding subjects' rankings and eye responses of other stimuli
% for the other experiments, such as faces / fractals.

% Ask if you want to use Eye Tracker (DIALOG BOX)
if use_eyetracker_general==1
    use_eyetracker = 0;
    ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'ScaleRanking', 'Yes', 'No', 'No');
    if strcmp(ask_if_want_eyetracker, 'Yes')
        use_eyetracker = 1;
    end
end

skip_synctest = 0;

keepTrying = 1;
while keepTrying < 10
    try
        ScaleRanking(subjectID, use_eyetracker,skip_synctest, debugging_mode)
        keepTrying = 10;
    catch
        sca;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - ScaleRanking!');
    end
end


questdlg('Thank you!!! Please call the experimenter','This part is over','Continue','Continue');

keepTrying = 1;
while keepTrying < 10
    try
        ScaleRanking_familiarity(subjectID,0,skip_synctest, debugging_mode)
        keepTrying = 10;
    catch
        sca;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - ScaleRanking!');
    end
end

questdlg('Thank you!!! Please call the experimenter','This part is over','Continue','Continue');


%======================================
%%   'Part 6 - recognition demo & full'
%======================================
% if use_eyetracker_general==1
%     use_eyetracker = 0;
%     ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'recognition demo', 'Yes', 'No', 'No');
%     if strcmp(ask_if_want_eyetracker, 'Yes')
%         use_eyetracker = 1;
%     end
% end

use_eyetracker = 0;

keepTrying = 1;
while keepTrying < 10
    try
        recognition_confidence_demo(isMRI,mainPath,use_eyetracker)
        keepTrying = 10;
    catch
        sca;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - RECOGNITION DEMO!');
    end
end

set(groot,'defaultUicontrolFontSize', 16)
demo_again = questdlg('Do you want more practice?','recognition demo','Yes','No', 'Yes');
if strcmp(demo_again, 'No')
    demo_again = 0;
else
    demo_again = 1;
end

if demo_again==1
    keepTrying = 1;
    while keepTrying < 10
        try
            recognition_confidence_demo(isMRI,mainPath,use_eyetracker)
            keepTrying = 10;
        catch
            sca;
            keepTrying = keepTrying + 1;
            disp('CODE HAD CRASHED - RECOGNITION DEMO!');
        end
    end
end

if use_eyetracker_general==1
    use_eyetracker = 0;
    ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'recognition', 'Yes', 'No', 'No');
    if strcmp(ask_if_want_eyetracker, 'Yes')
        use_eyetracker = 1;
    end
end

keepTrying = 1;
while keepTrying < 10
    try
        recognition_confidence(subjectID, isMRI, mainPath, order, sessionNum, use_eyetracker, debugging_mode);
        keepTrying = 10;
    catch
        sca;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - RECOGNITION!');
    end
end


%=================================
%%   'Part 7 - probe_demo & probe'
%=================================
numBlocks = 1; % Define how many blocks for the probe session. Each block includes all comparisons, one time each.
numRunsPerBlock = 2; % Define the required number of runs per block

% probe demo
% - - - - - -

% if use_eyetracker_general==1
%     use_eyetracker = 0;
%     ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'probe demo', 'Yes', 'No', 'No');
%     if strcmp(ask_if_want_eyetracker, 'Yes')
%         use_eyetracker = 1;
%     end
% end

use_eyetracker = 0;

keepTrying = 1;
while keepTrying < 10
    try
        probeDemo(subjectID, order, mainPath, isMRI, sessionNum, use_eyetracker);
        keepTrying = 10;
    catch
        sca;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - PROBE DEMO!');
    end
end

set(groot,'defaultUicontrolFontSize', 16)
demo_again = questdlg('Do you want more practice?','probe demo','Yes','No','Yes');
if strcmp(demo_again, 'No')
    demo_again = 0;
else
    demo_again = 1;
end

if demo_again==1
    keepTrying = 1;
    while keepTrying < 10
        try
            probeDemo(subjectID, order, mainPath, isMRI, sessionNum, use_eyetracker);
            keepTrying = 10;
        catch
            sca;
            keepTrying = keepTrying + 1;
            disp('CODE HAD CRASHED - PROBE DEMO!');
        end
    end
end

% probe full
% - - - - - -

% Run blocks. Before each block, stimuli need to be organized in a list and
% divided to the required number of runs

if use_eyetracker_general==1
    use_eyetracker = 0;
    ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'probe', 'Yes', 'No', 'No');
    if strcmp(ask_if_want_eyetracker, 'Yes')
        use_eyetracker = 1;
    end
end

for ind = 1:numBlocks
    block = (sessionNum-1)*numBlocks + ind;
    % Organize the stimuli for the probe
    % ===================================
    [trialsPerRun] = organizeProbe(subjectID, order, mainPath, block, numRunsPerBlock);
    for numRun = 1:numRunsPerBlock
        keepTrying = 1;
        while keepTrying < 10
            try
                probe(subjectID, order, mainPath, isMRI, sessionNum, block, numRun, numRunsPerBlock, trialsPerRun, numBlocks, use_eyetracker, debugging_mode);
                keepTrying=10;
            catch
                sca;
                keepTrying = keepTrying+1;
                disp('CODE HAD CRASHED - PROBE!');
            end
        end
        questdlg('Please call the experimenter.','','OK','OK');
    end % end for numRun = 1:numRunsPerBlock
end % end for block = 1:numBlocks


%=================================
%%   'Part 8 - familiarity'
%=================================

% % Ask if you want to use Eye Tracker (DIALOG BOX)
% if use_eyetracker_general==1
%     use_eyetracker = 0;
%     ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'Familiarity demo', 'Yes', 'No', 'No');
%     if strcmp(ask_if_want_eyetracker, 'Yes')
%         use_eyetracker = 1;
%     end
% end

% demo
keepTrying = 1;
while keepTrying < 10
    try
        familiarity(subjectID, sessionNum, 'demo', 0);
        keepTrying = 10;
    catch
        sca;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - familiarity demo!');
    end
end

% % Ask if you want to use Eye Tracker (DIALOG BOX)
% if use_eyetracker_general==1
%     use_eyetracker = 0;
%     ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'Familiarity', 'Yes', 'No', 'No');
%     if strcmp(ask_if_want_eyetracker, 'Yes')
%         use_eyetracker = 1;
%     end
% end


% full
keepTrying = 1;
while keepTrying < 10
    try
        familiarity(subjectID, sessionNum, 'familiarity', 0, debugging_mode);
        keepTrying = 10;
    catch
        sca;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - familiarity!');
    end
end


%==========================================================
%%   'post-task'
%==========================================================

% The experimenter needs to press '123' for the results to be shown
questdlg('Thank you!!! Please call the experimenter to finish the experiment.','','OK','OK');
%msgbox('Thank you!!! Please call the experimenter to finish the experiment.');

code_to_continue = 0;
while code_to_continue ~= 123
    code_to_continue = inputdlg('Enter the code when you are ready for the BDM lottery');
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

end % end function
