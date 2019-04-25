function run_bmem_snacks2_session1()
% function run_bmem_snacks2_session1()
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ================= by Rotem Botvinik Nezer March 2017 ===============
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% This function runs all the parts of session1 of the bmem_snacks2 experiment -session1
% BDM, training demo, training, personal details, BDM resolve
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
% % %   'disp_resolveBDM'
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
experimentName = 'bmem_snacks2'; % hard-coded since true for all subjects
sessionNum = 1; % this is the code for day1 - the first session
isMRI = 0; % this is a behavioral study. change if switch to MRI...
use_eyetracker_general = 1; % in this version we use eye-tracking
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
%fid_crash = fopen([outputPath '/' subjectID '_crashingLogs' num2str(sessionNum) '_' timestamp '.txt'], 'a');

% % =========================================================================
% % PART 1: BDM (python)
% % =========================================================================
do_bdm = 'Yes';
subject_bdm_files = dir([outputPath '/' subjectID '*BDM*']);
if ~isempty(subject_bdm_files)
    do_bdm = questdlg ('Do you want to do the BDM?', 'BDM demo and full', 'Yes', 'No', 'Yes');
end
if strcmp(do_bdm, 'Yes')
    system(['/usr/local/bin/python2.7 ' mainPath '/BDM_demo.py ' subjectID]);
    % % % % % % % % % %
    system(['/usr/local/bin/python2.7 ' mainPath '/BDM.py ' subjectID]);
end
%
% =========================================================================
% Sort stimuli according to the BDM ranking (BDM = PART 1)
% =========================================================================
sort_BDM(subjectID,order,outputPath);

% =========================================================================
% Training (including demo) - (Training = PART 2)
% =========================================================================

% Set number of runs
% -------------------------------------------------
total_num_runs_training = 16;
num_parts_of_training = 2; % we use eyetracking and want to calibrate and validat it once at the middle of training
parts = linspace(1,total_num_runs_training,num_parts_of_training+1);
runInd = ceil(parts(1:end-1));
last_run = [ceil(parts(2:end-1)-1) parts(end)];

% for debugging:
% total_num_runs_training = 4;

%runInd = 1; % The first run to run with the training function (should be changes if for example the first runs are with eye-tracking)

%% Ask if you want to use Eye Tracker (DIALOG BOX)
% =========================================================================
if use_eyetracker_general==1
    use_eyetracker = 0;
    ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'training demo', 'Yes', 'No', 'No');
    if strcmp(ask_if_want_eyetracker, 'Yes')
        use_eyetracker = 1;
    end
end

crashedDemoTraining = 0;
keepTrying = 1;
while keepTrying < 10
    try
        trainingDemo(subjectID,order,mainPath,isMRI,use_eyetracker);
        keepTrying = 10;
    catch
        sca;
        crashedDemoTraining = crashedDemoTraining + 1;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - TRAINING DEMO!');
    end
end
%fprintf(fid_crash,'training demo crashed:\t %d\n', crashedDemoTraining);

set(groot,'defaultUicontrolFontSize', 16)
demo_again = questdlg('Do you want more practice?','','Yes','No','No');
if strcmp(demo_again, 'No')
    demo_again = 0;
else
    demo_again = 1;
end

if demo_again == 1
    crashedDemoTraining = 0;
    keepTrying = 1;
    while keepTrying < 10
        try
            trainingDemo(subjectID,order,mainPath,isMRI,use_eyetracker);
            keepTrying = 10;
        catch
            sca;
            crashedDemoTraining = crashedDemoTraining + 1;
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
    crashedTraining = 0;
    keepTrying = 1;
    while keepTrying < 10
        try
            [Ladder1IN,Ladder2IN] = training(subjectID,order,mainPath,isMRI,runInd(part_num),total_num_runs_training,Ladder1IN,Ladder2IN,use_eyetracker,last_run(part_num));
            keepTrying = 10;
        catch
            sca;
            crashedTraining = crashedTraining + 1;
            keepTrying = keepTrying + 1;
            disp('CODE HAD CRASHED - TRAINING!');
        end
    end
    %fprintf(fid_crash,'training crashed:\t %d\n', crashedTraining);
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
%%   'post-task'
%==========================================================

% The experimenter needs to press '123' for the results to be shown

%questdlg('Thank you!!! Please call the experimenter to finish the experiment.','','OK','OK');

code_to_continue = 0;
while code_to_continue ~= 123
    code_to_continue = inputdlg('Enter the code when you are ready for the BDM lottery');
    code_to_continue = cell2mat(code_to_continue);
    code_to_continue = str2double(code_to_continue);
end

disp_resolveBDM(subjectID, mainPath, sessionNum);

%fclose(fid_crash);

end % end function
