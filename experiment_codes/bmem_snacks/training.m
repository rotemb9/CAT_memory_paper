function [] = training(subjectID,order,mainPath,isMRI,runInd,total_num_runs_training,Ladder1IN,Ladder2IN)

% function [] = training_Israel(subjectID,order,mainPath,isMRI,runInd,total_num_runs_training_trainingLadder1IN,Ladder2IN)

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ======================= by Rotem Botvinik Nezer May 2015 ======================
% ============== edited by Rotem Botvinik Nezer on June 2017 ==============
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function runs the cue-approach training session,
% in which the items are shown on the screen while some of them (GO items) are
% paired with a beep. The subject should press a predefined button (b) as fast
% as possible after hearing the beep.
% This session is composed of total_num_runs_training number of runs.
% After two runs there is a short break. If the subject was bad (less than
% %50 of in-time button pressing out of go trials in these two runs) there
% is a request to press faster (a feedback just for keeping the subjects
% aware if their responses shows they are not).

% This version is for training only 40 items!!!


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_trainingstim.txt' ---> The file for training 40 items,
%   created by the function sort_BDM_Israel
%   '/Onset_files/train_onset_' num2str(r(1)) '.mat''  where r=1-4
%   all the contents of 'stim/' food images
%   'Misc/soundfile.mat'
%   'CenterText.m'


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'training_run_' sprintf('%02d',runNum) '_' timestamp '.txt'
%


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID = 'bmem_snacks_999';
% subjectID = 'bmem_snacks_998'; % to test both order 1 and 2
% mainPath = pwd;
% isMRI = 0;
% runInd = 1;
% total_num_runs_training = 4;
% Ladder1IN = 750;
% Ladder2IN = 750;


tic

rng shuffle

%---------------------------------------------------------------
%%   'GLOBAL VARIABLES'
%---------------------------------------------------------------

outputPath = [mainPath '/Output'];

% about timing
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];

% about ladders
Step = 50;

% about timing
image_duration = 1; %because stim duration is 1.5 secs in opt_stop
baseline_fixation = 1;
afterrunfixation = 1;

% ladders
if nargin < 8
    Ladder1IN = 750;
    Ladder2IN = 750;
end

% -----------------------------------------------
%% Load Instructions
% -----------------------------------------------

% Load Hebrew instructions image files
Instructions = dir([mainPath '/Instructions/*training.JPG' ]);
Instructions_image = imread([mainPath '/Instructions/' Instructions(1).name]);

% -----------------------------------------------
%% 'INITIALIZE SCREEN'
%---------------------------------------------------------------

Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
screennum = max(Screen('Screens'));

pixelSize=32;
% [w] = Screen('OpenWindow',screennum,[],[0 0 640 480],pixelSize);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

%   colors
% - - - - - -
black = BlackIndex(w); % Should equal 0.
white = WhiteIndex(w); % Should equal 255.
Green = [0 255 0];

Screen('FillRect', w, black);
Screen('Flip', w);

%   text
% - - - - - -
theFont = 'Arial';
Screen('TextSize',w,36);
Screen('TextFont',w,theFont);
Screen('TextColor',w,white);


WaitSecs(1);
HideCursor;

% -------------------------------------------------------
%% 'Sound settings'
%%---------------------------------------------------------------

load('Misc/soundfile.mat');

wave = sin(1:0.25:1000);
freq = 22254;
nrchannels = size(wave,1);

deviceID = -1;
% Audio = audioplayer(wave,freq);

reqlatencyclass = 2; % class 2 empirically the best, 3 & 4 == 2

% Initialize driver, request low-latency preinit:
InitializePsychSound(1);

% % Open audio device for low-latency output:
pahandle = PsychPortAudio('Open', deviceID, [], reqlatencyclass, freq, nrchannels);
PsychPortAudio('RunMode', pahandle, 1);

%Play the sound
% play(Audio);
PsychPortAudio('FillBuffer', pahandle, wave);
PsychPortAudio('Start', pahandle, 1, 0, 0);
WaitSecs(1);

% Close the sound and open a new port for the next sound with low latency

% % PsychPortAudio('Stop', pahandle);
PsychPortAudio('Close', pahandle);
pahandle = PsychPortAudio('Open', deviceID, [], reqlatencyclass, freq, nrchannels);
PsychPortAudio('RunMode', pahandle, 1);
PsychPortAudio('FillBuffer', pahandle, wave);

%%---------------------------------------------------------------
%%  'FEEDBACK VARIABLES'
%%---------------------------------------------------------------

if isMRI == 1
    %     trigger = KbName('t');
    blue = KbName('b');
    yellow = KbName('y');
    %     green = KbName('g');
    %     red = KbName('r');
    %     LEFT = [98 5 10];   % blue (5) green (10)
    %     RIGHT = [121 28 21]; % yellow (28) red (21)
else
    BUTTON = 98; %[197];  %<
    %RIGHT = [110]; %[198]; %>
end; % end if isMRI == 1

%---------------------------------------------------------------
%%   'PRE-TRIAL DATA ORGANIZATION'
%---------------------------------------------------------------

%   'Reading in the sorted BDM list - defines which items will be GO/NOGO'
% - - - - - - - - - - - - - - -
file = dir([outputPath '/' subjectID '_stopGoList_trainingstim.txt']);
fid = fopen([outputPath '/' sprintf(file(length(file)).name)]);
vars = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sortbdm
fclose(fid);

%---------------------------------------------------------------
%% 'Display Main Instructions'
%---------------------------------------------------------------
Screen('TextSize', w, 24); %Set textsize

if isMRI
    Screen('PutImage',w,Instructions_image);
    Screen(w,'Flip');

    noresp = 1;
    while noresp,
        [keyIsDown,~,~] = KbCheck; %(experimenter_device);
        if keyIsDown && noresp,
            noresp = 0;
        end;
    end;
    WaitSecs(0.001);
    CenterText(w,'Please focus on the food item onscreen when you hear the sound.',white,0,-220);
    CenterText(w,'You will receive a bonus for looking at the images.',white,0,-170);
    CenterText(w,'Waiting for trigger...Get READY....',white, 0, -50);
    Screen('Flip',w);
    escapeKey = KbName('t');
    while 1
        [keyIsDown,secs,keyCode] = KbCheck(-1); %#ok<ASGLU>
        if keyIsDown && keyCode(escapeKey)
            break;
        end
    end
    
    KbQueueCreate(-1,keylist);%starting the response cue after the trigger
    
else
    Screen('PutImage',w,Instructions_image);
    Screen(w,'Flip');

    noresp = 1;
    while noresp,
        [keyIsDown,~,~] = KbCheck; %(experimenter_device);
        if keyIsDown && noresp,
            noresp = 0;
        end;
    end;
    WaitSecs(0.001);
end

DisableKeysForKbCheck(KbName('t')); % So trigger is no longer detected

KbQueueCreate;



%---------------------------------------------------------------
%%  'TRIAL PRESENTATION'
%---------------------------------------------------------------
%
%   trial_type definitions:
% - - - - - - - - - - - -
% 11 = High-Value GO
% 12 = High-Value NOGO
% 22 = Low-Value GO
% 24 = Low-Value NOGO
% 0 = not trained

% Setting the size of the variables for the loop
%---------------------------
shuff_names = cell(1,total_num_runs_training);
shuff_ind = cell(1,total_num_runs_training);
bidIndex = cell(1,total_num_runs_training);
shuff_bidIndex = cell(1,total_num_runs_training);
BDMtrialIndex = cell(1,total_num_runs_training);
shuff_BDMtrialIndex = cell(1,total_num_runs_training);
trialType = cell(1,total_num_runs_training);
shuff_trialType = cell(1,total_num_runs_training);
Audio_time = cell(1,total_num_runs_training);
respTime = cell(1,total_num_runs_training);
respInTime = cell(1,total_num_runs_training);
keyPressed = cell(1,total_num_runs_training);
Ladder1 = cell(1,total_num_runs_training);
Ladder2 = cell(1,total_num_runs_training);
actual_onset_time = cell(1,total_num_runs_training);
fix_time = cell(1,total_num_runs_training);
fixcrosstime = cell(1,total_num_runs_training);
Ladder1end = cell(1,total_num_runs_training);
Ladder2end = cell(1,total_num_runs_training);
correct = cell(1,total_num_runs_training);
numGoTrials = zeros(1, total_num_runs_training);
mean_RT = cell(1,total_num_runs_training);
bidValues = cell(1,total_num_runs_training);
shuff_bidValues = cell(1,total_num_runs_training);


anchor = GetSecs ; % (before baseline fixation) ;

%  for runNum = runInd:runInd+3; % for debugging 4 runs starting with runInd
for runNum = runInd:total_num_runs_training %this for loop allows all runs to be completed
    
    KbQueueFlush;
    
    %   'load onsets'
    %---------------------------
    r = Shuffle(1:4);
    load(['Onset_files/train_onset_' num2str(r(1)) '.mat']);
    
    
    %   'Write output file header'
    %---------------------------------------------------------------
    c = clock;
    hr = sprintf('%02d', c(4));
    minutes = sprintf('%02d', c(5));
    timestamp = [date,'_',hr,'h',minutes,'m'];
    
    
    fid1 = fopen([outputPath '/' subjectID '_training40_run_' sprintf('%02d',runNum) '_' timestamp '.txt'], 'a');
    fprintf(fid1,'subjectID\torder\trunNum\titemName\tonsetTime\tshuff_trialType\tRT\trespInTime\tAudioTime\tresponse\tfixationTime\tladder1\tladder2\tbidIndex\tBDMtrialIndex\tbidValue\n'); %write the header line
    
    
    %   'pre-trial fixation'
    %---------------------------
    
    firstOrSecond = mod(runNum,2);
    
    switch firstOrSecond
        case 1
            prebaseline = GetSecs;
            % baseline fixation - currently 10 seconds = 4*Volumes (2.5 TR)
            while GetSecs < prebaseline + baseline_fixation
                CenterText(w,'+', white,0,0);
                Screen('TextSize',w, 60);
                Screen(w,'Flip');
            end
        case 0
            prebaseline = GetSecs;
            % baseline fixation - currently 10 seconds = 4*Volumes (2.5 TR)
            while GetSecs < prebaseline + afterrunfixation
                CenterText(w,'+', white,0,0);
                Screen('TextSize',w, 60);
                Screen(w,'Flip');
            end
    end
    
    
    %   Reading everying from the sorted StopGo file - vars has everything
    %---------------------------
    [shuff_names{runNum},shuff_ind{runNum}] = Shuffle(vars{1});
    
    trialType{runNum} = vars{2};
    shuff_trialType{runNum} = trialType{runNum}(shuff_ind{runNum});
    
    bidIndex{runNum} = vars{3};
    shuff_bidIndex{runNum} = bidIndex{runNum}(shuff_ind{runNum});
    
    bidValues{runNum} = vars{4};
    shuff_bidValues{runNum} = bidValues{runNum}(shuff_ind{runNum});
    
    BDMtrialIndex{runNum} = vars{5};
    shuff_BDMtrialIndex{runNum} = BDMtrialIndex{runNum}(shuff_ind{runNum});
    
    %	pre-allocating matrices
    %---------------------------
    Audio_time{runNum}(1:length(shuff_trialType{runNum}),1) = 999;
    respTime{runNum}(1:length(shuff_trialType{runNum}),1) = 999;
    respInTime{runNum}(1:length(shuff_trialType{runNum}),1) = 999;
    keyPressed{runNum}(1:length(shuff_trialType{runNum}),1) = 999;
    
    %   reading in images
    %---------------------------
    food_items = cell(1, length(shuff_names{runNum}));
    for i = 1:length(shuff_names{runNum})
        food_items{i} = imread(sprintf('stim/%s',shuff_names{runNum}{i}));
    end
    
    %   Read in info about ladders
    % - - - - - - - - - - - - - - -
    
    if runNum == 1
        Ladder1{runNum}(1,1) = 750;
        Ladder2{runNum}(1,1) = 750;
    elseif runNum == runInd
        Ladder1{runNum}(1,1) = Ladder1IN;
        Ladder2{runNum}(1,1) = Ladder2IN;
    else % runNum > 1
        Ladder1{runNum}(1,1) = Ladder1end{runNum-1};
        Ladder2{runNum}(1,1) = Ladder2end{runNum-1};
        %         tmp2 = dir([outputPath '/' subjectID '_ladders_run_' num2str(runNum - 1) '.txt']);
        %         fid2 = fopen([outputPath '/' tmp2(length(tmp2)).name]);
        %         ladders = textscan(fid2, '%d %d %d', 'Headerlines',1);
        %         fclose(fid2);
        %         Ladder1{runNum}(1,1) = ladders{1};
        %         Ladder2{runNum}(1,1) = ladders{2};
        %         %         runInd = ladders{3}+1;
    end
    
    
    %   Loop through all trials in a run
    %---------------------------
    runStartTime = GetSecs - anchor;
    
%         for trialNum = 1:6; % shorter version for debugging
    for trialNum = 1:length(shuff_trialType{runNum})   % To cover all the items in one run.
        Screen('PutImage',w,food_items{trialNum});
        Screen('Flip',w,anchor+onsets(trialNum)+runStartTime); % display images according to Onset times
        image_start_time = GetSecs;
        actual_onset_time{runNum}(trialNum,1) = image_start_time - anchor;
        
        noresp = 1;
        notone = 1;
        KbQueueFlush; % this is important for the queue to be empty for the next trial (to prevent RT with negative values)
        KbQueueStart;
        
        %---------------------------------------------------
        %% 'EVALUATE RESPONSE & ADJUST LADDER ACCORDINGLY'
        %---------------------------------------------------
        while (GetSecs-image_start_time < image_duration)
            
            %High-Valued BEEP items
            %---------------------------
            if  shuff_trialType{runNum}(trialNum) == 11 && (GetSecs - image_start_time >= Ladder1{runNum}(length(Ladder1{runNum}),1)/1000) && notone % shuff_trialType contains the information if a certain image is a GO/NOGO trial
                % Beep!
                % %                 PsychPortAudio('FillBuffer', pahandle, wave);
                PsychPortAudio('Start', pahandle, 1, 0, 0);
                %                 play(Audio);
                notone = 0;
                Audio_time{runNum}(trialNum,1) = GetSecs-image_start_time;
                
                %   look for response
                [pressed, firstPress, ~, ~, ~] = KbQueueCheck;
                if pressed && noresp
                    firstKeyPressed = find(firstPress==min(firstPress(firstPress>0)));
                    if length(firstKeyPressed)>=2 % In case two keys were pressed together on the exact micro-sec
                        firstKeyPressed = firstKeyPressed(1);
                    end
                    respTime{runNum}(trialNum,1) = firstPress(firstKeyPressed)-image_start_time;
                    %                     findfirstPress = find(firstPress);
                    %                     respTime{runNum}(trialNum,1) = firstPress(findfirstPress(1))-image_start_time;
                    tmp = KbName(firstKeyPressed);
                    if ischar(tmp) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed{runnum} to be a char, so this converts it and takes the first key pressed
                        tmp = char(tmp);
                    end
                    keyPressed{runNum}(trialNum,1) = tmp(1);
                    
                    
                    % different response types in scanner or in testing room
                    if isMRI == 1
                        if keyPressed{runNum}(trialNum,1) == blue || keyPressed{runNum}(trialNum,1) == yellow
                            noresp = 0;
                            
                            if respTime{runNum}(trialNum,1) < Ladder1{runNum}(length(Ladder1{runNum}),1)/1000
                                respInTime{runNum}(trialNum,1) = 11; %was a GO trial with HV item but responded before SS
                            else
                                respInTime{runNum}(trialNum,1)= 110; %was a Go trial with HV item but responded after SS within 1000 msec
                            end
                        end
                    else
                        
                        if keyPressed{runNum}(trialNum,1) == BUTTON %| keyPressed{runnum}(trialnum,1)==RIGHT
                            noresp = 0;
                            
                            if respTime{runNum}(trialNum,1) < Ladder1{runNum}(length(Ladder1{runNum}),1)/1000
                                respInTime{runNum}(trialNum,1) = 11; %was a GO trial with HV item but responded before SS
                            else
                                respInTime{runNum}(trialNum,1) = 110; %was a Go trial with HV item and responded after SS within 1000 msec - good trial
                            end
                        end
                    end % if isMRI == 1
                    
                end
                
                %Low-Valued BEEP items
                %---------------------------
            elseif  shuff_trialType{runNum}(trialNum) == 22 && (GetSecs - image_start_time >= Ladder2{runNum}(length(Ladder2{runNum}),1)/1000) && notone %shuff_trialType contains the information if a certain image is a GO/NOGO trial
                
                % Beep!
                %                 PsychPortAudio('FillBuffer', pahandle, wave);
                PsychPortAudio('Start', pahandle, 1, 0, 0);
                %                 play(Audio);
                notone = 0;
                Audio_time{runNum}(trialNum,1) = GetSecs-image_start_time;
                
                
                %   look for response
                [pressed, firstPress, ~, ~, ~] = KbQueueCheck;
                if pressed && noresp
                    firstKeyPressed = find(firstPress==min(firstPress(firstPress>0)));
                    if length(firstKeyPressed)>=2 % In case two keys were pressed together on the exact micro-sec
                        firstKeyPressed = firstKeyPressed(1);
                    end
                    respTime{runNum}(trialNum,1) = firstPress(firstKeyPressed)-image_start_time;
                    %                     findfirstPress = find(firstPress);
                    %                     respTime{runNum}(trialNum,1) = firstPress(findfirstPress(1))-image_start_time;
                    tmp = KbName(firstKeyPressed);
                    if ischar(tmp) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed{runnum} to be a char, so this converts it and takes the first key pressed
                        tmp = char(tmp);
                    end
                    keyPressed{runNum}(trialNum,1) = tmp(1);
                    
                    %   different response types in scanner or in testing room
                    if isMRI == 1
                        if keyPressed{runNum}(trialNum,1) == blue || keyPressed{runNum}(trialNum,1) == yellow
                            noresp = 0;
                            if respTime{runNum}(trialNum,1) < Ladder2{runNum}(length(Ladder2{runNum}),1)/1000
                                respInTime{runNum}(trialNum,1) = 22; %was a GO trial with LV item but responded before SS
                            else
                                respInTime{runNum}(trialNum,1) = 220; %was a Go trial with LV item but responded after SS within 1000 msec
                            end
                        end
                    else
                        if keyPressed{runNum}(trialNum,1) == BUTTON %| keyPressed{runnum}(trialnum,1)==RIGHT
                            noresp = 0;
                            if respTime{runNum}(trialNum,1) < Ladder2{runNum}(length(Ladder2{runNum}),1)/1000
                                respInTime{runNum}(trialNum,1) = 22;  %was a GO trial with LV item but responded before SS
                            else
                                respInTime{runNum}(trialNum,1) = 220; %was a Go trial with LV item and responded after SS within 1000 msec - good trial
                            end
                        end
                    end % if isMRI == 1
                    
                end % end if pressed && noresp
                
                %No-BEEP
                %---------------------------
            elseif   mod(shuff_trialType{runNum}(trialNum),11) ~= 0 && noresp % these will now be the NOGO trials
                
                %   look for response
                [pressed, firstPress, ~, ~, ~] = KbQueueCheck;
                if pressed && noresp
                    firstKeyPressed = find(firstPress==min(firstPress(firstPress>0)));
                    if length(firstKeyPressed)>=2 % In case two keys were pressed together on the exact micro-sec
                        firstKeyPressed = firstKeyPressed(1);
                    end
                    respTime{runNum}(trialNum,1) = firstPress(firstKeyPressed)-image_start_time;
                    %                     findfirstPress = find(firstPress);
                    %                     respTime{runNum}(trialNum,1) = firstPress(findfirstPress(1))-image_start_time;
                    tmp = KbName(firstKeyPressed);
                    if ischar(tmp) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed{runnum} to be a char, so this converts it and takes the first key pressed
                        tmp = char(tmp);
                    end
                    keyPressed{runNum}(trialNum,1) = tmp(1);
                    
                    % different response types in scanner or in testing room
                    if isMRI == 1
                        if keyPressed{runNum}(trialNum,1) == blue || keyPressed{runNum}(trialNum,1) == yellow
                            noresp = 0;
                            if shuff_trialType{runNum}(trialNum) == 12
                                respInTime{runNum}(trialNum,1) = 12; % a stop trial but responded within 1000 msec HV item - not good but don't do anything
                            else
                                respInTime{runNum}(trialNum,1) = 24; % a stop trial but responded within 1000 msec LV item - not good but don't do anything
                            end
                        end
                    else
                        if keyPressed{runNum}(trialNum,1) == BUTTON %| keyPressed{runnum}(trialnum,1)==RIGHT
                            noresp = 0;
                            if shuff_trialType{runNum}(trialNum) == 12
                                respInTime{runNum}(trialNum,1) = 12; %% a stop trial but responded within 1000 msec HV item - not good but don't do anything
                            else
                                respInTime{runNum}(trialNum,1) = 24; %% a stop trial but responded within 1000 msec LV item - not good but don't do anything
                            end
                        end
                    end % end if test+comp == 1
                    
                end % end if pressed && noresp
            end %evaluate trial_type
            
        end %%% End big while waiting for response within 1000 msec
        
        
        %         %   Close the Audio port and open a new one
        %         %------------------------------------------
        % %                 PsychPortAudio('Stop', pahandle);
        PsychPortAudio('Close', pahandle);
        pahandle = PsychPortAudio('Open', deviceID, [], reqlatencyclass, freq, nrchannels);
        PsychPortAudio('RunMode', pahandle, 1);
        PsychPortAudio('FillBuffer', pahandle, wave);
        %
        
        %   Show fixation
        %---------------------------
        CenterText(w,'+', white,0,0);
        Screen('TextSize',w, 60);
        Screen(w,'Flip', image_start_time+1);
        fix_time{runNum}(trialNum,1) = GetSecs ;
        fixcrosstime{runNum} = GetSecs;
        
        
        if noresp == 1
            %---------------------------
            % these are additional 500msec to monitor responses
            
            while (GetSecs-fix_time{runNum}(trialNum,1) < 0.5)
                
                %   look for response
                [pressed, firstPress, ~, ~, ~] = KbQueueCheck;
                if pressed && noresp
                    firstKeyPressed = find(firstPress==min(firstPress(firstPress>0)));
                    if length(firstKeyPressed)>=2 % In case two keys were pressed together on the exact micro-sec
                        firstKeyPressed = firstKeyPressed(1);
                    end
                    respTime{runNum}(trialNum,1) = firstPress(firstKeyPressed)-image_start_time;
                    %                     findfirstPress = find(firstPress);
                    %                     respTime{runNum}(trialNum,1) = firstPress(findfirstPress(1))-image_start_time;
                    tmp = KbName(firstKeyPressed);
                    if ischar(tmp) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed{runnum} to be a char, so this converts it and takes the first key pressed
                        tmp = char(tmp);
                    end
                    keyPressed{runNum}(trialNum,1) = tmp(1);
                    
                    if isMRI == 1
                        if keyPressed{runNum}(trialNum,1) == blue || keyPressed{runNum}(trialNum,1) == yellow
                            noresp = 0;
                            switch shuff_trialType{runNum}(trialNum)
                                case 11
                                    if respTime{runNum}(trialNum,1) >= 1
                                        respInTime{runNum}(trialNum,1) = 1100; % a Go trial and  responded after 1000msec  HV item - make it easier decrease SSD
                                    elseif respTime{runNum}(trialNum,1) < 1
                                        respInTime{runNum}(trialNum,1) = 110;
                                    end
                                case 22
                                    if respTime{runNum}(trialNum,1) >= 1
                                        respInTime{runNum}(trialNum,1) = 2200; % a Go trial and  responded after 1000msec  HV item - make it easier decrease SSD
                                    elseif respTime{runNum}(trialNum,1) < 1
                                        respInTime{runNum}(trialNum,1) = 220;
                                    end
                                case 12
                                    respInTime{runNum}(trialNum,1) = 12; % a stop trial and responded after 1000 msec  HV item - don't touch
                                case 24
                                    respInTime{runNum}(trialNum,1) = 24; % % a stop trial and  responded after 1000 msec HV item - don't touch
                            end
                        end
                        
                    else
                        
                        if keyPressed{runNum}(trialNum,1) == BUTTON % | keyPressed{runnum}(trialnum,1)==RIGHT
                            noresp = 0;
                            switch shuff_trialType{runNum}(trialNum)
                                case 11
                                    if respTime{runNum}(trialNum,1) >= 1
                                        respInTime{runNum}(trialNum,1) = 1100;% a Go trial and responded after 1000msec  HV item  - make it easier decrease SSD
                                    elseif respTime{runNum}(trialNum,1) < 1
                                        respInTime{runNum}(trialNum,1) = 110;% a Go trial and responded before 1000msec  HV item  -  make it harder increase SSD/3
                                    end
                                case 22
                                    if respTime{runNum}(trialNum,1) > 1
                                        respInTime{runNum}(trialNum,1) = 2200;% a Go trial and responded after 1000msec  LV item - - make it easier decrease SSD
                                    elseif respTime{runNum}(trialNum,1) < 1
                                        respInTime{runNum}(trialNum,1) = 220;% a Go trial and responded before 1000msec  LV item - - make it harder increase SSD/3
                                    end
                                case 12
                                    respInTime{runNum}(trialNum,1) = 12;% a NOGO trial and didnt respond on time HV item - don't touch
                                case 24
                                    respInTime{runNum}(trialNum,1) = 24;% a NOGO trial and didnt respond on time LV item - don't touch
                                    
                            end
                        end
                    end % end if isMRI == 1
                end % end if pressed && noresp
            end % End while of additional 500 msec
        else % the subject has already responded during the first 1000 ms
            WaitSecs(0.5);
        end  % end if noresp
        
        %%	This is where its all decided !
        %---------------------------
        if noresp
            switch shuff_trialType{runNum}(trialNum)
                case 11
                    respInTime{runNum}(trialNum,1) = 1; %unsuccessful Go trial HV - didn't press a button at all - trial too hard - need to decrease ladder
                case 22
                    respInTime{runNum}(trialNum,1) = 2; % unsuccessful Go trial LV - didn't press a button at all - trial too hard - need to decrease ladder
                case 12
                    respInTime{runNum}(trialNum,1) = 120; % ok NOGO trial didn't respond after 1500 msec in NOGO trial HV
                case 24
                    respInTime{runNum}(trialNum,1) = 240; % ok NOGO trial didn't respond after 1500 msec in NOGO trial LV
            end
        end
        
        
        switch respInTime{runNum}(trialNum,1)
            case 1 % didn't respond even after 1500 msec on HV GO trial - make it easier decrease SSD by step
                if (Ladder1{runNum}(length(Ladder1{runNum}),1)<0.001)
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1);
                else
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1)-Step;
                end;
                
            case 2 % didn't respond even after 1500 msec on LV GO trial - make it easier decrease SSD by step
                if (Ladder2{runNum}(length(Ladder2{runNum}),1)<0.001)
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1);
                else
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1)-Step;
                end;
                
                
            case 1100 %  responded after 1500 msec on HV GO trial - make it easier decrease SSD by step
                if (Ladder1{runNum}(length(Ladder1{runNum}),1)<0.001)
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1);
                else
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1)-Step;
                end;
                
            case 2200 %  responded after 1500 msec on LV GO trial - make it easier decrease SSD by step
                if (Ladder2{runNum}(length(Ladder2{runNum}),1)<0.001)
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1);
                else
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1)-Step;
                end;
                
                
                
            case 11
                if (Ladder1{runNum}(length(Ladder1{runNum}),1) > 910); %was a GO trial with HV item but responded before SS make it harder - increase SSD by Step/3
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1);
                else
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1)+Step/3;
                end;
                
            case 22
                if (Ladder2{runNum}(length(Ladder2{runNum}),1) > 910); %was a GO trial with LV item but responded before SS make it harder - - increase SSD by Step/3
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1);
                else
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1)+Step/3;
                end;
                
            case 110 % pressed after Go signal but below 1000 - - increase SSD by Step/3 - these are the good trials!
                if (Ladder1{runNum}(length(Ladder1{runNum}),1) > 910);
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1);
                else
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1)+Step/3;
                end;
                
            case 220 % pressed after Go signal but below 1000 - - increase SSD by Step/3 - these are the good trials!
                if (Ladder2{runNum}(length(Ladder2{runNum}),1) > 910);
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1);
                else
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1)+Step/3;
                end;
                
        end % end switch respInTime{runNum}(trialNum,1)
        
        
        %   'Save data'
        %---------------------------
        
        fprintf(fid1,'%s\t %d\t %d\t %s\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %.2f\t %.2f\t %.2f\t %d\t %.2f\t \n', subjectID, order, runNum, shuff_names{runNum}{trialNum}, actual_onset_time{runNum}(trialNum,1), shuff_trialType{runNum}(trialNum), respTime{runNum}(trialNum,1)*1000, respInTime{runNum}(trialNum,1), Audio_time{runNum}(trialNum,1)*1000, keyPressed{runNum}(trialNum,1),   fix_time{runNum}(trialNum,1)-anchor, Ladder1{runNum}(length(Ladder1{runNum})), Ladder2{runNum}(length(Ladder2{runNum})), shuff_bidIndex{runNum}(trialNum,1), shuff_BDMtrialIndex{runNum}(trialNum,1),shuff_bidValues{runNum}(trialNum,1));
             
    end; %	End the big trialNum loop showing all the images in one run.
    
    KbQueueFlush;
    
    Ladder1end{runNum} = Ladder1{runNum}(length(Ladder1{runNum}));
    Ladder2end{runNum} = Ladder2{runNum}(length(Ladder2{runNum}));
    correct{runNum}(1) = 0;
    % Correct trials are when the subject pressed the button on a go trial,
    % either before (11,22) or after (110,220)the beep (but before the
    % image disappeared)
    correct{runNum}(1) = length(find(respInTime{runNum} == 11 | respInTime{runNum} == 110 | respInTime{runNum} == 22 | respInTime{runNum} == 220 ));
    numGoTrials(runNum) = length(find(trialType{runNum} == 11 | trialType{runNum} == 22));
    mean_RT{runNum} = mean(respTime{runNum}(respInTime{runNum} == 110 | respInTime{runNum} == 220));
    
    % afterrun fixation
    % ---------------------------
    postexperiment = GetSecs;
    while GetSecs < postexperiment+afterrunfixation;
        CenterText(w,'+', white,0,0);
        Screen('TextSize',w, 60);
        Screen(w,'Flip');
        
    end
    
    if runNum ~= total_num_runs_training && firstOrSecond == 0 && runNum ~= runInd % if this is not the last run, but it is an even run and there was a run before (runNum~=runInd), display instructions again for the next run
        goodTrials = correct{runNum-1} + correct{runNum};
        goTrials = numGoTrials(runNum-1) + numGoTrials(runNum);
        Screen('TextSize', w, 40); %Set textsize
        if goodTrials < goTrials/2
            %             CenterText(w,strcat(sprintf('You responded on %.2f', ((correct{runnum-1}(1)+correct{runnum}(1)))/20*100), '% of Go trials'), white, 0,-270);
            CenterText(w,sprintf('Please try to respond faster, as fast as you can'), white, 0,-300);
        else
            CenterText(w,sprintf('This is a short break'), white, 0,-300);
            %             CenterText(w,strcat(sprintf('You responded on %.2f', (correct{runNum}(1))/20*100), '% of BEEP trials'), white, 0,-300);
        end % end if correct{runNum}(1) < numGoItems/2
        CenterText(w,sprintf('Press any key to continue to the next run') ,white,0,-100);
        Screen('Flip',w);
        
        noresp = 1;
        while noresp,
            [keyIsDown,~,~] = KbCheck;
            if keyIsDown && noresp,
                noresp = 0;
            end;
        end;
        WaitSecs(0.5);
        
        CenterText(w,'Press the button b when you hear the sound.',white,0,-270);
        CenterText(w,'Press the button as FAST as you can.',white,0,-170);
        Screen('Flip',w);
        
        noresp=1;
        while noresp,
            [keyIsDown,~,~] = KbCheck;
            if keyIsDown && noresp,
                noresp=0;
            end;
        end;
        WaitSecs(0.5);
    end
    
    %   write Ladders info to txt file
    % ------------------------------
    %     fid2 = fopen([outputPath '/' subjectID sprintf('_ladders_run_%d.txt', runNum)], 'w');
    %     fprintf(fid2,'Ladder1\t Ladder2\t runnum\t \n'); %write the header line
    %     fprintf(fid2, '%d\t %d\t %d\t \n', Ladder1{runNum}(length(Ladder1{runNum})), Ladder2{runNum}(length(Ladder2{runNum})),runNum);
    %     fprintf(fid2, '\n');
    %     fclose(fid2);
    
end % End the run loop to go over all the runs

%---------------------------------------------------------------
%%   save data to a .mat file & close out
%---------------------------------------------------------------
% outfile = strcat(outputPath, '/', subjectID,'_training_run', sprintf('%02d',runInd),'_to_run', sprintf('%02d',runNum), '_eyetracking_', timestamp,'.edf');
outfile = strcat(outputPath, '/', subjectID,'_training_run', sprintf('%02d',runInd),'_to_run', sprintf('%02d',runNum),'_', timestamp,'.mat');
% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;
clear food_items Instructions*;

save(outfile);

% Close the audio device:
PsychPortAudio('Close', pahandle);

%   outgoing msg & closing
% ------------------------------
Screen('TextSize',w,36);
Screen('TextFont',w,'Ariel');

CenterText(w,'Great Job. Thank you!',Green, 0,-270);
CenterText(w,'The next part will begin soon', white, 0, -180);
CenterText(w,'Press any key to continue',white, 0,-120);
Screen('Flip',w);

noresp = 1;
while noresp
    [keyIsDown,~,~] = KbCheck;%(experimenter_device);
    if keyIsDown && noresp,
        noresp = 0;
    end;
end;
WaitSecs(0.2);
KbQueueFlush;

Screen('CloseAll');
ShowCursor;


end % end function

