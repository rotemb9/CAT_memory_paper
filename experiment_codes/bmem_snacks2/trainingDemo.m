function [] = trainingDemo(subjectID,order,mainPath,isMRI, use_eyetracker)

% function [] = trainingDemo_Israel(subjectID,order,mainPath,isMRI,use_eyetracker)
%%= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Rotem Botvinik November 2014 ====================
%                 Edited on March 2017 to include eyetracker
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% This function runs the demo of the cue-approach training session,
% in which the demo items are shown on the screen while some of them (GO items) are
% paired with a beep. The subject should press a pre-defined button ('b') as fast
% as possible after hearing the beep.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'stopGoList_allstim_order%d.txt', order
%   '/Onset_files/train_onset_%d' num2str(r(1)) '.mat']  where r=1-4
%    all the contents of '/stim/' food images
%   'Misc/soundfile.mat'
%   'Misc/demo_items.txt'
%   'Misc/demo_go.mat'
%   'CenterText.m'
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'training_demo_' timestamp '.txt'
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID = 'bmem_snacks_999';
% subjectID = 'bmem_snacks_998'; % to test both order 1 and 2
% mainPath = pwd;
% isMRI = 0;

tic

rng shuffle

if nargin < 5
    use_eyetracker = 1;
end

if nargin < 4
    isMRI = 0;
end

%=========================================================================
%% Training task code
%=========================================================================

outputPath = [mainPath '/Output'];
LADDER1IN = 750;
LADDER2IN = 750;

% -----------------------------------------------
%% Load Instructions
% -----------------------------------------------

% Load Hebrew instructions image files
Instructions = dir([mainPath '/Instructions/*training_demo.JPG' ]);
Instructions_image = imread([mainPath '/Instructions/' Instructions(1).name]);

%---------------------------------------------------------------
%% 'INITIALIZE SCREEN'
%---------------------------------------------------------------

Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
screennum = min(Screen('Screens'));

pixelSize = 32;
% [w] = Screen('OpenWindow',screennum,[],[0 0 640 480],pixelSize);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

%   colors
% - - - - - -
black = BlackIndex(w); % Should equal 0.
white = WhiteIndex(w); % Should equal 255.
Green = [0 255 0];

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);

%   text
% - - - - - -
theFont = 'Arial';
Screen('TextSize',w,36);
Screen('TextFont',w,theFont);
Screen('TextColor',w,white);

HideCursor;

%---------------------------------------------------------------
%%   'GLOBAL VARIABLES'
%---------------------------------------------------------------

c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];

%	timing variables
% - - - - - - - - -
Step = 50;
image_duration = 1; %because stim duration is 1.5 secs in opt_stop
baseline_fixation = 1;

runNum = 1;
Ladder1{runNum}(1,1) = LADDER1IN;
Ladder2{runNum}(1,1) = LADDER2IN;

%---------------------------------------------------------------
%%   'PRE-TRIAL DATA ORGANIZATION'
%---------------------------------------------------------------

%   'Reading from the sorted BDM list - defines which items will be GO\NOGO'
% - - - - - - - - - - - - - - -
tmp = dir([outputPath '/' subjectID sprintf('_stopGoList_allstim_order%d.txt', order)]);
fid = fopen([outputPath '/' tmp(length(tmp)).name]);
vars = textscan(fid, '%s %d %d %f %d') ;% these contain everything from the sortbdm
fid_demoNames = fopen('Misc/demo_items.txt');
demoNames = textscan(fid_demoNames,'%s');
demoNames = demoNames{1};
fclose(fid);

%---------------------------------------------------------------
%%% FEEDBACK VARIABLES
%---------------------------------------------------------------
if isMRI == 1
    %     trigger = KbName('t');
    blue = KbName('b');
    yellow = KbName('y');
    %     green = KbName('g');
    %     red = KbName('r');
    %     LEFT = [98 5 10];   %blue (5) green (10)
    %     RIGHT = [121 28 21]; %yellow (28) red (21)
else
    BUTTON = 98; %[197];  %<
    %RIGHT = [110]; %[198]; %>
end;

%-----------------------------------------------------------------
%% Initializing eye tracking system %
%-----------------------------------------------------------------
% use_eyetracker=0; % set to 1/0 to turn on/off eyetracker functions
if use_eyetracker
    dummymode=0;
    
    % STEP 2
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    el=EyelinkInitDefaults(w);
    % Disable key output to Matlab window:
    
    el.backgroundcolour = black;
    el.backgroundcolour = black;
    el.foregroundcolour = white;
    el.msgfontcolour    = white;
    el.imgtitlecolour   = white;
    el.calibrationtargetcolour = el.foregroundcolour;
    EyelinkUpdateDefaults(el);
    
    % STEP 3
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit(dummymode, 1)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        return;
    end;
    
    [~,vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,HREF,AREA');
    
    % open file to record data to
    task = GenFlags.TrainingDemo.str; % change to current task - with GenFlags
    edfFile='TrainDem.edf';
    Eyelink('Openfile', edfFile);
    
    % STEP 4
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
    % do a final check of calibration using driftcorrection
    EyelinkDoDriftCorrection(el);
    
    %     % STEP 5
    %     % start recording eye position
    %     Eyelink('StartRecording');
    %     % record a few samples before we actually start displaying
    %     WaitSecs(0.1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Finish Initialization %
    %%%%%%%%%%%%%%%%%%%%%%%%%
end

WaitSecs(1);
%---------------------------------------------------------------
%% Sound settings
%---------------------------------------------------------------
%%%% Psychportaudio
% respInTime{runNum} = shuff_stop;

load('Misc/soundfile.mat');

wave = sin(1:0.25:1000);

freq = 22254;
nrchannels = size(wave,1);

deviceID = -1;

% Audio = audioplayer(wave,freq);

reqlatencyclass = 2; % class 2 empirically the best, 3 & 4 == 2
% Initialize driver, request low-latency preinit:
InitializePsychSound(1);

% Open audio device for low-latency output:
pahandle = PsychPortAudio('Open', deviceID, [], reqlatencyclass, freq, nrchannels);
PsychPortAudio('RunMode', pahandle, 1);

% Play the sound
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

numRuns = 1;

% Setting the size of the cell matrices:
shuff_names = cell(1,numRuns);
shuff_ind = cell(1,numRuns);
bidIndex = cell(1,numRuns);
shuff_bidIndex = cell(1,numRuns);
itemnameIndex = cell(1,numRuns);
shuff_itemnameIndex = cell(1,numRuns);
trialType = cell(1,numRuns);
Audio_time = cell(1,numRuns);
respTime = cell(1,numRuns);
keyPressed = cell(1,numRuns);
demofood_items = cell(1,length(demoNames));
actual_onset_time = cell(1,numRuns);
respInTime = cell(1,numRuns);
fix_time = cell(1,numRuns);
fixcrosstime = cell(1,numRuns);
Ladder1end = cell(1,numRuns);
Ladder2end = cell(1,numRuns);
correct = cell(1,numRuns);

%---------------------------------------------------------------
%% 'Display Main Instructions'
%---------------------------------------------------------------

if isMRI
    Screen('PutImage',w,Instructions_image);
    Screen(w,'Flip');
    noResponse=1;
    while noResponse,
        [keyIsDown,~,~] = KbCheck; %(experimenter_device);
        if keyIsDown && noResponse,
            noResponse=0;
        end;
    end;
    WaitSecs(0.001);
    Screen('TextSize',w, 40);
    CenterText(w,'Please focus on the food item onscreen when you hear the sound.',white,0,-220);
    CenterText(w,'You will receive a bonus for looking at the images.',white,0,-170);
    CenterText(w,'Waiting for trigger...Get READY....',white, 0, -50);
    
    Screen('Flip',w);
    
    escapeKey = KbName('t');
    while 1
        [keyIsDown,~,keyCode] = KbCheck(-1); % #ok<ASGLU>
        if keyIsDown && keyCode(escapeKey)
            break;
        end
    end
    
else
    Screen('PutImage',w,Instructions_image);
    Screen(w,'Flip');
    noResponse=1;
    while noResponse,
        [keyIsDown,~,~] = KbCheck; %(experimenter_device);
        if keyIsDown && noResponse,
            noResponse=0;
        end;
    end;
    WaitSecs(0.001);
end

DisableKeysForKbCheck(KbName('t')); % So trigger is no longer detected

anchor = GetSecs; % reference point for the task start point (after instructions)

KbQueueCreate;

if use_eyetracker
    % start recording eye position
    %---------------------------
    Eyelink('StartRecording');
    %   Eyelink MSG
    % ---------------------------
    % messages to save on each trial ( trial number, onset and RT)
    WaitSecs(.05);
end

%---------------------------------------------------------------
%%  'TRIAL PRESENTATION
%---------------------------------------------------------------

for runNum = 1:numRuns
    
    %     runStartTime = GetSecs - anchor;
    runStartTime = GetSecs;
    
    if use_eyetracker
        Eyelink('Message', Eventflag(GenFlags.RunStart.str,task,runNum,1,anchor)); % mark start time in file
    end
    
    r = Shuffle(1:4);
    load(['Onset_files/train_onset_' num2str(r(1)) '.mat']);
    %% 'Write output file header'
    %---------------------------------------------------------------
    c = clock;
    hr = sprintf('%02d', c(4));
    minutes = sprintf('%02d', c(5));
    timestamp = [date,'_',hr,'h',minutes,'m'];
    
    
    fid1 = fopen([outputPath '/' subjectID '_training_demo_' timestamp '.txt'], 'a');
    fprintf(fid1,'subjectID\trunNum\titemname\tonsettime\ttrialtype\tRT\trespInTime\tAudiotime\tresponse\tfixationtime\tladder1\tladder2\n'); %write the header line
    
    prebaseline = GetSecs;
    
    %-----------------------------------------------------------------
    % baseline fixation - currently 10 seconds = 4*Volumes (2.5 TR)
    while GetSecs < prebaseline + baseline_fixation
        %    Screen(w,'Flip', anchor);
        CenterText(w,'+', white,0,0);
        Screen('TextSize',w, 60);
        Screen(w,'Flip');
        
    end
    
    if use_eyetracker
        %   Eyelink MSG
        % ---------------------------
        % messages to save on each trial ( trial number, onset and RT)
        Eyelink('Message', Eventflag(GenFlags.Fixation.str,task,runNum,1,anchor)); % mark start time in file
    end
    %-----------------------------------------------------------------
    
    [shuff_names{runNum},shuff_ind{runNum}] = Shuffle(vars{1});
    trialType{runNum} = vars{2};
    bidIndex{runNum} = vars{4};
    shuff_bidIndex{runNum} = bidIndex{runNum}(shuff_ind{runNum});
    itemnameIndex{runNum} = vars{5};
    shuff_itemnameIndex{runNum} = itemnameIndex{runNum}(shuff_ind{runNum});
    
    load('Misc/demo_go.mat');
    
    
    Audio_time{runNum}(1:length(trialType{1}),1) = 999;
    respTime{runNum}(1:length(trialType{1}),1) = 999;
    keyPressed{runNum}(1:length(trialType{1}),1) = 999;
    
    for i = 1:length(demoNames)
        demofood_items{i} = imread(sprintf('stim/demo/%s',demoNames{i}));
    end
    
    
    
    if runNum>1
        
        Ladder1{runNum}(1,1) = Ladder1end{runNum-1};
        Ladder2{runNum}(1,1) = Ladder2end{runNum-1};
    end
    
    
    
    for trial = 1:length(demoNames) %length(stop{runNum}),   % To cover all the items in one run.
        
        Screen('PutImage',w,demofood_items{trial});
        Screen('Flip',w,onsets(trial)+runStartTime); % display images according to Onset times
        image_start_time = GetSecs;
        actual_onset_time{runNum}(trial,1) = image_start_time - anchor;
        
        if use_eyetracker
            Eyelink('Message', Eventflag(GenFlags.TrialStart.str,task,runNum,trial,anchor)); % mark start time in file
        end
        
        %-----------------------------------------------------------------
        % get response for all trial types
        noresp = 1;
        notone = 1;
        KbQueueStart;
        
        while (GetSecs-image_start_time < image_duration)
            
            if  demo_go(trial) == 11 && (GetSecs - image_start_time >= Ladder1{runNum}(length(Ladder1{runNum}),1)/1000) && notone %demo_go contains the information if a certain image is a demo_go trial or not
                % Beep!
                %                 PsychPortAudio('FillBuffer', pahandle, wave);
                PsychPortAudio('Start', pahandle, 1, 0, 0);
                %                 play(Audio);
                notone = 0;
                Audio_time{runNum}(trial,1) = GetSecs - image_start_time;
                
                %   Eyelink MSG
                % ---------------------------
                if use_eyetracker
                    Eyelink('Message',Eventflag(GenFlags.CueStart.str,task,runNum,trial,anchor));
                end
                
                % look for response
                [pressed, firstPress, ~, ~, ~] = KbQueueCheck;
                if pressed && noresp
                    findfirstPress = find(firstPress);
                    respTime{runNum}(trial,1) = firstPress(findfirstPress(1))-image_start_time;
                    tmp = KbName(firstPress);
                    if ischar(tmp) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed{runNum} to be a char, so this converts it and takes the first key pressed
                        tmp = char(tmp);
                    end
                    keyPressed{runNum}(trial,1) = tmp(1);
                    
                    %   Eyelink MSG
                    % ---------------------------
                    if use_eyetracker
                        Eyelink('Message',Eventflag(GenFlags.Response.str,task,runNum,trial,anchor));
                    end
                    
                    % different response types in scanner or in testing room
                    if isMRI == 1
                        if keyPressed{runNum}(trial,1) == blue || keyPressed{runNum}(trial,1) == yellow
                            noresp = 0;
                            
                            if respTime{runNum}(trial,1) < Ladder1{runNum}(length(Ladder1{runNum}),1)/1000
                                respInTime{runNum}(trial,1) = 11; %was a GO trial with HV item but responded before SS
                            else
                                respInTime{runNum}(trial,1) = 110; %was a Go trial with HV item but responded after SS within 1000 msec
                            end
                        end
                    else
                        if keyPressed{runNum}(trial,1) == BUTTON %|| keyPressed{runNum}(trial,1) == RIGHT
                            noresp = 0;
                            
                            if respTime{runNum}(trial,1) < Ladder1{runNum}(length(Ladder1{runNum}),1)/1000
                                respInTime{runNum}(trial,1) = 11; %was a GO trial with HV item but responded before SS
                            else
                                respInTime{runNum}(trial,1) = 110; %was a Go trial with HV item and responded after SS within 1000 msec - good trial
                            end
                        end
                    end
                end
                
                %Low-Valued BEEP items
                %---------------------------
            elseif  demo_go(trial) == 22 && (GetSecs - image_start_time >= Ladder2{runNum}(length(Ladder2{runNum}),1)/1000) && notone % demo_go contains the information if a certain image is a GO trial or not
                
                % Beep!
                %                 PsychPortAudio('FillBuffer', pahandle, wave);
                PsychPortAudio('Start', pahandle, 1, 0, 0);
                %                 play(Audio);
                notone = 0;
                Audio_time{runNum}(trial,1) = GetSecs-image_start_time;
                %   Eyelink MSG
                % ---------------------------
                if use_eyetracker
                    Eyelink('Message',Eventflag(GenFlags.CueStart.str,task,runNum,trial,anchor));
                end
                
                %   look for response
                [pressed, firstPress, ~, ~, ~] = KbQueueCheck;
                if pressed && noresp
                    findfirstPress = find(firstPress);
                    respTime{runNum}(trial,1) = firstPress(findfirstPress(1))-image_start_time;
                    
                    tmp = KbName(firstPress);
                    if ischar(tmp) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed{runNum} to be a char, so this converts it and takes the first key pressed
                        tmp = char(tmp);
                    end
                    keyPressed{runNum}(trial,1) = tmp(1);
                    %   Eyelink MSG
                    % ---------------------------
                    if use_eyetracker
                        Eyelink('Message',Eventflag(GenFlags.Response.str,task,runNum,trial,anchor));
                    end
                    
                    %   different response types in scanner or in testing room
                    if isMRI == 1
                        if keyPressed{runNum}(trial,1) == blue || keyPressed{runNum}(trial,1) == yellow
                            noresp = 0;
                            if respTime{runNum}(trial,1) < Ladder2{runNum}(length(Ladder2{runNum}),1)/1000
                                respInTime{runNum}(trial,1) = 22; %was a GO trial with LV item but responded before SS
                            else
                                respInTime{runNum}(trial,1) = 220; %was a Go trial with LV item but responded after SS within 1000 msec
                            end
                        end
                    else
                        if keyPressed{runNum}(trial,1) == BUTTON %| keyPressed{runNum}(trial,1)==RIGHT
                            noresp = 0;
                            if respTime{runNum}(trial,1) < Ladder2{runNum}(length(Ladder2{runNum}),1)/1000
                                respInTime{runNum}(trial,1) = 22;  %was a GO trial with LV item but responded before SS
                            else
                                respInTime{runNum}(trial,1) = 220; %was a Go trial with LV item and responded after SS within 1000 msec - good trial
                            end
                        end
                    end
                end % end if pressed && noresp
                
                %No-BEEP
                %---------------------------
            elseif mod(demo_go(trial),11) ~= 0 && noresp % these will now be the NOGO trials
                
                %   look for response
                [pressed, firstPress, ~, ~, ~] = KbQueueCheck;
                if pressed && noresp
                    findfirstPress = find(firstPress);
                    respTime{runNum}(trial,1) = firstPress(findfirstPress(1))-image_start_time;
                    tmp = KbName(firstPress);
                    if ischar(tmp) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed{runNum} to be a char, so this converts it and takes the first key pressed
                        tmp = char(tmp);
                    end
                    keyPressed{runNum}(trial,1) = tmp(1);
                    %   Eyelink MSG
                    % ---------------------------
                    if use_eyetracker
                        Eyelink('Message',Eventflag(GenFlags.Response.str,task,runNum,trial,anchor));
                    end
                    
                    % different response types in scanner or in testing room
                    if isMRI == 1
                        if keyPressed{runNum}(trial,1) == blue || keyPressed{runNum}(trial,1) == yellow
                            noresp = 0;
                            if demo_go(trial) == 12
                                respInTime{runNum}(trial,1) = 12; % a NOGO trial but responded within 1000 msec HV item - not good but don't do anything
                            else
                                respInTime{runNum}(trial,1) = 24; % a NOGO trial but responded within 1000 msec LV item - not good but don't do anything
                            end
                        end
                    else
                        if keyPressed{runNum}(trial,1) == BUTTON %| keyPressed{runNum}(trial,1)==RIGHT
                            noresp = 0;
                            if demo_go(trial) == 12
                                respInTime{runNum}(trial,1) = 12; %% a NOGO trial but responded within 1000 msec HV item - not good but don't do anything
                            else
                                respInTime{runNum}(trial,1) = 24; %% a NOGO trial but responded within 1000 msec LV item - not good but don't do anything
                            end
                        end
                    end % end if isMRI == 1
                end % end if pressed && noresp
            end % end if  demo_go(trial) == 11 && (GetSecs - image_start_time >= Ladder1{runNum}(length(Ladder1{runNum}),1)/1000) && notone
            
        end % End big while waiting for response within 1000 msec
        
        
        %   Close the Audio port
        %---------------------------
        %         PsychPortAudio('Stop', pahandle); % Close the Audio part
        
        %   Show fixation
        %---------------------------
        CenterText(w,'+', white,0,0);
        Screen('TextSize',w, 60);
        Screen(w,'Flip', image_start_time+1);
        fix_time{runNum}(trial,1) = GetSecs ;
        fixcrosstime{runNum} = GetSecs;
        
        if use_eyetracker
            Eyelink('Message',Eventflag(GenFlags.Fixation.str,task,runNum,trial,anchor));
        end
        
        
        if noresp == 1
            %---------------------------
            % these are additional 500msec to monitor responses
            while (GetSecs-fix_time{runNum}(trial,1) < 0.5)
                
                [pressed, firstPress, ~, ~, ~] = KbQueueCheck;
                if pressed && noresp
                    findfirstPress = find(firstPress);
                    respTime{runNum}(trial,1) = firstPress(findfirstPress(1))-image_start_time;
                    
                    tmp = KbName(firstPress);
                    if ischar(tmp) == 0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed{runNum} to be a char, so this converts it and takes the first key pressed
                        tmp = char(tmp);
                    end
                    keyPressed{runNum}(trial,1) = tmp(1);
                    %   Eyelink MSG
                    % ---------------------------
                    if use_eyetracker
                        Eyelink('Message',Eventflag(GenFlags.Response.str,task,runNum,trial,anchor));
                    end
                    if isMRI == 1
                        if keyPressed{runNum}(trial,1) == blue || keyPressed{runNum}(trial,1) == yellow
                            noresp = 0;
                            switch demo_go(trial)
                                case 11
                                    if respTime{runNum}(trial,1) >= 1
                                        respInTime{runNum}(trial,1) = 1100; % a Go trial and  responded after 1000msec  HV item - make it easier decrease SSD
                                    elseif respTime{runNum}(trial,1) < 1
                                        respInTime{runNum}(trial,1) = 110;
                                    end
                                case 22
                                    if respTime{runNum}(trial,1) >= 1
                                        respInTime{runNum}(trial,1) = 2200; % a Go trial and  responded after 1000msec  HV item - make it easier decrease SSD
                                    elseif respTime{runNum}(trial,1) < 1
                                        respInTime{runNum}(trial,1) = 220;
                                    end
                                case 12
                                    respInTime{runNum}(trial,1) = 12; % a NOGO trial and responded after 1000 msec  HV item - don't touch
                                case 24
                                    respInTime{runNum}(trial,1) = 24; % % a NOGO trial and  responded after 1000 msec HV item - don't touch
                            end
                        end
                        
                    else
                        if keyPressed{runNum}(trial,1) == BUTTON % | keyPressed{runNum}(trial,1)==RIGHT
                            noresp = 0;
                            switch demo_go(trial)
                                case 11
                                    
                                    if respTime{runNum}(trial,1) >= 1
                                        respInTime{runNum}(trial,1) = 1100; % a Go trial and  responded after 1000msec  HV item - make it easier decrease SSD
                                    elseif respTime{runNum}(trial,1) < 1
                                        respInTime{runNum}(trial,1) = 110;
                                    end
                                    
                                case 22
                                    
                                    if respTime{runNum}(trial,1) >= 1
                                        respInTime{runNum}(trial,1) = 2200; % a Go trial and  responded after 1000msec  HV item - make it easier decrease SSD
                                    elseif respTime{runNum}(trial,1) < 1
                                        respInTime{runNum}(trial,1) = 220;
                                    end
                                    
                                case 12
                                    respInTime{runNum}(trial,1) = 12;% a NOGO trial and didnt respond on time HV item - don't touch
                                case 24
                                    respInTime{runNum}(trial,1) = 24;% a NOGO trial and didnt respond on time LV item - don't touch
                                    
                            end
                        end
                    end
                end
            end % End while of additional 500 msec
        else % the subject has already responded during the first 1000 ms
            WaitSecs(0.5);
        end % end if noresp
        
        
        %%%%% This is where its all decided !
        if noresp
            switch demo_go(trial)
                case 11
                    respInTime{runNum}(trial,1) = 1; %unsuccessful Go trial HV - didn't press a button at all - trial too hard - need to decrease ladder
                case 22
                    respInTime{runNum}(trial,1) = 2; % unsuccessful Go trial LV - didn't press a button at all - trial too hard - need to decrease ladder
                case 12
                    respInTime{runNum}(trial,1) = 120; % ok NOGO trial didn't respond after 1500 msec in NOGO trial HV
                case 24
                    respInTime{runNum}(trial,1) = 240; % ok NOGO trial didn't respond after 1500 msec in NOGO trial LV
            end
        end
        
        
        switch respInTime{runNum}(trial,1)
            case 1 % didn't respond even after 1500 msec on go trial - make it easier decrease SSD by step
                if (Ladder1{runNum}(length(Ladder1{runNum}),1) < 0.001)
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1);
                else
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1)-Step;
                end;
                
            case 2 % didn't respond even after 1500 msec on go trial - make it easier decrease SSD by step
                if (Ladder2{runNum}(length(Ladder2{runNum}),1) < 0.001)
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1);
                else
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1)-Step;
                end;
                
                
            case 1100
                if respTime{runNum}(trial,1)>1%  responded after 1500 msec on go trial - make it easier decrease SSD by step
                    if (Ladder1{runNum}(length(Ladder1{runNum}),1) < 0.01)
                        Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1);
                    else
                        Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1)-Step;
                    end;
                end
                
            case 2200 %  responded after 1500 msec on go trial - make it easier decrease SSD by step
                if respTime{runNum}(trial,1)>1
                    if (Ladder2{runNum}(length(Ladder2{runNum}),1) < 0.01)
                        Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1);
                    else
                        Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1)-Step;
                    end;
                end
                
                
            case 11
                if (Ladder1{runNum}(length(Ladder1{runNum}),1) > 999.999); %was a GO trial with HV item but responded before SS make it harder - increase SSD by Step/3
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1);
                else
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1)+Step/3;
                end;
                
            case 22
                if (Ladder2{runNum}(length(Ladder2{runNum}),1) > 999.999); %was a GO trial with LV item but responded before SS make it harder - - increase SSD by Step/3
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1);
                else
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1)+Step/3;
                end;
                
            case 110 % pressed after Go signal but below 1000 HV item - - increase SSD by Step/3 - these are the good trials!
                if (Ladder1{runNum}(length(Ladder1{runNum}),1) > 999.999);
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1);
                else
                    Ladder1{runNum}(length(Ladder1{runNum})+1,1) = Ladder1{runNum}(length(Ladder1{runNum}),1)+Step/3;
                end;
                
            case 220 % pressed after Go signal but below 1000 LV item - - increase SSD by Step/3 - these are the good trials!
                if (Ladder2{runNum}(length(Ladder2{runNum}),1) > 999.999);
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1);
                else
                    Ladder2{runNum}(length(Ladder2{runNum})+1,1) = Ladder2{runNum}(length(Ladder2{runNum}),1)+Step/3;
                end;
                
        end % end switch respInTime{runNum}(pos,1)
        
        KbQueueFlush;
        
        % save data to .txt file
        
        fprintf(fid1,'%s\t %d\t %s\t %d\t %d\t %d\t %d\t %d\t %d\t %.2f\t %d\t %d\t \n', subjectID, runNum, demoNames{trial}, actual_onset_time{runNum}(trial,1), demo_go(trial), respTime{runNum}(trial,1)*1000, respInTime{runNum}(trial,1), Audio_time{runNum}(trial,1)*1000, keyPressed{runNum}(trial,1),   fix_time{runNum}(trial,1)-anchor, Ladder1{runNum}(length(Ladder1{runNum})), Ladder2{runNum}(length(Ladder2{runNum})));
        
        
    end; % %% End the big trial loop showing all the images in one run.
    
    Ladder1end{runNum} = Ladder1{runNum}(length(Ladder1{runNum}));
    Ladder2end{runNum} = Ladder2{runNum}(length(Ladder2{runNum}));
    correct{runNum}(1) = 0;
    
    
    correct{runNum}(1) = length(find(respInTime{runNum} == 110 | respInTime{runNum} == 220 | respInTime{runNum} == 1100 | respInTime{runNum} == 2200 | respInTime{runNum} == 11 | respInTime{runNum} == 22 ) );
    
    
    
end % end for runNum=1:1

fclose(fid1);

%   Eyelink MSG
% ---------------------------
if use_eyetracker
    Eyelink('Message',Eventflag(GenFlags.RunEnd.str,task,runNum,trial,anchor));
end

% mean_RT{runNum} = mean(respTime{runNum}(find(respInTime{runNum} == 110 | respInTime{runNum} == 220 | respInTime{runNum} == 1100 | respInTime{runNum} == 2200 ) ));

Screen('TextSize', w, 30); %Set textsize
CenterText(w,strcat(sprintf('You responded on %.2f', ((correct{runNum}(1))/2)*100), '% of Go trials'), white, 0,-270);

Screen('Flip',w);



% save variables to a file, with date\time stamp
outfile = strcat(outputPath, sprintf('/%s_training_demo_%s.mat', subjectID, timestamp));

% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;
clear food_items demofood_items Instructions*;

save(outfile);

%%   Finishing eye tracking  %
    %---------------------------------------------------------------
    if use_eyetracker
        
        %---------------------------
        % finish up: stop recording eye-movements,
        % close graphics window, close data file and shut down tracker
        Eyelink('StopRecording');
        %   Eyelink MSG
        % ---------------------------
        Eyelink('Message',['Eyetracking_closeTime: ',num2str(GetSecs-anchor)]);
        WaitSecs(.1);
        Eyelink('CloseFile');
        
        
        % download data file
        try
            fprintf('Receiving data file ''%s''\n', edfFile );
            status=Eyelink('ReceiveFile');
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2==exist(edfFile, 'file')
                fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
            end
        catch rdf
            fprintf('Problem receiving data file ''%s''\n', edfFile );
            rdf;
        end
        
        
        if dummymode==0
            movefile(edfFile,['./Output/', subjectID,'_',task,'_eyetracking_', timestamp,'.edf']);
        end;
    end

% Close the audio device:
% PsychPortAudio('Close', pahandle);
%rethrow(lasterror);
Screen('TextSize',w,40);
Screen('TextFont',w,'Ariel');
WaitSecs(2);

Screen('TextSize',w, 36);
CenterText(w,'The Demo is done. Questions?', Green,0,-170);
CenterText(w,sprintf('Press any key to continue') ,Green,0,0);
Screen('Flip',w);

noresp = 1;
while noresp,
    [keyIsDown,~,~] = KbCheck; %(experimenter_device);
    if keyIsDown && noresp,
        noresp = 0;
    end;
end;
WaitSecs(0.001);
KbQueueFlush;


Screen('CloseAll');
ShowCursor;


end % end function

