function recognition_confidence_demo_new(isMRI,mainPath,use_eyetracker)

% function recognition_confidence_demo(subjectID,isMRI,mainPath,order, sessionNum,use_eyetracker)
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Rotem Botvinik November 2015 ====================
% edited on March 2017 to include eye-tracking
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%
% This function runs a demo for the recognition session of the cue-approach task.
% Subjects are presented with demo stimuli, and should answer whether they
% recognize each stimuli from the previous sessions or not.
% In this version, subjects also indicate their level of confidence
% (high/low/uncertain).
% Then, for every item, they are immediately asked whether this item was
% paired with a beep during training, while again indicating their level of confidence.
%
%
% Demo stimuli should be located in the folder [mainPath'/Stim/demo']
%
% This version of the function fits the boost version with training only 40
% items!
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   none...
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   none...
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID = 'test1';
% order = 1;
% isMRI = 4;
% sessionNum = 1;
% mainPath = '/Users/schonberglabimac1/Documents/BMI_BS_40';


tic

rng shuffle

if nargin < 3
    use_eyetracker = 1;
end

% %---------------------------------------------------------------
% %% Ask if you want to use Eye Tracker (DIALOG BOX)
% % =========================================================================
% if use_eyetracker==1
%     use_eyetracker = 0;
%     ask_if_want_eyetracker = questdlg ('Do you want to use Eye Tracker?', 'Yes', 'Yes', 'No', 'No');
%     if strcmp(ask_if_want_eyetracker, 'Yes')
%         use_eyetracker = 1;
%     end
% end

%==========================================================
%% 'INITIALIZE Screen variables to be used in each task'
%==========================================================

Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
screennum = max(Screen('Screens'));

pixelSize = 32;
% [w] = Screen('OpenWindow',screennum,[],[0 0 640 480],pixelSize);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

[wWidth, wHeight] = Screen('WindowSize', w);
xcenter = wWidth/2;
ycenter = wHeight/2;

sizeFactor = 0.8;
stimW = 576*sizeFactor;
stimH = 432*sizeFactor;
rect = [xcenter-stimW/2 ycenter-stimH/2 xcenter+stimW/2 ycenter+stimH/2];


% Colors settings
black = BlackIndex(w); % Should equal 0.
white = WhiteIndex(w); % Should equal 255.
% green = [0 255 0];

Screen('FillRect', w, black);
Screen('Flip', w);

% set up screen positions for stimuli
[wWidth, wHeight] = Screen('WindowSize', w);
xcenter = wWidth/2;
ycenter = wHeight/2;


% text settings
theFont = 'Arial';
Screen('TextFont',w,theFont);
Screen('TextSize',w, 40);

HideCursor;

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
    task = GenFlags.MemoryDemo.str;
    edfFile='RecoDemo.edf';
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

% -----------------------------------------------
%% Load Instructions
% -----------------------------------------------

% Load Hebrew instructions image files
Instructions = dir([mainPath '/Instructions/recognition_confidence_demo_new.JPG' ]);
Instructions_image = imread([mainPath '/Instructions/' Instructions(1).name]);

% -----------------------------------------------
%% Load Questions Images
% -----------------------------------------------

% Load Hebrew questions image files
goNoGoQuestion = dir([mainPath '/Instructions/recognition_goNoGo_question.JPG' ]);
newOldQuestion = dir([mainPath '/Instructions/recognition_oldNew_question.JPG' ]);
goNoGoQuestion_image = imread([mainPath '/Instructions/' goNoGoQuestion(1).name]);
newOldQuestion_image = imread([mainPath '/Instructions/' newOldQuestion(1).name]);
% size_goNoGoQuestion_image = size(goNoGoQuestion_image);
% size_goNoGoQuestion_image = size_goNoGoQuestion_image(1:2);
question1_height = size(newOldQuestion_image,1);
question1_width = size(newOldQuestion_image,2);
question2_height = size(goNoGoQuestion_image,1);
question2_width = size(goNoGoQuestion_image,2);

question1_location = [xcenter-question1_width/2 0.3*ycenter-question1_height/2 xcenter+question1_width/2 0.3*ycenter+question1_height/2];
question2_location = [xcenter-question2_width/2 0.3*ycenter-question2_height/2 xcenter+question2_width/2 0.3*ycenter+question2_height/2];

% -----------------------------------------------
%% Load Answers Images
% -----------------------------------------------

% Load Hebrew answers image files
recognition_answers = dir([mainPath '/Instructions/recognition_answers.jpg' ]);
recognition_answers_image = imread([mainPath '/Instructions/' recognition_answers(1).name]);
answer_images = cell(1,5);
for answerInd = 1:5
    answers = dir([mainPath '/Instructions/recognition_' num2str(answerInd) '.jpg' ]);
    answer_images{answerInd} = imread([mainPath '/Instructions/' answers(1).name]);
end

% the dimensions of the answers' images are 1288 X 142
answers_width = 1342;
answers_height = 154;
answers_location = [xcenter-answers_width/2 1.75*ycenter-answers_height/2 xcenter+answers_width/2 1.75*ycenter+answers_height/2];

%---------------------------------------------------------------
%% 'GLOBAL VARIABLES'
%---------------------------------------------------------------
% timeNow = clock;
% hours = sprintf('%02d', timeNow(4));
% minutes = sprintf('%02d', timeNow(5));
% timestamp = [date,'_',hours,'h',minutes,'m'];
% 
% outputPath = [mainPath '/Output'];

%---------------------------------------------------------------
%% 'Assign response keys'
%---------------------------------------------------------------
KbName('UnifyKeyNames');

if isMRI == 1
%     leftresp = 'b';
%     rightresp = 'y';
    %     badresp = 'x';
else
%     leftresp = 'u';
%     rightresp = 'i';
    %     badresp = 'x';
    key1 = '1'; % high confidence yes / beep
    key2 = '2'; % low confidence yes / beep
    key3 = '3'; % uncertain
    key4 = '4'; % low confidence no / no beep
    key5 = '5'; % high confidence no / no beep
end

%---------------------------------------------------------------
%% 'LOAD image arrays'
%---------------------------------------------------------------

demoStimName = dir([mainPath '/Stim/demo/*.bmp']); % Read demo stimuli

% Read old images to a cell array
imgArraysDemo = cell(1,length(demoStimName));
for i = 1:length(demoStimName)
    imgArraysDemo{i} = imread([mainPath '/Stim/demo/' demoStimName(i).name]);
end

%---------------------------------------------------------------
%% 'Display Main Instructions'
%---------------------------------------------------------------

Screen('TextSize',w, 40);

if isMRI == 1
   Screen('PutImage',w,Instructions_image_fmri);    
    
else
   Screen('PutImage',w,Instructions_image);    
end
Screen(w,'Flip');


noresp = 1;
while noresp,
    [keyIsDown] = KbCheck(-1); % deviceNumber=keyboard
    if keyIsDown && noresp,
        noresp = 0;
    end;
end

WaitSecs(0.001);
runStart = GetSecs;

if use_eyetracker
    % start recording eye position
    %---------------------------
    Eyelink('StartRecording');
    WaitSecs(.05);
    Eyelink('Message',Eventflag(GenFlags.RunStart.str,task,1,1,runStart));
end

Screen('TextSize',w, 60);
Screen('DrawText', w, '+', xcenter, ycenter, white);
Screen(w,'Flip');
if use_eyetracker
    Eyelink('Message',Eventflag(GenFlags.Fixation.str,task,1,1,runStart));
end
WaitSecs(1);

KbQueueCreate;

%---------------------------------------------------------------
%% 'Run Trials'
%---------------------------------------------------------------

% runStart = GetSecs;

shuffleNames = randperm(length(demoStimName));

for trial = shuffleNames(1:8)  
    % isOld part
    if use_eyetracker
        Eyelink('Message',Eventflag(GenFlags.TrialStart.str,task,1,trial,runStart));
    end    
    %-----------------------------------------------------------------
    % display image
    % 3 secs per answer
    
    Screen('PutImage',w, imgArraysDemo{trial},rect); % display item
    Screen('TextSize',w, 40);
    Screen('PutImage',w, newOldQuestion_image,question1_location); % display question
    Screen('PutImage',w, recognition_answers_image,answers_location); % display answers
    Screen(w,'Flip');
    StimOnset_isOld = GetSecs;
    if use_eyetracker
        Eyelink('Message',Eventflag(GenFlags.PresentIsOld.str,task,1,trial,runStart));
    end
    
    KbQueueStart;
    %-----------------------------------------------------------------
    % get response
    
    noresp = 1;
    while noresp && (GetSecs - StimOnset_isOld < 3)
        % check for response
        [keyIsDown, firstPress] = KbQueueCheck(-1);
        
        if keyIsDown && noresp
            findfirstPress = find(firstPress);
            % respTime_isOld = firstPress(findfirstPress(1))-StimOnset_isOld;
            tmp = KbName(findfirstPress);
            if ischar(tmp)==0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                tmp = char(tmp);
            end
            response_isOld = tmp(1);
            if use_eyetracker
                Eyelink('Message',Eventflag(GenFlags.ResponseIsOld.str,task,1,trial,runStart));
            end
            if response_isOld==key1||response_isOld==key2||response_isOld==key3||response_isOld==key4||response_isOld==key5 % A valid response is only 1,2,3,4 or 5
            noresp = 0;
            end
            
        end % end if keyIsDown && noresp
        
    end % end while noresp
    
    %-----------------------------------------------------------------
    if noresp
        % display a message for responding faster
        Screen('TextSize',w, 40);
        CenterText(w,sprintf('You must respond faster!') ,white,0,0);
        Screen(w,'Flip');
        if use_eyetracker
            Eyelink('Message',Eventflag(GenFlags.RespondFaster.str,task,1,trial,runStart));
        end
    else
        % redraw text output with the appropriate colorchanges to highlight
        % response
        Screen('PutImage',w,imgArraysDemo{trial}, rect); % display item
        Screen('TextSize',w, 40);
        Screen('PutImage',w, newOldQuestion_image,question1_location); % display question
        Screen('PutImage',w, answer_images{str2double(response_isOld)},answers_location); % display answers
        Screen(w,'Flip');
    end
    
    WaitSecs(0.5);

    %-----------------------------------------------------------------
%     % show fixation ITI
%     Screen('TextSize',w, 60);
%     Screen('DrawText', w, '+', xcenter, ycenter, white);
%     Screen(w,'Flip');
%     WaitSecs(1);

    KbQueueFlush;    
    
    
    % Go\NoGo question
    Screen('PutImage',w, imgArraysDemo{trial},rect); % display item
    Screen('TextSize',w, 40);
    Screen('PutImage',w, goNoGoQuestion_image,question2_location); % display question2
    Screen('PutImage',w, recognition_answers_image,answers_location); % display answers
    Screen(w,'Flip');
    StimOnset_isGo = GetSecs;
    if use_eyetracker
        Eyelink('Message',Eventflag(GenFlags.PresentIsGo.str,task,1,trial,runStart));
    end
    
    KbQueueStart;
    %-----------------------------------------------------------------
    % get response
    
    
    noresp = 1;
    while noresp && (GetSecs - StimOnset_isGo < 3)
        % check for response
        [keyIsDown, firstPress] = KbQueueCheck(-1);
        
        if keyIsDown && noresp
            findfirstPress = find(firstPress);
            % respTime_isGo = firstPress(findfirstPress(1))-StimOnset_isGo;
            tmp = KbName(findfirstPress);
            if ischar(tmp)==0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                tmp = char(tmp);
            end
            response_isGo = tmp(1);
            if use_eyetracker
                Eyelink('Message',Eventflag(GenFlags.ResponseIsGo.str,task,1,trial,runStart));
            end
            if response_isGo==key1||response_isGo==key2||response_isGo==key3||response_isGo==key4||response_isGo==key5 % A valid response is only 1,2,3,4 or 5
                noresp = 0;
            end
            
        end % end if keyIsDown && noresp
        
    end % end while noresp
    
    
    %-----------------------------------------------------------------
    % redraw text output with the appropriate colorchanges to highlight
    % response
    if noresp
        % display a message for responding faster
        Screen('TextSize',w, 40);
        CenterText(w,sprintf('You must respond faster!') ,white,0,0);
        Screen(w,'Flip');
        if use_eyetracker
            Eyelink('Message',Eventflag(GenFlags.RespondFaster.str,task,1,trial,runStart));
        end
    else
        % redraw text output with the appropriate colorchanges to highlight
        % response
        Screen('PutImage',w,imgArraysDemo{trial}, rect); % display item
        Screen('TextSize',w, 40);
        Screen('PutImage',w, goNoGoQuestion_image,question2_location); % display question
        Screen('PutImage',w, answer_images{str2double(response_isGo)},answers_location); % display answers
        Screen(w,'Flip');
    end
    
    WaitSecs(0.5);
    
    %-----------------------------------------------------------------
    % show fixation ITI
    Screen('TextSize',w, 60);
    Screen('DrawText', w, '+', xcenter, ycenter, white);
    Screen(w,'Flip');
    WaitSecs(1);
    
    
    KbQueueFlush;
    
end % end loop for trial = 1:length(food_images);

if use_eyetracker
    Eyelink('Message',Eventflag(GenFlags.RunEnd.str,task,1,trial,runStart));
end

% End of session screen

Screen('TextSize',w, 40);
CenterText(w,'Thank you!', white,0,-50);
CenterText(w,'The demo is done', white,0,50);
Screen(w,'Flip');

% Closing

WaitSecs(3);
toc
ShowCursor;
Screen('CloseAll');

if use_eyetracker
    %---------------------------------------------------------------
    %%   Finishing eye tracking system %
    %---------------------------------------------------------------
    % STEP 7
    %---------------------------
    % finish up: stop recording eye-movements,
    % close graphics window, close data file and shut down tracker
    Eyelink('StopRecording');
    %   Eyelink MSG
    % ---------------------------
    Eyelink('Message', ['Eyetracking close time ',num2str(GetSecs-runStart)]); % mark start time in file
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
    
    
   % if dummymode==0
        %movefile(edfFile,['./Output/', subjectID,'_recognitionDemo_eyetracking_' timestamp,'.edf']);
   % end
end % end closing eye-tracker

end % end function