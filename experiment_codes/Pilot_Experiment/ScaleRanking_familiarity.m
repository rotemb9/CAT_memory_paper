function ScaleRanking_familiarity(subjectID,use_eyetracker,skip_synctest, debugging_mode)
%=========================================================================
% Scale Raniking Code - This code will present a set of stimuli (randomly
% ordered. It reads all images in stim/ScaleRanking/ directory. Eacth type
% of stimuli should be arranges in independent directory according to the
% following coding 'Stim_XX':
% 1 - Snacks, 2 - Faces, 3 - Fractals, 4 - PosIAPS, 5 - NegIAPS, 6 - FamiliarFaces
%
% Coded by: Tom Salomon and Gilad Levy.
% Based on a BDM code by Rani Gera
% July, 2017
%=========================================================================

if nargin < 4
    debugging_mode = 0;
end

if nargin<3
    skip_synctest=0; % set to 1/0 to turn on/off Screen sync test
    if nargin <2
        use_eyetracker=0; % set to 1/0 to turn on/off eyetracker functions
    end
end

% Don't use eyetracker
use_eyetracker=0;

c=clock;
hr=num2str(c(4));
minute=num2str(c(5));
timestamp=[date,'_',hr,'h',minute,'m'];

rng shuffle

tic

%--------------------------------------------------------------------------
%% PARAMETERS FOR THE RANKING AXIS
%--------------------------------------------------------------------------
% Parameters for RANKING RANGE
RankingMin = 0;
RankingMax = 10;

% Parameters for Creatining the RANKING AXIS:
RelativeSizeOfRankingAxisFromScreenSize = 1/3;
YaxisRelativeMiddlePoint = 0.9;
YaxisWidthFactor = 0.0102;

% Parameters for the MOVING INDICATOR on the ranking axis:
penWidth = 3;
AdditionalMovingIndicatorLengthFactor = 0.018; % Extension from Each side of the ranking axis.

% Parameters for FIGURES PRESENTATION:
TextSizeForFiguresOnAxis = 30;
DistanceOfFiguresFactor = 0.0065;

% Fixation cross fix:
FixForFixationCrossLocation = -33.5; % A fix for fixation cross on center text to be in center on the Y axis. Relevant for text size 60.


SamplingRate = 0.001; %Sampling Rate of Mouse Position
MaxTimeRank = 4;
TrialTime = 0;

%---------------------------------------------------------------

%% Load chosen language:
%---------------------------------------------------------------
%TextBlock = Language(ChosenLanguage, 'BDM_Fractals');% mfilename gets the name of the running function;

%---------------------------------------------------------------
%% 'INITIALIZE Screen variables'
%---------------------------------------------------------------

if skip_synctest==1
    Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference', 'VisualDebuglevel', 0); %No PTB intro screen
    
else
    Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
end
screennum = min(Screen('Screens'));

pixelSize=32;
% [w] = Screen('OpenWindow',screennum,[],[0 0 800 800],pixelSize);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

% Here Be Colors
black=BlackIndex(w); % Should equal 0.
white=WhiteIndex(w); % Should equal 255.
green=[0 255 0];


% set up Screen positions for stimuli
[wWidth, wHeight]=Screen('WindowSize', w);
xcenter=wWidth/2;
ycenter=wHeight/2;

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);

% text stuffs
theFont='Arial';
Screen('TextFont',w,theFont);

instrSZ=40;
betsz=60;

%--------------------------------------------------------------------------
%% SETTINGS FOR THE RANKING AXIS
%--------------------------------------------------------------------------
% Settings for Creatining the RANKING AXIS:
AxisFromX = wWidth*RelativeSizeOfRankingAxisFromScreenSize; % Starting point on X axis
AxisToX = wWidth*(1-RelativeSizeOfRankingAxisFromScreenSize); % Ending point on X axis
AxisFromY = round(wHeight * (YaxisRelativeMiddlePoint - (YaxisRelativeMiddlePoint * YaxisWidthFactor))); % Starting point on Y axis
AxisToY = round(wHeight * (YaxisRelativeMiddlePoint + (YaxisRelativeMiddlePoint * YaxisWidthFactor))); % Ending point on Y axis

% Settings for the MOVING INDICATOR on the ranking axis:
CenterOfMovingIndicator = mean([AxisFromY AxisToY]);
AdditionToYAxisFromEachSide = wHeight*AdditionalMovingIndicatorLengthFactor;

% Settings for FIGURES PRESENTATION:
RankingIntegers = RankingMax - RankingMin + 1;
SpotsForIndicatorsOnAxis = linspace(AxisFromX, AxisToX, RankingIntegers);
DistanceOfFiguresFromAxis = round(wHeight * DistanceOfFiguresFactor) + AdditionToYAxisFromEachSide ;
FixForFiguresOnXaxis = round(wHeight * 0.05);
FixForFiguresOnXaxis = round(4*AdditionToYAxisFromEachSide);
% FixForFiguresOnXaxis=100 ;

ScaleText=cell(1,length(SpotsForIndicatorsOnAxis));
ScaleText{1}='not familiar';
ScaleText{end}='very familiar';
%--------------------------------------------------------------------------
%% Locations, File Types & Names PARAMETERS
%--------------------------------------------------------------------------
%Stimuli_index - 0 for politicians, 1 for neutral faces.If Stimuli_Index is not provided - default folder is ./stim/BDM/.

% Stimuli:
stimuli_dirs=dir( './stim/ScaleRanking/Stim_6*');
StimTypeInd_list=[];
StimType_list=cell(0);
imageArrays=cell(0);
stimuli_images=[];
for stim_ind = 1:length(stimuli_dirs)
    StimLocation = ['./stim/ScaleRanking/',stimuli_dirs(stim_ind).name,'/'];
    
    
    StimTypeInd=str2double(stimuli_dirs(stim_ind).name(6:end)); % 1 - Snacks, 2 - Faces, 3 - Fractals, 4 - PosIAPS, 5 - NegIAPS, 6 - FamiliarFaces
    switch StimTypeInd
        case 1
            StimType = 'Snacks';
            StimFileType = 'bmp';
        case 2
            StimType = 'Faces';
            StimFileType = 'jpg';
        case 3
            StimType = 'Fractals';
            StimFileType = 'jpg';
        case 4
            StimType = 'PosIAPS';
            StimFileType = 'jpg';
        case 5
            StimType = 'NegIAPS';
            StimFileType = 'jpg';
        case 6
            StimType = 'FamiliarFaces';
            StimFileType = 'jpg';
    end
    
    % 'LOAD image arrays'
    %---------------------------------------------------------------
    images2load=dir([StimLocation '*.' StimFileType]);
    
    
    StimTypeInd_list(end+1:end+length(images2load))=StimTypeInd;
    for i=1:length(images2load)
        imageArrays{end+1}=imread([StimLocation,images2load(i).name]);
        StimType_list{end+1}=StimType;
    end
    stimuli_images=[stimuli_images;images2load];
end

% Shuffle stimuli
shuffledlist=Shuffle(1:length(imageArrays));
imageArrays=imageArrays(shuffledlist);
StimTypeInd_list=StimTypeInd_list(shuffledlist);
StimType_list=StimType_list(shuffledlist);
stimuli_images=stimuli_images(shuffledlist);

% Set Output folder
OutputFolder = 'Output/';

MousePositionX = cell(length(imageArrays),1);
MousePositionY = cell(length(imageArrays),1);
%---------------------------------------------------------------
%% 'Write output file header'
%---------------------------------------------------------------
fid1=fopen([OutputFolder subjectID '_Scale_Familiarity_Ranking_' timestamp '.txt'], 'a');
fprintf(fid1,'subjectID\truntrial\tonsettime\tName\tRanking\tRT\tStimType\tStimTypeInd\n'); %write the header line

%---------------------------------------------------------------
%% Initializing eye tracking system %
%-----------------------------------------------------------------
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
    
    [v vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,HREF,AREA');
    
    % open file to record data to
    edfFile='Rank.edf';
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
%==============================================
%% 'Display Main Instructions'
%---------------------------------------------------------------

Screen('TextSize',w, round(instrSZ*1.5));
CenterText(w,'Ranking Familiarity', green,0,-260);

Screen('TextSize',w, instrSZ);
CenterText(w,'You are about to see a set of images', white,0,-100);
CenterText(w,'For each image, please indicate how much the person in the image is familiar to you', white,0,-20);
CenterText(w,'using the following scale. Use the mouse cursor', white,0,60);
CenterText(w,'Press any key to start', green,0,220);
HideCursor;
Screen('Flip', w);
WaitSecs(0.05); % prevent key spillover
noresp=1;
while noresp
    [keyIsDown] = KbCheck;
    if keyIsDown && noresp
        noresp = 0;
    end
end

Screen('TextSize',w, betsz);
CenterText(w,'+', white,0,FixForFixationCrossLocation);
Screen(w,'Flip');
WaitSecs(0.3);

runStart=GetSecs;

if use_eyetracker
    % start recording eye position
    %---------------------------
    Eyelink('StartRecording');
    %   Eyelink MSG
    % ---------------------------
    WaitSecs(.05);
end

%---------------------------------------------------------------
%% 'Run Trials'
%---------------------------------------------------------------

% for trial=1:5 % debugging
num_trials = length(stimuli_images);

if debugging_mode
    num_trials = 4;
end

for trial=1:num_trials
    
    TrialStartTime = GetSecs;
    
    bid=[];
    noresp=1;
    Screen('TextSize',w,TextSizeForFiguresOnAxis);
    eventTime=[];
    HideCursor;
    SetMouse(xcenter,ycenter);
    if isempty(eventTime) % recording the presentation start time
        eventTime = GetSecs-runStart;
    end
    
    if use_eyetracker
        %   Eyelink MSG
        % ---------------------------
        % messages to save on each trial ( trial number, onset and RT)
        Eyelink('Message', Eventflag(GenFlags.TrialStart.str,GenFlags.BDM.str,1,trial,runStart)); % mark start time in file
    end
    
    while (TrialStartTime + TrialTime) > GetSecs
        Screen('PutImage',w,imageArrays{trial});
        Screen(w,'Flip');
    end
    % Show scale
    if use_eyetracker
        %   Eyelink MSG
        % ---------------------------
        % messages to save on each trial ( trial number, onset and RT)
        Eyelink('Message', Eventflag(GenFlags.ScaleStart.str,GenFlags.BDM.str,1,trial,runStart)); % mark start time in file
    end
    ScaleAppearTime = GetSecs;
    SetMouse(xcenter,ycenter);
    MouseTrajectoryTime = 0;
    
    while (ScaleAppearTime + MaxTimeRank) > GetSecs
        Screen('PutImage',w,imageArrays{trial});
        
        if noresp==1
            Screen('FillRect', w ,[211 211 211] ,[AxisFromX, AxisFromY,  AxisToX, AxisToY]);
            for i = 1:length(SpotsForIndicatorsOnAxis)
                DrawFormattedText(w, ScaleText{i}, SpotsForIndicatorsOnAxis(i)- FixForFiguresOnXaxis, CenterOfMovingIndicator+DistanceOfFiguresFromAxis, [255 255 255]);
            end
            Screen(w,'Flip');
            ShowCursor;
            % Track cursor movement and check for response
            [CurrentX,CurrentY,buttons] = GetMouse(w);
            if  GetSecs > (MouseTrajectoryTime + SamplingRate) %% Get mouse trajectory.
                MousePositionX{trial}(end+1,:) = CurrentX;
                MousePositionY{trial}(end+1,:) = CurrentY;
                MouseTrajectoryTime = GetSecs;
            end
            
            if CurrentX >= AxisFromX && CurrentX <= AxisToX && CurrentY >= AxisFromY - AdditionToYAxisFromEachSide && CurrentY <= AxisToY + AdditionToYAxisFromEachSide
                Screen('PutImage',w,imageArrays{trial});
                Screen('FillRect', w ,[211 211 211] ,[AxisFromX, AxisFromY,  AxisToX, AxisToY]);
                for i = 1:length(SpotsForIndicatorsOnAxis)
                    DrawFormattedText(w, ScaleText{i}, SpotsForIndicatorsOnAxis(i)- FixForFiguresOnXaxis, CenterOfMovingIndicator+DistanceOfFiguresFromAxis, [255 255 255]);
                end
                Screen('DrawLine', w ,[0 0 255], CurrentX, CenterOfMovingIndicator+AdditionToYAxisFromEachSide, CurrentX, CenterOfMovingIndicator-AdditionToYAxisFromEachSide ,penWidth);
                Screen(w,'Flip');
                if buttons(1) == 1
                    bid = (CurrentX - AxisFromX) / (AxisToX - AxisFromX) * (RankingMax - RankingMin) + RankingMin; % Number of pixels from X axis beggining / Length of the axis * Units + Beggining of units.
                    noresp = 0;
                    HideCursor;
                    Screen('PutImage',w,imageArrays{trial});
                    Screen(w,'Flip');
                    while any(buttons) % wait for release
                        [~,~,buttons] = GetMouse;
                    end
                    respTime = GetSecs - TrialStartTime;
                    if use_eyetracker
                        %   Eyelink MSG
                        % ---------------------------
                        % messages to save on each trial ( trial number, onset and RT)
                        Eyelink('Message', Eventflag(GenFlags.Response.str,GenFlags.BDM.str,1,trial,runStart));
                    end
                end
            end
        end
    end
    
    if noresp == 1
        bid = 999;
        respTime = 999;
        if use_eyetracker
            %   Eyelink MSG
            % ---------------------------
            Eyelink('Message', Eventflag(GenFlags.RespondFaster.str,GenFlags.BDM.str,1,trial,runStart));
        end
        RespFasterMSGTime = GetSecs;
        while GetSecs < (RespFasterMSGTime + 1.5)
            CenterText(w,sprintf('You must respond faster!') ,white,0,0);
            %       Screen('DrawText', w, 'You must respond faster!', xcenter-450, ycenter, white);
            Screen(w,'Flip');
        end
        noresp = 0;
    end
    
    
    %-----------------------------------------------------------------
    % show fixation ITI
    FixationStartTime = GetSecs;
    if use_eyetracker
        %   Eyelink MSG
        % ---------------------------
        Eyelink('Message', Eventflag(GenFlags.Fixation.str,GenFlags.BDM.str,1,trial,runStart));
    end
    while GetSecs < (FixationStartTime + 0.5)
        Screen('TextSize',w, betsz);
        CenterText(w,'+', white,0,FixForFixationCrossLocation);
        Screen(w,'Flip');
    end
    
    
    
    %-----------------------------------------------------------------
    % write to output file
    
    fprintf(fid1,'%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\n', subjectID, trial, eventTime, stimuli_images(trial).name, bid, respTime,StimType_list{trial},StimTypeInd_list(trial));
end

if use_eyetracker
    %   Eyelink MSG
    % ---------------------------
    Eyelink('Message', Eventflag(GenFlags.RunEnd.str,GenFlags.BDM.str,1,trial,runStart));
end

HideCursor;
CenterText(w,'Thank you!', green,0,20);
Screen('Flip', w, 0,1);

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
        movefile(edfFile,[OutputFolder subjectID '_Scale_Familiarity_Ranking_' timestamp '.edf']);
    end
end

WaitSecs(3); % prevent key spillover

fclose(fid1);
toc

ShowCursor;
Screen('closeall');
