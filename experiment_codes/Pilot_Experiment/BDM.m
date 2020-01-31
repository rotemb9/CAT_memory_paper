function BDM(subjectID, sessionNum, type, use_eyetracker, debugging_mode)
%=========================================================================
% BDM task
%=========================================================================
% created by Rani, edited by Roni and then Rotem on June 2017
% function BDM(subjectID, sessionNum, type, use_eyetracker)
% creates the output file: [subjectID '_BDM' num2str(sessionNum) '_' timestamp '.txt']

% types: BDM/demo/faces/fractals etc.

c=clock;
hr=num2str(c(4));
minute=num2str(c(5));
timestamp=[date,'_',hr,'h',minute,'m'];
mainPath=pwd;

rng shuffle

tic

if nargin < 5
    debugging_mode = 0;
end

%--------------------------------------------------------------------------
%% Locations, File Types & Names PARAMETERS
%--------------------------------------------------------------------------
% Stimuli:
switch type
    case 'BDM'
        stim_location = 'Stim/';
        StimuliName = 'snacks';
    case 'demo'
        stim_location = 'Stim/demo/';
        StimuliName = 'snacks_demo';
    otherwise
        stim_location = ['Stim/' type '/'];
        StimuliName = type;
end

% StimFileTypes = {'jpg','bmp','png'};

% for output file
OutputFolder = 'Output/';

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

%--------------------------------------------------------------------------
%% PARAMETERS FOR THE STIMULI SIZE AND LOCATION
%--------------------------------------------------------------------------
% StartingPointOnX = 0;
% StartingPointOnY = 0;
% StimuliWidthScaleFactor = 1;
% StimuliHeightScaleFactor = 0.8361;

%---------------------------------------------------------------

%% 'INITIALIZE Screen variables'
%---------------------------------------------------------------
Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
screennum = min(Screen('Screens'));

pixelSize=32;

%[w] = Screen('OpenWindow',screennum,[],[0 0 1000 800],pixelSize);% %debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

% Here Be Colors
black=BlackIndex(w); % Should equal 0.
white=WhiteIndex(w); % Should equal 255.
%yellow=[0 255 0];


% set up Screen positions for stimuli
[screenXpixels, screenYpixels]=Screen('WindowSize', w);
xcenter=screenXpixels/2;
ycenter=screenYpixels/2;

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);

% text stuffs
theFont='Arial';
Screen('TextFont',w,theFont);

instrSZ=45;
betsz=60;

%--------------------------------------------------------------------------
%% SETTINGS FOR THE RANKING AXIS
%--------------------------------------------------------------------------
% Settings for Creatining the RANKING AXIS:
AxisFromX = screenXpixels*RelativeSizeOfRankingAxisFromScreenSize; % Starting point on X axis
AxisToX = screenXpixels*(1-RelativeSizeOfRankingAxisFromScreenSize); % Ending point on X axis
AxisFromY = round(screenYpixels * (YaxisRelativeMiddlePoint - (YaxisRelativeMiddlePoint * YaxisWidthFactor))); % Starting point on Y axis
AxisToY = round(screenYpixels * (YaxisRelativeMiddlePoint + (YaxisRelativeMiddlePoint * YaxisWidthFactor))); % Ending point on Y axis

% Settings for the MOVING INDICATOR on the ranking axis:
CenterOfMovingIndicator = mean([AxisFromY AxisToY]);
AdditionToYAxisFromEachSide = screenYpixels*AdditionalMovingIndicatorLengthFactor;

% Settings for FIGURES PRESENTATION:
RankingIntegers = RankingMax - RankingMin + 1;
SpotsForIndicatorsOnAxis = linspace(AxisFromX, AxisToX, RankingIntegers);
DistanceOfFiguresFromAxis = round(screenYpixels * DistanceOfFiguresFactor) + AdditionToYAxisFromEachSide ;
FixForFiguresOnXaxis = round(screenYpixels * 0.0074);

%--------------------------------------------------------------------------
%% SETTINGS FOR THE STIMULI SIZE AND LOCATION
%--------------------------------------------------------------------------
% PictureSizeOnX = round(screenXpixels * StimuliWidthScaleFactor);
% PictureSizeOnY = round(screenYpixels * StimuliHeightScaleFactor);
% if strcmp(type, 'faces_BDM')
%     PictureLocationVector = [];
% else
%     PictureLocationVector = [StartingPointOnX, StartingPointOnY, StartingPointOnX + PictureSizeOnX, StartingPointOnY + PictureSizeOnY];
% end


%  PictureLocationVector = []; % Activate it to draw the picture in the original size and in the center.

%---------------------------------------------------------------
%% 'LOAD image arrays'
%---------------------------------------------------------------
stimuli_images = dir([stim_location '*.*' ]);
stimuli_images = stimuli_images(4:end); % remove '.' '..' '.DS_Store'
stimuli_images = Shuffle(stimuli_images);
images_names = extractfield(stimuli_images,'name');
images_names = strcat(stim_location,images_names);
imageArrays = cellfun(@imread,images_names,'UniformOutput',0);

%---------------------------------------------------------------
%% 'Write output file header'
%---------------------------------------------------------------
switch type
    case 'BDM'
    Task=GenFlags.BDM.str;
    output_filename = [OutputFolder subjectID '_' type num2str(sessionNum) '_' timestamp '.txt'];
    case 'demo'
    Task=GenFlags.BDMDemo.str;
    output_filename = [OutputFolder subjectID '_BDM_demo_' timestamp '.txt'];
    otherwise
    output_filename = [OutputFolder subjectID '_BDM_' StimuliName '_' timestamp '.txt'];
end

fid1 = fopen(output_filename, 'a');
fprintf(fid1,'subjectID\truntrial\tonsettime\tName\tBid\tRT\tfirst_mouse_movement\n'); %write the header line

%% Initializing eye tracking system %
%-----------------------------------------------------------------
% use_eyetracker=0; % set to 1/0 to turn on/off eyetracker functions
if use_eyetracker
    dummymode=0;
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
    
    [~, vs]=Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    
    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,HREF,AREA');
    
    % open file to record data to
    edfFile=[Task,'.edf'];
    
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

%---------------------------------------------------------------
%% 'Display Main Instructions'
%---------------------------------------------------------------

% Load Hebrew instructions image file
switch type
    case 'BDM'
        Instructions=dir([mainPath '/Instructions/BDM.JPG' ]);
    case 'demo'
        Instructions=dir([mainPath '/Instructions/BDM_Demo.JPG' ]);
    otherwise
        Instructions=dir([mainPath '/Instructions/BDM_' type '.JPG' ]);    
end

Instructions_name = struct2cell(rmfield(Instructions,{'date','bytes','isdir','datenum'}));
Instructions_image = imread([mainPath '/Instructions/' sprintf(Instructions_name{1})]);

KbQueueCreate;
Screen('PutImage',w,Instructions_image);
Screen(w,'Flip');

noresp=1;
while noresp,
    [keyIsDown] = KbCheck(-1); % deviceNumber=keyboard
    if keyIsDown && noresp,
        noresp=0;
    end;
end;


%---------------------------------------------------------------
%% 'Run Trials'
%---------------------------------------------------------------
runStart = GetSecs;
KbQueueFlush;
ShowCursor;
KbQueueStart;

if use_eyetracker
    % start recording eye position
    %---------------------------
    Eyelink('StartRecording');
    %   Eyelink MSG
    % ---------------------------
    % messages to save on each trial ( trial number, onset and RT)
    Eyelink('Message', Eventflag(GenFlags.RunStart.str,Task,1,1,runStart)); % mark start time in file
    
end

if strcmp(type, 'demo')
    NumTrials = 4 ;% shorter version for demo
else
    NumTrials = length(stimuli_images);
end

if debugging_mode
   NumTrials = 4; 
end

for trial = 1:NumTrials
    %bid = [];
    noresp = 1;
    Screen('TextSize',w,TextSizeForFiguresOnAxis);
    eventTime=[];
    ShowCursor;
    SetMouse(xcenter,ycenter);
    if use_eyetracker
        Eyelink('Message',Eventflag(GenFlags.TrialStart.str,Task,1,trial,runStart));
    end
    % Intialize variables for checking when the mouse first moved:
    FirstlyMoved = 0;
    while noresp
        %[keyIsDown, firstPress] = KbQueueCheck;
        
        % Track cursor movement and check for response
        [CurrentX,CurrentY,buttons] = GetMouse(w);
        if FirstlyMoved == 0 && CurrentX ~= screenXpixels/2 && CurrentY ~= screenYpixels/2
            FirstMouseMovement = GetSecs - runStart - eventTime;
            FirstlyMoved = 1;
        end
        if CurrentX >= AxisFromX && CurrentX <= AxisToX && CurrentY >= AxisFromY - AdditionToYAxisFromEachSide && CurrentY <= AxisToY + AdditionToYAxisFromEachSide
            Screen('PutImage',w,imageArrays{trial});
            Screen('FillRect', w ,[211 211 211] ,[AxisFromX, AxisFromY,  AxisToX, AxisToY]);
            for i = 1:length(SpotsForIndicatorsOnAxis)
                DrawFormattedText(w, num2str(i-1 + RankingMin), SpotsForIndicatorsOnAxis(i)-FixForFiguresOnXaxis, CenterOfMovingIndicator+DistanceOfFiguresFromAxis, [255 255 255]);
            end
            Screen('DrawLine', w ,[0 0 255], CurrentX, CenterOfMovingIndicator+AdditionToYAxisFromEachSide, CurrentX, CenterOfMovingIndicator-AdditionToYAxisFromEachSide ,penWidth);
            Screen(w,'Flip');
            if buttons(1) == 1
                bid = (CurrentX - AxisFromX) / (AxisToX - AxisFromX) * (RankingMax - RankingMin) + RankingMin; % Number of pixels from X axis beggining / Length of the axis * Units + Beggining of units.
                respTime = GetSecs - runStart - eventTime;
                noresp = 0;
                while any(buttons) % wait for release
                    [~,~,buttons] = GetMouse;
                end
                if use_eyetracker
                    Eyelink('Message',Eventflag(GenFlags.Response.str,Task,1,trial,runStart));
                end
            end
        else
            Screen('PutImage',w,imageArrays{trial});
            Screen('FillRect', w ,[211 211 211] ,[AxisFromX, AxisFromY,  AxisToX, AxisToY]);
            for i = 1:length(SpotsForIndicatorsOnAxis)
                DrawFormattedText(w, num2str(i-1 + RankingMin), SpotsForIndicatorsOnAxis(i)- FixForFiguresOnXaxis, CenterOfMovingIndicator+DistanceOfFiguresFromAxis, [255 255 255]);
            end
            Screen(w,'Flip');
            if isempty(eventTime) % recording the presentation start time
                eventTime = GetSecs-runStart;
            end
        end
    end
    
    %-----------------------------------------------------------------
    % show fixation ITI
    Screen('TextSize',w, betsz);
    CenterText(w,'+', white,0,FixForFixationCrossLocation);
    Screen(w,'Flip');
    if use_eyetracker
        Eyelink('Message',Eventflag(GenFlags.Fixation.str,Task,1,trial,runStart));
    end
    WaitSecs(0.3);
    
    %-----------------------------------------------------------------
    % write to output file
    
    fprintf(fid1,'%s\t%d\t%d\t%s\t%d\t%d\t%d \n', subjectID, trial, eventTime, stimuli_images(trial).name, bid, respTime, FirstMouseMovement);
end

%HideCursor;
Screen('TextSize',w, instrSZ);

Screen('Flip', w);
WaitSecs(3); % prevent key spillover

fclose(fid1);
toc
%-----------------------------------------------------------------
%	display end of part message
%-----------------------------------------------------------------

Screen('FillRect', w, black);
Screen('TextSize',w, 40);
CenterText(w,sprintf('Thank you! Please call the experimenter.') ,white,0,-170);
Screen('Flip',w);


%%   Finishing eye tracking  %
%---------------------------------------------------------------
if use_eyetracker
    
    %---------------------------
    % finish up: stop recording eye-movements,
    % close graphics window, close data file and shut down tracker
    Eyelink('StopRecording');
    %   Eyelink MSG
    % ---------------------------
    Eyelink('Message',['Eyetracking_closeTime: ',num2str(GetSecs-runStart)]);
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
        disp(rdf);
    end
    
    
    if dummymode==0
        movefile(edfFile,['./Output/', subjectID,'_',Task,'_eyetracking_', timestamp,'.edf']);
    end;
end



ShowCursor;
WaitSecs(5);
Screen('closeall');

end
