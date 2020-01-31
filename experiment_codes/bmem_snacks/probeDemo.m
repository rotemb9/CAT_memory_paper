function [] = probeDemo(subjectID, order, mainPath, isMRI, sessionNum)
% function [] = probeDemo(subjectID, order, mainPath, isMRI, sessionNum)

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
% =============== Created based on the previous boost codes ===============
% ==================== by Rotem Botvinik November 2014 ====================
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% This function runs the demo of the probe session of the cue-approach task.

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % --------- Exterior files needed for task to run correctly: ----------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   'Misc/demo_items.txt'


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % % ------------------- Creates the following files: --------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   ''probe_demo_' timestamp '.txt''
%


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % ------------------- dummy info for testing purposes -------------------
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% subjectID = 'bmem_snacks_999';
% subjectID = 'bmem_snacks_998'; % it is important to test both orders
% order = 1;
% order = 2;
% isMRI = 0;
% sessionNum = 1;
% mainPath = pwd;

tic

rng shuffle


%---------------------------------------------------------------
%% 'GLOBAL VARIABLES'
%---------------------------------------------------------------

%   'timestamp'
% - - - - - - - - - - - - - - - - -
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp=[date,'_',hr,'h',minutes,'m'];

outputPath = [mainPath '/Output'];

fid_demoNames=fopen('Misc/demo_items.txt');
demoNames=textscan(fid_demoNames,'%s');
demoNames = demoNames{1};
fclose(fid_demoNames);


%   'set phase times'
% - - - - - - - - - - - - - - - - -
maxtime = 1.5;      % 1.5 second limit on each selection
baseline_fixation_dur = 2; % Need to modify based on if first few volumes are saved or not
% afterrunfixation = 6;

tic

% -----------------------------------------------
%% Load Instructions
% -----------------------------------------------

% Load Hebrew instructions image files
Instructions = dir([mainPath '/Instructions/probe_demo.JPG' ]);
Instructions_image = imread([mainPath '/Instructions/' Instructions(1).name]);


%---------------------------------------------------------------
%% 'INITIALIZE Screen variables's
%---------------------------------------------------------------
Screen('Preference', 'VisualDebuglevel', 3); %No PTB intro screen
screennum = max(Screen('Screens'));

pixelSize = 32;
% [w] = Screen('OpenWindow',screennum,[],[0 0 640 480],pixelSize/4); % debugging screensize
[w] = Screen('OpenWindow',screennum,[],[],pixelSize);

% Here Be Colors
% - - - - - - - - - - - - - - - - -
black = BlackIndex(w); % Should equal 0.
white = WhiteIndex(w); % Should equal 255.
yellow = [0 255 0];



% set up Screen positions for stimuli
% - - - - - - - - - - - - - - - - -
[wWidth, wHeight] = Screen('WindowSize', w);
xcenter = wWidth/2;
ycenter = wHeight/2;

Screen('FillRect', w, black);  % NB: only need to do this once!
Screen('Flip', w);


% setting the text
% - - - - - - - - - - - - - - - - -
theFont = 'Arial';
Screen('TextFont',w,theFont);
Screen('TextSize',w, 28);


% stimuli locations
% - - - - - - - - - - - - - - - - -
stimW = 576;
stimH = 432;

distcent = 300;
leftRect = [xcenter-stimW-distcent ycenter-stimH/2 xcenter-distcent ycenter+stimH/2];
rightRect = [xcenter+distcent ycenter-stimH/2 xcenter+stimW+distcent ycenter+stimH/2];

penWidth = 10;

HideCursor;


%---------------------------------------------------------------
%% 'Assign response keys'
%---------------------------------------------------------------

KbName('UnifyKeyNames');

if isMRI == 1
    leftstack = 'b';
    rightstack = 'y';
    badresp = 'x';
else
    leftstack = 'u';
    rightstack = 'i';
    badresp = 'x';
    
end % end if isMRI == 1


%---------------------------------------------------------------
%% 'Load image arrays'
%---------------------------------------------------------------

demofood_items = cell(1, length(demoNames));
for i = 1:length(demoNames)
    demofood_items{i} = imread([mainPath sprintf('/stim/demo/%s',demoNames{i})]);
end
%


%   'load onsets'
% - - - - - - - - - - - - - - -

r = Shuffle(1:2); %re-shuffled done every run
onsetlist = load(['Onset_files/probe_demo_onset_' num2str(r(1)) '.mat']);
onsetlist = onsetlist.onsetlist;


%---------------------------------------------------------------
%%   'PRE-TRIAL DATA ORGANIZATION'
%---------------------------------------------------------------

% Determine stimuli to use
% - - - - - - - - - - - - - - - - -


HH_HS = [1 2 3 4];
HH_HG = [5 6 7 8];


leftGo = cell(1,length(HH_HS)+1);
stimnum1 = cell(1,length(HH_HS)+1);
stimnum2 = cell(1,length(HH_HS)+1);
leftname = cell(1,length(HH_HS)+1);
rightname = cell(1,length(HH_HS)+1);
pairtype = cell(1,length(HH_HS)+1);

for block = 1:1
    pairtype{block} = [1 1 1 1 1 1 1 1];
    leftGo{block} = Shuffle([1 1 0 0 1 1 0 0]);
end;

HH=1;
% LL=1;
% HL_S=1;
% HL_G=1;


for i=1:length(HH_HS) % trial num within block
    stimnum1{block}(i)=HH_HS(HH);
    stimnum2{block}(i)=HH_HG(HH);
    HH=HH+1;
    if leftGo{block}(i)==1
        leftname{block}(i)=demoNames(stimnum1{block}(i));
        rightname{block}(i)=demoNames(stimnum2{block}(i));
    else
        leftname{block}(i)=demoNames(stimnum1{block}(i));
        rightname{block}(i)=demoNames(stimnum2{block}(i));
    end
    
end % end for i=1:4

%ListenChar(2); %suppresses terminal ouput
KbQueueCreate;
%---------------------------------------------------------------
%% 'Write output file header'
%---------------------------------------------------------------

fid1=fopen([outputPath '/' subjectID sprintf('_probe_demo_') timestamp '_session' num2str(sessionNum) '.txt'], 'a');
fprintf(fid1,'subjectID\tscanner\torder\truntrial\tonsettime\tImageLeft\tImageRight\tTypeLeft\tTypeRight\tIsleftGo\tResponse\tPairType\tOutcome\tRT\n'); %write the header line

%---------------------------------------------------------------
%% 'Display Main Instructions'
%---------------------------------------------------------------

Screen('TextSize',w, 40);

%%% While they are waiting for the trigger
if isMRI    
    Screen('PutImage',w,Instructions_image_fmri);
    Screen(w,'Flip');
    
    % wait for the subject to press the button
    noresp = 1;
    while noresp,
        [keyIsDown,~,~] = KbCheck; %(experimenter_device);
        if keyIsDown && noresp,
            noresp = 0;
        end;
    end;
    CenterText(w,'This is JUST A DEMO', white, 0, -50);
    CenterText(w,'Waiting for trigger...GET READY....', white, 0, 50);
    Screen(w,'Flip');
    
    escapeKey = KbName('t');
    while 1
        [keyIsDown,~,keyCode] = KbCheck(-1);
        if keyIsDown && keyCode(escapeKey)
            break;
        end
    end
    
else
    Screen('PutImage',w,Instructions_image);
    Screen('Flip',w);
    
    noresp = 1;
    while noresp,
        [keyIsDown] = KbCheck(-1); % deviceNumber=keyboard
        if keyIsDown && noresp,
            noresp = 0;
        end;
    end;
    
    WaitSecs(0.001);
    
end % end if isMRI == 1

prebaseline = GetSecs;
KbQueueCreate;


%-----------------------------------------------------------------

while GetSecs < prebaseline + baseline_fixation_dur
    CenterText(w,'+', white,0,0);
    Screen('TextSize',w, 60);
    Screen(w,'Flip');
    
end
%-----------------------------------------------------------------
% postbaseline = GetSecs;
% baseline_fixation = postbaseline - prebaseline;


%---------------------------------------------------------------
%% 'Run Trials'
%---------------------------------------------------------------
runtrial = 1;
runStart = GetSecs;



for block=1:1
    for trial=1:length(HH_HS)
        
        colorleft=black;
        colorright=black;
        out=999;
        %-----------------------------------------------------------------
        % display images
        if leftGo{block}(trial)==1
            Screen('PutImage',w,demofood_items{stimnum1{block}(trial)}, leftRect);
            Screen('PutImage',w,demofood_items{stimnum2{block}(trial)}, rightRect);
        else
            Screen('PutImage',w,demofood_items{stimnum2{block}(trial)}, leftRect);
            Screen('PutImage',w,demofood_items{stimnum1{block}(trial)}, rightRect);
        end
        CenterText(w,'+', white,0,0);
        Screen(w,'Flip', runStart+onsetlist(runtrial));
        StimOnset = GetSecs;
        KbQueueFlush;
        KbQueueStart;
        
        %-----------------------------------------------------------------
        % get response
        
        
        noresp = 1;
        goodresp = 0;
        while noresp
            % check for response
            [keyIsDown, firstPress] = KbQueueCheck;
            
            
            if keyIsDown && noresp
                keyPressed = KbName(firstPress);
                if ischar(keyPressed)==0 % if 2 keys are hit at once, they become a cell, not a char. we need keyPressed to be a char, so this converts it and takes the first key pressed
                    keyPressed = char(keyPressed);
                    keyPressed = keyPressed(1);
                end
                switch keyPressed
                    case leftstack
                        respTime = firstPress(KbName(leftstack))-StimOnset;
                        noresp = 0;
                        goodresp = 1;
                    case rightstack
                        respTime = firstPress(KbName(rightstack))-StimOnset;
                        noresp = 0;
                        goodresp = 1;
                end
            end % end if keyIsDown && noresp
            
            
            % check for reaching time limit
            if noresp && GetSecs-runStart >= onsetlist(runtrial)+maxtime
                noresp = 0;
                keyPressed = badresp;
                respTime = maxtime;
            end
        end % end while noresp
        
        
        %-----------------------------------------------------------------
        
        
        % determine what bid to highlight
        
        switch keyPressed
            case leftstack
                colorleft=yellow;
                if leftGo{block}(trial)==0
                    out=0;
                else
                    out=1;
                end
            case rightstack
                colorright=yellow;
                if leftGo{block}(trial)==1
                    out=0;
                else
                    out=1;
                end
        end % end switch keyPressed
        
        if goodresp==1
            if leftGo{block}(trial)==1
                Screen('PutImage',w,demofood_items{stimnum1{block}(trial)}, leftRect);
                Screen('PutImage',w,demofood_items{stimnum2{block}(trial)}, rightRect);
            else
                Screen('PutImage',w,demofood_items{stimnum2{block}(trial)}, leftRect);
                Screen('PutImage',w,demofood_items{stimnum1{block}(trial)}, rightRect);
            end
            Screen('FrameRect', w, colorleft, leftRect, penWidth);
            Screen('FrameRect', w, colorright, rightRect, penWidth);
            CenterText(w,'+', white,0,0);
            Screen(w,'Flip',runStart+onsetlist(trial)+respTime);
            
        else % didn't respond in time
            Screen('DrawText', w, 'You must respond faster!', xcenter-400, ycenter, white);
            Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime);
        end
        
        
        %-----------------------------------------------------------------
        % show fixation ITI
        CenterText(w,'+', white,0,0);
        Screen(w,'Flip',runStart+onsetlist(runtrial)+respTime+.5);
        
        if goodresp ~= 1
            respTime=999;
        end
        
        %-----------------------------------------------------------------
        % save data to .txt file
        fprintf(fid1,'%s\t %d\t %d\t %d\t %d\t %s\t %s\t %d\t %d\t %d\t %s\t %d\t %d\t %.2f\t \n', subjectID, isMRI, order, runtrial, onsetlist(runtrial), char(leftname{block}(trial)), char(rightname{block}(trial)), stimnum1{block}(trial), stimnum2{block}(trial), leftGo{block}(trial), keyPressed, pairtype{block}(trial), out, respTime*1000);
        runtrial = runtrial+1;
%         KbQueueFlush;
    end % end for trial = 1:4
end % end for block = 1:1

WaitSecs(2);
Screen('FillRect', w, black);
Screen('TextSize',w, 40);

% if sessionNum == 1
%     CenterText(w,sprintf('Please read the Part 5 instructions and continue on your own.') ,white,0,-270);
% else
%     CenterText(w,sprintf('Please read the Part 2 instructions and continue on your own.') ,white,0,-270);
% end

if isMRI == 1 % edit this part so that in the fMRI the begining of the next part would depend on the experimenter and not on the subject
    CenterText(w,sprintf('The demo is done. Questions?') ,white,0,-170);
    Screen('Flip',w);
    WaitSecs(5);
else
    CenterText(w,sprintf('The demo is done.') ,white,0,-270);
    CenterText(w,sprintf('Press any key when you are ready to continue.') ,white,0,-170);
    Screen('Flip',w);
    noresp=1;
    while noresp,
        [keyIsDown,~,~] = KbCheck;
        if keyIsDown && noresp,
            noresp=0;
        end;
    end;
    
end % end if isMRI == 1



fclose(fid1);

Screen('CloseAll');
ShowCursor;

outfile=strcat(outputPath, sprintf('/%s_probe_demo_%s_session%d.mat', subjectID, timestamp, sessionNum));

% create a data structure with info about the run
run_info.subject = subjectID;
run_info.date = date;
run_info.outfile = outfile;

run_info.script_name = mfilename;
clear demofood_items Instruction*;
save(outfile);

toc

end % end function

