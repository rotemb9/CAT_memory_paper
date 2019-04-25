function disp_resolveBDM(subjectID, mainPath, sessionNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  compute and display BDM raffle - by Nadav Aridan
% Modified by Rotem Botvinik Nezer, March 2016
%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               .-------.
%                               |Jackpot|
%                   ____________|_______|____________
%                  |  __    __    ___  _____   __    |
%                  | / _\  / /   /___\/__   \ / _\   |
%                  | \ \  / /   //  //  / /\ \\ \ 25c|
%                  | _\ \/ /___/ \_//  / /  \/_\ \ []|
%                  | \__/\____/\___/   \/     \__/ []|
%                  |===_______===_______===_______===|
%                  ||*|\_     |*| _____ |*|\_     |*||
%                  ||*|| \ _  |*||     ||*|| \ _  |*||
%                  ||*| \_(_) |*||*BAR*||*| \_(_) |*||
%                  ||*| (_)   |*||_____||*| (_)   |*|| __
%                  ||*|_______|*|_______|*|_______|*||(__)
%                  |===_______===_______===_______===| ||
%                  ||*| _____ |*|\_     |*|  ___  |*|| ||
%                  ||*||     ||*|| \ _  |*| |_  | |*|| ||
%                  ||*||*BAR*||*| \_(_) |*|  / /  |*|| ||
%                  ||*||_____||*| (_)   |*| /_/   |*|| ||
%                  ||*|_______|*|_______|*|_______|*||_//
%                  |===_______===_______===_______===|_/
%                  ||*|  ___  |*|   |   |*| _____ |*||
%                  ||*| |_  | |*|  / \  |*||     ||*||
%                  ||*|  / /  |*| /_ _\ |*||*BAR*||*||
%                  ||*| /_/   |*|   O   |*||_____||*||
%                  ||*|_______|*|_______|*|_______|*||
%                  |lc=___________________________===|
%                  |  /___________________________\  |
%                  |   |                         |   |
%                 _|    \_______________________/    |_
%                (_____________________________________)
%
% run from main experiment folder
% requires "stim" and "output" folders

%% pre settings
if nargin < 1
    subjectID = input('Subject full code:','s'); % the BDM1 file name without the "_BDM1" part
end

if nargin < 2
    mainPath = pwd; % Change if you don't run from the experimnet folder
end

if nargin < 3
    sessionNum = 1;
end

outputPath = [mainPath '/Output'];
BackColor = 0; % Black background color
TextColor = 255; % White Text color

% open screen
keepTrying = 1;
PresentScreen = max(Screen('Screens'));
while keepTrying < 10
    try
        Window = Screen('OpenWindow', PresentScreen ,BackColor); % Opens the screen with the chosen backcolor
        Screen ('fillRect',Window,BackColor);
        Screen ('flip',Window);
        Screen('TextSize',Window,28);
        Screen('TextFont',Window,'Arial');
        keepTrying = 10;
    catch
        sca;
        keepTrying = keepTrying + 1;
        disp('CODE HAD CRASHED - cant open screen!');
    end
end



% Open data from BDM`
fid = fopen([outputPath '/' subjectID '_BDM' num2str(sessionNum) '.txt']);
BDMdata = textscan(fid, '%d %s %f %f', 'Headerlines',1);
fclose(fid);

total_num_items = length(BDMdata{1});
% load images
%! stimuli=dir([pwd '/stim/*.bmp' ]);
%! stimname=struct2cell(rmfield(stimuli,{'date','bytes','isdir','datenum'}));
stimname=BDMdata{2};
Images=cell(length(stimname),1);
for i=1:length(stimname)
    Images{i}=imread([ pwd '/stim/' stimname{i}]);
end


%%
% instruct: "the snack chosen is"
DrawFormattedText(Window,'the snack chosen is:','center',150,TextColor);
% Screen('PutImage',Window,instrct_snack);
Screen(Window,'Flip');

WaitSecs(2);
%% Choose a random item

% instruct: "the snack chosen is"
DrawFormattedText(Window,'the snack chosen is:','center',150,TextColor);
% Screen('PutImage',Window,instrct_snack);

taim= GetSecs;
while (GetSecs - taim) <= 2,
    
    DrawFormattedText(Window,'the snack chosen is:','center',150,TextColor);
    % Screen('PutImage',Window,instrct_snack);
    
    
    randomOrder = Shuffle(1:total_num_items);
    item_ind = randomOrder(1);
    
    % show the snack chosen
    Screen('PutImage',Window, Images{item_ind} );
    %  "your bid was:"
    
    item_chosen=char(stimname(item_ind));
    
    
    Screen('Flip',Window);
    WaitSecs(0.025)
    
end;



%%

% instruct: "the snack chosen is"
DrawFormattedText(Window,'the snack chosen is:','center',150,TextColor);
% Screen('PutImage',Window,instrct_snack);

% show the snack chosen
stimArray=imread([mainPath '/Stim/' item_chosen]);
Screen('PutImage',Window,stimArray );
Screen(Window,'Flip');

WaitSecs(1);


%%

% instruct: "the snack chosen is"
DrawFormattedText(Window,'the snack chosen is:','center',150,TextColor);
% Screen('PutImage',Window,instrct_snack);

% show the snack chosen
stimArray=imread([mainPath '/Stim/' item_chosen]);
Screen('PutImage',Window,stimArray );
%  "your bid was:"

item_bid = BDMdata{3}(item_ind);
item_bid=num2str(item_bid,'%.2f');

DrawFormattedText(Window,['your bid was: ' item_bid ''],'center',800,TextColor);
Screen(Window,'Flip');

WaitSecs(2);


%%
taim= GetSecs;
while (GetSecs - taim) <= 2,
    
    DrawFormattedText(Window,'the snack chosen is:','center',150,TextColor);
    % Screen('PutImage',Window,instrct_snack);
    
    % show the snack chosen
    stimArray=imread([mainPath '/Stim/' item_chosen]);
    Screen('PutImage',Window,stimArray );
    %  "your bid was:"
    
    
    DrawFormattedText(Window,['your bid was: ' item_bid ''],'center',800,TextColor);
    % show subject bid
    
    % Compute a random number
    possibleNumbers = Shuffle(0:0.5:10);
    chosenNumber = possibleNumbers(1);
    chosenNumber=num2str(chosenNumber,'%.2f');
    
    
    
    % show computer's bid
    DrawFormattedText(Window,'the computer bid is: ' ,'center',850,TextColor);
    DrawFormattedText(Window,[chosenNumber '' ],1100,850,TextColor);
    
    
    Screen('Flip',Window);
    WaitSecs(0.01)
    
end;

%%

DrawFormattedText(Window,'the snack chosen is:','center',150,TextColor);
% Screen('PutImage',Window,instrct_snack);

% show the snack chosen
stimArray=imread([mainPath '/Stim/' item_chosen]);
Screen('PutImage',Window,stimArray );
%  "your bid was:"


DrawFormattedText(Window,['your bid was: ' item_bid ''],'center',800,TextColor);
% show subject bid


% show computer's bid
DrawFormattedText(Window,['the computer bid is: ' chosenNumber ''],'center',850,TextColor);
Screen('Flip',Window);

WaitSecs(1);

%%
DrawFormattedText(Window,'the snack chosen is:','center',150,TextColor);
% Screen('PutImage',Window,instrct_snack);

% show the snack chosen
stimArray=imread([mainPath '/Stim/' item_chosen]);
Screen('PutImage',Window,stimArray );
%  "your bid was:"


DrawFormattedText(Window,['your bid was: ' item_bid ''],'center',800,TextColor);
% show subject bid


% show computer's bid
DrawFormattedText(Window,['the computer bid is: ' chosenNumber ''],'center',850,TextColor);


if str2double(chosenNumber) <= str2double(item_bid)
    
    DrawFormattedText(Window,['so you can buy it for ' chosenNumber ''],'center',900,[0 255 0]);
    Screen(Window,'Flip');
    
else
    DrawFormattedText(Window,'so you can not buy it','center',900,[255 0 0]);
    Screen(Window,'Flip');
end
WaitSecs(1);
%%
fid2 = fopen([outputPath '/' subjectID '_BDM_resolve_session_' num2str(sessionNum) '.txt'],'a');

% save results

if str2double(chosenNumber) <= str2double(item_bid)
    % The subject buys the item for chosenNumber NIS
    fprintf(fid2, 'You may buy item %s for price %s  your bid was %s\n', item_chosen, chosenNumber, item_bid);
else % the subject cannot buy the item
    fprintf(fid2, 'You do not recieve item %s Random number %s is higher than your bid %s\n', item_chosen, chosenNumber, item_bid);
end % end if chosenNumber <= item_bid

fclose(fid2);
%%
DrawFormattedText(Window,'the snack chosen is:','center',150,TextColor);
% Screen('PutImage',Window,instrct_snack);

% show the snack chosen
stimArray=imread([mainPath '/Stim/' item_chosen]);
Screen('PutImage',Window,stimArray );
%  "your bid was:"


DrawFormattedText(Window,['your bid was: ' item_bid ''],'center',800,TextColor);
% show subject bid


% show computer's bid
DrawFormattedText(Window,['the computer bid is: ' chosenNumber ''],'center',850,TextColor);


if str2double(chosenNumber) <= str2double(item_bid)
    
    DrawFormattedText(Window,['so you can buy it for ' chosenNumber ''],'center',900,[0 255 0]);
    
else
    
    DrawFormattedText(Window,'so you can not buy it','center',900,[255 0 0]);
end
DrawFormattedText(Window,'press any key to close','center',980,TextColor);
Screen(Window,'Flip');

while KbCheck; end
noresp=1;
while noresp,
    [keyIsDown] = KbCheck(-1); % deviceNumber=keyboard
    if keyIsDown && noresp,
        noresp=0;
    end;
end;
while KbCheck; end
sca
end
