function personal_details(subjectID, order, outputPath, sessionNum)

% function personalDetails(subjectID, order, mainPath)
%   This function gets a few personal details from the subject and saves it
%   to file named subjectID '_personalDetails' num2str(sessionNum) '.txt'

% get time and date
c = clock;
hr = sprintf('%02d', c(4));
minutes = sprintf('%02d', c(5));
timestamp = [date,'_',hr,'h',minutes,'m'];


% open a txt file for the details
fid1 = fopen([outputPath '/' subjectID '_personalDetails' num2str(sessionNum) '_' timestamp '.txt'], 'a');
fprintf(fid1,'subjectID\torder\tdate\tgender(1-female, 2-male)\tage\tdominant hand (1-right, 2-left)\theight\tweight\toccupation\tlast food\thunger level\n'); %write the header line

% ask the subject for the details

%set text size of the dialod box
set(groot,'defaultUicontrolFontSize', 16)
Gender = questdlg('Please select your gender:','Gender','Female','Male','Female');
% options.Resize='on';


while isempty(Gender)
    Gender = questdlg('Please select your gender:','Gender','Female','Male','Female');
end
if strcmp(Gender,'Male')
    Gender = 2;
else
    Gender = 1;
end


%set text size of the dialod box
set(groot,'defaulttextfontsize',18);
Age = myinputdlg('Please enter your age: ','Age',1);
while isempty(Age) || isempty(Age{1})
    Age = myinputdlg ('Only integers between 18 and 40 are valid. Please enter your age: ','Age',1);
end
Age = cell2mat(Age);
Age = str2double(Age);
while mod(Age,1) ~= 0 || Age < 18 || Age > 40
    Age = myinputdlg ('Only integers between 18 and 40 are valid. Please enter your age: ','Age',1);
    Age = cell2mat(Age);
    Age = str2double(Age);
end


%set text size of the dialod box
set(groot,'defaultUicontrolFontSize', 16)
DominantHand = questdlg('Please select your domoinant hand:','Dominant hand','Left','Right','Right');
while isempty(DominantHand)
    DominantHand = questdlg('Please select your domoinant hand:','Dominant hand','Left','Right','Right');
end
if strcmp(DominantHand,'Left')
    DominantHand = 2;
else
    DominantHand = 1;
end

%set text size of the dialod box
set(groot,'defaultUicontrolFontSize', 16)
Height = myinputdlg ('Please type your height in Cm (for example 168): ','Height',1);
while isempty(Height) || isempty(Height{1})
    Height = myinputdlg ('Please correct. Enter your height: ','Height',1);
end
Height = cell2mat(Height);
Height = str2double(Height);
while mod(Height,1) ~= 0 || Height < 100 || Height > 230
    Height = myinputdlg ('Please correct. Enter your height: ','Height',1);
    Height = cell2mat(Height);
    Height = str2double(Height);
end

%set text size of the dialod box
set(groot,'defaultUicontrolFontSize', 16)
Weight = myinputdlg ('Please type your weight in Kg (for example 58):','Weight',1);
while isempty(Weight) || isempty(Weight{1})
    Weight = myinputdlg ('Please correct. Enter your Weight: ','Weight',1);
end
Weight = cell2mat(Weight);
Weight = str2double(Weight);
while mod(Height,1) ~= 0 || Weight < 30 || Weight > 200
    Weight = myinputdlg ('Please correct. Enter your Weight: ','Weight',1);
    Weight = cell2mat(Weight);
    Weight = str2double(Weight);
end

%set text size of the dialod box
set(groot,'defaultUicontrolFontSize', 16)
Occupation = myinputdlg('Please type your occupation (for example- a student for Psychology): ','Occupation',1);
Occupation = cell2mat(Occupation);


%set text size of the dialod box
set(groot,'defaultUicontrolFontSize', 16)
Eat = myinputdlg('How many hours ago did you eat? ','Eat',1);
Eat = cell2mat(Eat);


%set text size of the dialod box
set(groot,'defaultUicontrolFontSize', 16)
Hungry = myinputdlg('How hungry are you right now from 1 to 10? 1-not at all, 10- very much): ','Hungry',1);
while isempty(Hungry) || isempty(Hungry{1})
    Hungry = myinputdlg ('Please correct. How hungry are you right now from 1 to 10? 1-not at all, 10- very much): ','Hungry',1);
end
Hungry = cell2mat(Hungry);
Hungry = str2double(Hungry);
while mod(Hungry,1) ~= 0 || Hungry < 0 || Hungry > 10
    Hungry = myinputdlg ('Please correct. How hungry are you right now? 1-not at all, 10- very much): ','Hungry',1);
    Hungry = cell2mat(Hungry);
    Hungry = str2double(Hungry);
end


%% finish this part
%set text size of the dialod box
set(groot,'defaultUicontrolFontSize', 24)
questdlg('Thank you!!! Please call the experimenter','This part is over','Continue','Continue');
WaitSecs(1);

% close experiment window, and then start the ranking manually
Screen('CloseAll');
ShowCursor;

% Write details to file
fprintf(fid1,'%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\n', subjectID, order, timestamp, Gender, Age, DominantHand,Height, Weight, Occupation, Eat, Hungry);
fclose(fid1);
end