classdef GenFlags
    enumeration
        % Tasks
        BDM
        BDMDemo
        Training
        TrainingDemo
        Probe
        ProbeDemo
        Memory
        MemoryDemo
        BinaryRanking
        BinaryRankingDemo
        ResponseToStimuli
        ResponseToStimuliDemo
        FaceLocalizer
        FamiliarityDemo
        Familiarity
        ScaleRanking
        
        % for all:
        RunStart
        TrialStart
        Fixation
        Response
        RunEnd
        
        % for training:
        CueStart
        
        % for BDM Scale Ranking
        ScaleStart
       
        % for probe and recognition:
        RespondFaster
        
        % for recognition:
        PresentIsOld
        ResponseIsOld
        PresentIsGo
        ResponseIsGo
                
    end
    
    methods
        function text = str(GenFlags)
            text = char(GenFlags);
        end
    end
end


%%
