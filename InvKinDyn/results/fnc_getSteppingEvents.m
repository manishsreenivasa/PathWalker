function [l_step,r_step] = fnc_getSteppingEvents(c3d_fileName,boolID)

h1 = btkReadAcquisition(c3d_fileName);
events = btkGetEvents(h1);
samplingFreq = btkGetPointFrequency(h1);
try
    l_step = ceil(events.Left_Foot_Strike*samplingFreq)-btkGetFirstFrame(h1);
    if length(l_step)<2
        l_step = [NaN NaN];
    end
    if boolID
        try
            btkGetPoint(h1,'LGroundReactionForce');
        catch
            l_step = [NaN NaN];
        end
    end
catch
    l_step = [NaN NaN];
end

try
    r_step = ceil(events.Right_Foot_Strike*samplingFreq)-btkGetFirstFrame(h1);
    if length(r_step)<2
        r_step = [NaN NaN];
    end
    if boolID
        try
            btkGetPoint(h1,'RGroundReactionForce');
        catch
            r_step = [NaN NaN];
        end
    end
catch
    r_step = [NaN NaN];
end