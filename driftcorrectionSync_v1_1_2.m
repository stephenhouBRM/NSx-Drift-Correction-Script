%%% Drift correction/alignment script for non-synchronized datafiles
% Workfile:         driftcorrectionSync_v1-1.m
% Author:           Stephen Hou
% Version:          1.1.2
% Purpose:          Determine offset/alignment off of sync pulses into two
%                   NSPs, then align and correct drift in data files
% Changelog:        v 1.0.0 - Initial Release
%                   v 1.1.0  - remove filesaving/buffering features,
%                   refined sync pulse alignment, add post-correction 
%                   segment length validation
%                   v 1.1.1- fixed a few edge cases regarding calculating
%                   total drift
%                   v 1.1.2- rewrote how overall drift is calculated
%% edit these.
%Ensure that the sync1 and signal1 refer to files from NSP1 and
%sync2/signal2 refer to files from NSP2

sync1filepath = 'D:\case5\QaHkESkEkA-20181029-115113-INST0.ns3';
sync2filepath = 'D:\case5\QaHkESkEkA-20181029-115113-INST1.ns3';
signal1filepath = 'D:\case5\QaHkESkEkA-20181029-115113-INST0.ns3';
signal2filepath = 'D:\case5\QaHkESkEkA-20181029-115113-INST1.ns3';

% sync1filepath = 'D:\case4\YmCiLtiHHid-20180919-135000-INST0.ns3';
% sync2filepath = 'D:\case4\YmCiLtiHHid-20180919-135000-INST1.ns3';
% signal1filepath = 'D:\case4\YmCiLtiHHid-20180919-135000-INST0.ns3';
% signal2filepath = 'D:\case4\YmCiLtiHHid-20180919-135000-INST1.ns3';
 
% sync1filepath = 'D:\180901_1135_ieeg1.ns3';
% sync2filepath = 'D:\180901_1135_ieeg2.ns3';
% signal1filepath = 'D:\180901_1135_ieeg1.ns3';
% signal2filepath = 'd:\180901_1135_ieeg2.ns3';


syncChan = 4;  %INPUT WHAT ROW THE SYNC SIGNAL IS IN

%% Open files, check validity
NSP1 = openNSx(sync1filepath,strcat('c:',int2str(syncChan)));
NSP2 = openNSx(sync2filepath,strcat('c:',int2str(syncChan)));

syncFs = NSP1.MetaTags.SamplingFreq;
if syncFs ~= NSP2.MetaTags.SamplingFreq
    fprintf('ERROR: Sync pulse files at different sampling frequencies');
    return
end    

% sync1 = NSP1.Data;
% sync2 = NSP2.Data;

signalFs = checkFs(signal1filepath);
if signalFs > syncFs
    fprintf('ERROR: Sync pulse sampled at lower frequency than neural data')
    return
end

sync1 = downsample(NSP1.Data, syncFs/signalFs);
sync2 = downsample(NSP2.Data, syncFs/signalFs);

clear NSP1 NSP2

%% identify sync pulse

edgeThresh = 5000;
%find first rising edge in both sync files
idxRisingEdge1 = checkEdge(sync1,edgeThresh,'rising');
idxRisingEdge2 = checkEdge(sync2,edgeThresh,'rising');
%find last falling edge in both sync files
idxFallingEdge1 = checkEdge(sync1,edgeThresh,'falling');
idxFallingEdge2 = checkEdge(sync2,edgeThresh,'falling');

%% identify the file that has a faster clock rate, open the file that has a slower clock rate
if idxFallingEdge1-idxRisingEdge1 > idxFallingEdge2 - idxRisingEdge2
    fasterFile = 1;
    signalToCorrect = openNSx(signal2filepath);
else
    fasterFile = 2;
    signalToCorrect = openNSx(signal1filepath);
end
    
%% identify the file that started recording first

if idxRisingEdge1 > idxRisingEdge2
    firstFile=1;
elseif idxRisingEdge1 < idxRisingEdge2
    firstFile=2;
end

%% set a flag if the faster clock NSP did not start recording first
if firstFile ~= fasterFile
    driftOverlap=1;
else
    driftOverlap=0;
end

%% Validate sync pulse alignment

%check if the segments corresponding to the edge indices are identical 
segLength = 1240;
rRise = corrcoef(double(sync2(idxRisingEdge2:idxRisingEdge2+segLength)),double(sync1(idxRisingEdge1:idxRisingEdge1+segLength)));
rFall = corrcoef(double(sync2(idxFallingEdge2-segLength:idxFallingEdge2)),double(sync1(idxFallingEdge1-segLength:idxFallingEdge1)));

if rFall(1,2) && rRise(1,2) < .999
    fprintf('ERROR: Mismatched sync pulses');
    return
end

%% if necessary, convert indices 

% if signalToCorrect.MetaTags.SamplingFreq ~=syncFs
%     [idxRisingEdge1, idxRisingEdge2,idxFallingEdge1,idxFallingEdge2] = convertFs(idxRisingEdge1, idxRisingEdge2,idxFallingEdge1,idxFallingEdge2, syncFs, signalToCorrect.MetaTags.SamplingFreq);
% end

%% correct drift
corrected = correctDrift(signalToCorrect.Data,idxRisingEdge1, idxRisingEdge2,idxFallingEdge1,idxFallingEdge2,driftOverlap);

%% open faster-clock file, check the sync-pulse segment is the same length as the corrected signal
switch fasterFile
    case 1
        templateSignal = openNSx(signal1filepath);
        uncorrected = templateSignal.Data(:,idxRisingEdge1:idxFallingEdge1);
    case 2
        templateSignal = openNSx(signal2filepath);
        uncorrected = templateSignal.Data(:,idxRisingEdge2:idxFallingEdge2);
end

if max(size(uncorrected))~= max(size(corrected))
    fprintf('ERROR MISMATCHED SEGMENT LENGTH');
% else
%     fprintf('CORRELATION COEFFICIENT FOR CORRECTED SYNC CHANNEL:\n')
%     disp(corrcoef(double(corrected(4,:)),double(uncorrected(4,:))))
end


%% helper functions

%%correct drift
function [corrected] = correctDrift(data,idxRE1,idxRE2,idxFE1,idxFE2,driftOverlap)
startOffset = abs(idxRE1-idxRE2);
endOffset = abs(idxFE1-idxFE2);
seg1Length = idxFE1-idxRE1;
seg2Length = idxFE2-idxRE2;

drift = abs(seg1Length-seg2Length);

numChans = min(size(data));
%startBuf = zeros(1,startOffset);
corrected = zeros(numChans,min(idxFE1-idxRE1,idxFE2-idxRE2)+drift+1);

correctionInterval = floor(min(seg1Length,seg2Length)/drift);

repN = ones(1,min(seg1Length,seg2Length)+1);

idx = zeros(length(drift));
for i=1:drift
    idx(i) = (correctionInterval+1)+correctionInterval*(i-1);
end

repN(idx) = 2;

if seg1Length < seg2Length
    slowFile=1;
    segLength = seg1Length;
else
    slowFile=2;
    segLength = seg2Length;
end

if slowFile == 1
    for i=1:numChans
        corrected(i,:) = repelem(data(i,idxRE1:idxRE1+segLength),repN);
    end
elseif slowFile ==2
    for i=1:numChans
        corrected(i,:) = repelem(data(i,idxRE2:idxRE2+segLength),repN);
    end
end
% for i=1:numChans
%     %temp = repelem(data(i,min(idxRE2,idxRE1):end),repN);
%     corrected(i,:) =  repelem(data(i,min(idxRE2,idxRE1):min(idxFE2,idxFE1)),repN);
%   %  corrected(i,:) =  repelem(data(i,min(idxRE2,idxRE1):max(idxFE2,idxFE1)),repN);
% 
%     % if length(corrected(i,:)) ~= length(horzcat(data(i,1:min(idxRE2,idxRE1)-1),temp))
% %        fprintf('ERROR');
%  %     return;
% end 
% corrected(i,:) = horzcat(data(i,1:min(idxRE2,idxRE1)-1),temp);
end


%check for rising/falling edges
function [index] = checkEdge(data, thresh, edgeType)

if strcmp(edgeType, 'falling') 
    for i=length(data):-1:1
        if data(i) > thresh
            index = i+1;
            return;
        end
    end
end
if strcmp(edgeType, 'rising')
    for i=1:1:length(data)
        if data(i) > thresh
            index = i-1;
            return;
        end
    end
end
end

function [Fs] = checkFs(filename)
    ext = filename(length(filename)-3:length(filename));
    if strcmp('.ns1',ext)
        Fs = 500;
    elseif strcmp('.ns2',ext)
        Fs = 1000;
    elseif strcmp('.ns3',ext)
        Fs = 2000;
    elseif strcmp('.ns4',ext)
        Fs = 10000;
    elseif strcmp('.ns5',ext)
        Fs = 30000;
    elseif strcmp('.ns6', ext)
        Fs= 30000;
    end
end
% 
% function [signalRE1, signalRE2, signalFE1, signalFE2] = convertFs(RE1, RE2, FE1, FE2, syncFs, signalFs)
%     signalRE1 = floor(RE1/syncFs * signalFs);
%     signalRE2 = floor(RE2/syncFs * signalFs);
%     signalFE1 = floor(FE1/syncFs * signalFs);
%     signalFE2 = floor(FE2/syncFs * signalFs);
% end