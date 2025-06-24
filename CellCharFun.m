%% CellCharFun - Patch Clamp Data Analysis Script for step current spiking data
% This script processes electrophysiological data from CFS files.
% It extracts membrane properties, spike features and adaptation metrics,
% and includes analysis of initial burst frequencies used in identification
% of CA3 pyramidal neuron subclasses.
%
% This code has been written for analysis of recording files with specific
% structure, and is not intended as a general analysis script. I hope that
% you find it helpful, but please ensure all parameters and analyses are
% performed correctly and appropriately for your own data structure.
%
% This script assumes two negative current injections, a sweep at 0 nA injection,
% and subsequent positive current injections at a regular step size. It
% also restricts analysis to the first 11 sweeps for current injection
% range of -100 to +400 pA in 50 pA steps.
% Should alternative current step sizes or numbers be used, 
% the code should be thoroughly checked for data correspondence at
% all stages. If current injections of length very different to 1 s, all
% analyses should be thoroughly checked.
%
% Four output variables are produced:
% D - recorded data structured into D.param (Header paramters) and D.data
%
% t - A vector of time, for plotting D.data with a time axis (in seconds)
%
% Out_Spiking - This is a matrix of Initial Burst Frequency data. Each column
% contains the data from one channel, with the first row noting the recording
% channel number (for multicellular recordings).
% Rows 2-12 are the frequency of firing calculated from the first 3 spikes of
% each sweep. Rows 13-23 are the frequency of firing calculated from the
% first 2 spikes of each sweep.
%
% Out_Transposed - This is a matrix of calculated cell parameters, with one
% column for each channel. Data are arranged as follows
% Out_Transposed Row 1: Channel number, extracted from header information
% Out_Transposed Row 2: Membrane potential (median value of zero current input). Units - mV
% Out_Transposed Row 3: Max firing frequency achieved across full 1 s inputs. Units - Hz
% Out_Transposed Row 4: Current injection required for spike generation. Units - pA
% Out_Transposed Row 5: Input resistance. Units - MΩ
% Out_Transposed Row 6: Spike frequency adaptation. No units
% Out_Transposed Row 7-17: The spike frequency across 1 s current
% injections for sweeps 1-11 of recorded data.
%
% This code requires the cfs2mat script. This can be downloaded from: https://ced.co.uk/upgrades/sigimpcfs
%
% An example file (20220209_000.cfs) is included on the Github repository
% for this script. This file has 7 channels and 20 datasweeps. 
%
% Written by Jake Watson (2025), for Watson et al. Cell Reports 2025.

%% USER INPUT
% Please adjust step size to the magnitude of current injection in nA. Default is 0.05 nA.
stepsize = 0.05;
% Please insert the length of current injection step in seconds.
steplength = 1;
% Please insert the start time of current injection step in seconds.
stepstart = 0.25;

%% --- Initialization and Data Import
% Import CFS and convert to MAT format
cfs2mat
outNamemat = strcat(outName,'.mat');
open(outNamemat)

% Extract time vector
samplength = D.param.xScale(1);
samprate = 1/samplength;
t = 0:samplength:(samplength * size(D.data,1)-samplength);

%% --- Extracting channel names and numbers
channelNames = D.param.channelName';
fnamelist = repmat({fName}, D.param.channels, 1);
channelNums = str2double(regexp(channelNames, '[0-9]+', 'match', 'once'));

%% --- Step Parameter Generation
vmsteps = -(stepsize*2):stepsize:((D.param.dataSections - 3) * stepsize); % voltage step levels
Vm = squeeze(median(D.data(:,3,:),1)); % median voltage of zero injection sweep is assumed to be resting potential

%% --- Input Resistance Calculation
inbase = squeeze(median(D.data(1:(stepstart*samprate),[1:2 4],:),1)); % baseline is set as the median of pre-injection trace
startIdx = round((steplength + stepstart - 0.35) * samprate);
endIdx   = round((steplength + stepstart - 0.05) * samprate);
intest   = squeeze(median(D.data(startIdx:endIdx, [1:2 4], :), 1)); % response is set as the median of 300 ms towards the end of current injection period
inres = median(-(inbase - intest) ./ vmsteps([1:2 4])', 1)'; % Estimated input resistance from 2 negative and one postitive current injection sweep.

%% --- Spike Detection and Frequency Metrics
% Initialise matrices for spike analysis
% D.param.dataSections is the number of sweeps, extracted from recording header.
% D.param.channels is the number of reocording channels for multi-cell data.
% spikefreq - measures the total number of spikes across the injection pulse
% initialfreq - measures the frequency of the first three spikes 
% initialfreq2 - measures the frequency of the first two spikes on injection
% sfa - reports Spike Frequency Adaptation from the first sweep with 6 spikes

spikefreq       = zeros(D.param.dataSections,D.param.channels);
initialfreq     = zeros(D.param.dataSections,D.param.channels);
initialfreq2    = zeros(D.param.dataSections,D.param.channels);
sfa             = zeros(D.param.dataSections,D.param.channels);

% Analysis runs by nested loops of Channels and Sweeps
% Spikes here are detected as peaks rising above 'MinPeakHeight', with a
% refractory period of 'MinPeakDistance' between peaks. Please ensure that
% both of these parameters are appropriate for every recording. If current
% plateau on current injection rises above 'MinPeakHeight', spurious spike
% frequencies will be generated.

for i = 1:D.param.channels
    for j = 1:D.param.dataSections
        startIdx = round(stepstart * samprate);
        endIdx   = round((stepstart + steplength) * samprate);
        [spikemax, spikepos] = findpeaks(D.data(startIdx:endIdx, j, i), ...
        samprate,'MinPeakHeight', -10,'MinPeakDistance', 0.005);
        spikefreq(j,i) = size(spikepos,1);
        if size(spikepos,1) >= 3
            initialfreq(j,i) = 1/(spikepos(3)-spikepos(1));
        else
            initialfreq(j,i) = 0;
        end
             
        if size(spikepos,1) >= 2
            initialfreq2(j,i) = 1/(spikepos(2)-spikepos(1));
        else
            initialfreq2(j,i) = 0;
        end
        
        if size(spikepos,1) >=5
            earlysfa = ((spikepos(3) - spikepos(2)) + (spikepos(2) - spikepos(1))) ./2;
            latesfa = ((spikepos(end) - spikepos(end-1)) + (spikepos(end-1) - spikepos(end-2))) ./2;
            
            sfa(j,i) = earlysfa/latesfa;
        else
            sfa(j,i) = 0;
        end
     end
end

maxfreq = max(spikefreq(1:11,:),[],1)';

%% --- First Spike Metrics
% minI - reports the minimum current required for spike generation in pA
spikethere = zeros(D.param.channels,1);
minI = NaN(D.param.channels,1);
spikethresh = NaN(D.param.channels,1);

for n = 1:D.param.channels
    if stepsize == 0.05
        % Find first sweep with spikes (starting from step 3)
        idx = find(any(spikefreq(3:end,n) ~= 0, 2), 1);
        spikethere(n) = idx + 2; % offset by 2
    end
                
    if spikethere(n) ~= 0
        minI(n) = vmsteps(spikethere(n)) * 1000;
    end
end

%% --- SFA Metrics Summary
sfafirst = zeros(D.param.channels,1);

for p = 1:D.param.channels
    if any(sfa(:,p))
        firstIdx = find(sfa(:,p), 1);
        sfafirst(p) = sfa(firstIdx, p);
    end
end

%% --- Tabulate Stepwise Results
vmstep_labels = compose('%.2f', vmsteps);
spikeLabels = strcat("SpikeFreq_", vmstep_labels);
initLabels  = strcat("InitFreq_", vmstep_labels);

SpikeFreqT = array2table(spikefreq', 'VariableNames', spikeLabels);
InitFreqT = array2table(initialfreq', 'VariableNames', initLabels);

%% --- Output Matrix Preparation
Out_Transposed1 = [channelNums Vm maxfreq minI inres sfafirst spikefreq(1:11,:)']';
[~, idx] = sort(Out_Transposed1(1,:));
Out_Transposed = Out_Transposed1(:,idx);

%% --- Initial Spiking Output
sw = min(j, 11); % restrict to 11 steps max
Out_Spiking = [channelNums'; initialfreq(1:sw,:); initialfreq2(1:sw,:)];

%% --- Final Cleanup
clearvars -except Out_Transposed D t Out_Spiking  
