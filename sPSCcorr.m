%% sPSCcorr Analysis of spontaneous event correlation in multicellular patch-clamp data
% This script detects spontaneous postsynaptic current (sPSC) events in multichannel
% patch-clamp recording data using a template-based approach (based on Jonas et al., 1993). Event
% detection is performed using 'minidet', developed by Alois Schlögl, IST
% Austria and available as part of the Biosig package: https://biosig.sourceforge.net/
%
%     sPSCcorr.m Copyright (C) 2025 Jake Watson
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%     You should have received a copy of the GNU General Public LicenseAdd commentMore actions
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% This code has been written for analysis of recording files (CFS file type) with specific
% structure, and is not intended as a general analysis script. I hope that
% you find it helpful, but please ensure all parameters and analyses are
% performed correctly and appropriately for your own data structure.
%
% Correlation analysis is made currently for 3 cell types. Cells must be
% assigned to either class A, B, or C. These codes can be changed for
% relevance to your cell types, but this must be thoroughly checked
% throughout the script.
%
% This code requires Biosig (https://biosig.sourceforge.net/) and the NaN toolbox.
%
% User input is required to input the location and identity of the
% recording file for analysis, and the identity of each recorded cell, for
% within and between class correlation analysis. This information is input
% at the beginning of the script, and further user input is not required.
% Some parameters may need adjustment for best analysis, for example the
% template detection parameters of minidet.
%
% There is the option to remove sections of data from analysis through
% manual annotation. For example, this allows manual curation of data to
% remove rare artefacts, unwanted events, or regions of unstable recording.
% Such regions should be manually curated in Sigviewer
% (https://sigviewer.sourceforge.net/index.html), marked as events with 
% Event Type 6, and event file (EVT) saved with the same
% filename as recording data, with '_artwipe.evt' appended. This allows
% correspondence between recording and artefact event file. Please note,
% should this functionality be used, selecting regions on only one
% recording channel will remove data only from that recording channel, and
% therefore may influence the detection of correlated events during this
% time period with other channels. This caveat may affect correlation
% conclusions and event frequencies unless time windows for removal are marked on all recorded
% channels (recommended).
%
% Results are output into a structure 'RESULTS', where:
%
% RESULTS.EvtRate - reports the frequency of detected events during the
% recording, excluding the duration of test pulse, and artefact selected
% areas. Units = Hz.
%
% RESULTS.CorrelationArray - outputs a matrix reporting the percentage of
% correlated events for each cell pair. As the event rate will differ
% between each cell the percentage is calculated in both directions,
% reporting the percentage of events of Cell 1 correlated with Cell 2, and
% the percentage of events of Cell 2 correlated with Cell 1. Data is
% structured as rows = Cell for which the percentage of its events
% correlated with each cell are noted in each column. e.g. the percentage
% of events of Cell 1 correlated with Cell 2 would be found in (1,2).
% Note, the diagonal will not contain analysed data due to irrelevance of
% self correlations. Units = %.
%
% RESULTS.EvtList - outputs a list of all detected correlated events.
% Column 1 includes event time, column 2 the reference channel, and column
% 3 the channel in which the event was detected.
%
% RESULTS.SIMMED - outputs a matrix equivalent to RESULTS.CorrelationArray
% but for randomly shuffled event times.
%
% RESULTS.PSCFreq_ABC - outputs the results of RESULTS.EvtRate sorted into
% columns based on cell type, where column 1 is a list of event frequencies
% for cell type A, columns 2 for B, and column 3 for C. Units = Hz.
%
% RESULTS.FracMulti - outputs a matrix displaying the fraction of
% correlated events shared across more than just pairs of channels. This
% data is structured as one row per channel, with column 1 indicating the
% recording channel, column 2 the fraction of correlated events coincident
% with just 1 other channel, column 3 the fraction with two other channels,
% and 4 the fraction coincident with 3 or more other channels.
%
% Written by Jake Watson (2025), based on code by Alois Schlögl, for Watson et al. Cell Reports 2025.
%
%% Input Data Information

DATADIR   = 'X:\XXX\'; % Change to directory with recording files
DATAFILES = {'XXX.cfs'}; % Change to the name of the recording file
CellID    = ['A', 'A', 'B', 'B', 'C'];  % Input code for each cell type

rsPulseAmp = -0.002; % size of test pulse in Volts
t1 = 1.4; % start..
t2 = 1.48; % .. and end time of baseline region for series resistance calculation (seconds)
t3 = 1.48; % start..
t4 = 1.52; % .. and end time of test pulse negative going search region (seconds)
t5 = 1.68; % time after test pulse for artefact removal (seconds)

Rmthresh = 20; % Insert the threshold for median series resistance during recording
% for inclusion in further analysis. Note that it is important to manually
% assess the array of series resistance values across long recordings to
% determine the degree of change in this parameter, which can seriously
% influence the conclusions drawn. For further information on the impact of
% event detection on sPSC analysis, please refer to Greger & Watson,
% bioRxiv 2024 (https://doi.org/10.1101/2024.10.26.620084).

span = 0.001; % Set the window for consideration of a 'correlated event' here.
% The default value is 0.001 s. This corresponds to events within a window of 0.5 ms
% before or after an event for consideration as 'correlated timing'.

%% Event detection
[data, HDR] = mexSLOAD(fullfile(DATADIR, DATAFILES{1}));  % Load data
[~, f, e]   = fileparts(DATAFILES{1});
HDR.FILE.Name = f; HDR.FILE.Ext = e;

[~, yScale] = physicalunits(HDR.PhysDimCode);
Fs = round(HDR.SampleRate);                    % Extract sampling rate of recording from Header

% 'selpos' contains the location of sweep transition marks extracted from CFS
% file header, allowing reshaping to sweep data structure if required. A sweep
% transition mark is encoded at the beginning and end of recording,
% therefore the number of sweeps is the size of selpos - 1.
selpos = [1; HDR.EVENT.POS(HDR.EVENT.TYP == hex2dec('7ffe')); size(data,1) + 1];
channels = size(data, 2);
Sweeps = length(selpos) - 1;
ResArray = NaN(Sweeps, channels);
D = data;  % Working data matrix to avoid altering imported data. 'data' remains unaltered during this script, while D can be manipulated.

% Compute series resistance for exclusion of cells, and remove test
% pulses from recording
for k2 = 1:Sweeps
    x = data(selpos(k2):selpos(k2+1)-1, :);

    % Define baseline and pulse windows
    v2 = median(x(round(t1*Fs):round(t2*Fs), :));
    minC = min(x(round(t3*Fs):round(t4*Fs), :));

    ResArray(k2, :) = (rsPulseAmp ./ (minC - v2)) * 1e6;  % Resistance in MOhms

    % Remove test pulse segment from analysis
    D(selpos(k2) + round(t3*Fs):selpos(k2) + round(t5*Fs), :) = NaN;
end

dur1 = sum(~isnan(D(:,1))) / Fs;  % Total duration of event search area

%% Optional Artefact Removal.
% Artefacts can be wiped from the data to prevent detection and
% inclusion in final datasets. Artefacts must be marked using Sigviewer
% and saved as a separate EVT file under type '6'. Please see script
% header for full details.

if exist(fullfile(DATADIR, [HDR.FILE.Name, '_artwipe.evt']), 'file')
    [~, EVT1] = mexSLOAD(fullfile(DATADIR, [HDR.FILE.Name, '_artwipe.evt']));
    delloc = EVT1.EVENT.POS(EVT1.EVENT.TYP == 6);
    deldur = EVT1.EVENT.DUR(EVT1.EVENT.TYP == 6);

    for s = 1:length(delloc)
        D(delloc(s):(delloc(s) + deldur(s)), :) = NaN;
    end
end

dur2 = sum(~isnan(D(:,1))) / Fs; % Updated duration of event search area after artefact wiping.

%% Event Detection using minidet.
posns = [];
chancode = [];

for i = 1:channels
    lpdata = D(:, i);  % Channel data

    % Event detection using template search (minidet)
    refract = 0.004 * Fs;  % Refractory period is default set at 0.004 s.
    rCrit = 0.5; aCrit = 10; % Template search parameters. It is recommended
    % that these parameters be adjusted for optimal detection of recorded events.
    % Default values are 'rCrit = 0.5; aCrit = 10'
    res = minidet(lpdata, [], rCrit, aCrit, refract);
    pos1 = res.tEventList;


    if isempty(pos1), continue; end

    % The output of minidet is converted to two vectors. 'posns(x)' notes
    % the event time of each detected event. 'chancode(x)' contains the
    % corresponding channel for each event. Note that event times are
    % output from minidet as the midpoint of event rise.
    posns = [posns; pos1(:)];
    chancode = [chancode; repmat(i, length(pos1), 1)];
end

%% Save Events to File.
% Events are saved to a separate EVT file, and can be viewed for
% accuracy of detection in Sigviewer, by importing both recorded data
% and the detected event EVT.

EVT.TYPE = 'GDF'; EVT.VERSION = 3.0;
EVT.FileName = fullfile(HDR.FILE.Path, [HDR.FILE.Name, '.evtdet.evt']);
EVT.NS = 0;
EVT.EVENT.SampleRate = Fs;
EVT.EVENT.POS = [posns; selpos(Sweeps)+round(19.5*Fs)-1];
EVT.EVENT.TYP = [repmat(hex2dec('0214'), length(posns), 1); hex2dec('0213')];
EVT.EVENT.CHN = [chancode; 0];
EVT.EVENT.DUR = [zeros(length(posns), 1); round(0.1*Fs)];
EVT.CodeDesc{1} = 'RS pulse';

EVT = sopen(EVT, 'w');
EVT = sclose(EVT);


%% Preparation for event analysis.

clearvars -except ResArray posns chancode channels Fs selpos Sweeps data CellID dur1 dur2 D

%% Calculate Event Rates and Channel-wise Metrics

AvRm = median(ResArray, 1);  % Median series resistance per channel
keepidx = find(AvRm < Rmthresh);   % Thresholded for poor series resistance

RESULTS = struct();
RESULTS.EvtRate = zeros(1, channels);

for i = 1:channels
    fieldset = sprintf('Set_%d', i);
    tidx = find(chancode == i);
    RESULTS.(fieldset) = posns(tidx) / Fs;
    RESULTS.EvtRate(i) = length(tidx) / dur2;
end

%% Cross-channel Correlation Calculation

setlist = [nchoosek(1:channels,2); flip(nchoosek(1:channels,2),2)];
RESULTS.CorrelationArray = NaN(channels);
RESULTS.EvtList = zeros(0, 3);

range = span / 2;

for i = 1:size(setlist, 1)
    pre = setlist(i,1); post = setlist(i,2);
    t_pre = RESULTS.(sprintf('Set_%d', pre));
    t_post = RESULTS.(sprintf('Set_%d', post));

    matches = [];
    for j = 1:length(t_pre)
        match = find((t_post > t_pre(j)-range) & (t_post < t_pre(j)+range));
        matches = [matches; match];
    end

    evts = t_post(matches);
    if isempty(evts)
        RESULTS.CorrelationArray(pre, post) = 0;
    else
        RESULTS.CorrelationArray(pre, post) = length(evts) / length(t_pre);
        RESULTS.EvtList = [RESULTS.EvtList; [evts, repmat(pre, length(evts), 1), repmat(post, length(evts), 1)]];
    end
end

%% Simulated Random Correlation
% The script continues by calculating the expected percentage correlation
% between recorded cells, should events occur at the recorded frequency,
% but with random event timing.

RESULTS.SIMMED = NaN(channels);

for i = 1:size(setlist,1)
    pre = setlist(i,1); post = setlist(i,2);
    sz1 = length(RESULTS.(sprintf('Set_%d', pre)));
    sz2 = length(RESULTS.(sprintf('Set_%d', post)));

    maxtime_s = dur2;
    simrange = range * Fs;
    runs = 100; % Change this value to alter the number of simulation runs performed per cell pair. Default = 100.
    Tempres = [];

    for k = 1:runs
        rand1 = sort(randi(round(maxtime_s * Fs), sz1, 1));
        rand2 = sort(randi(round(maxtime_s * Fs), sz2, 1));
        match = [];

        for j = 1:sz1
            m = find((rand2 > rand1(j)-simrange) & (rand2 < rand1(j)+simrange));
            match = [match; m];
        end

        Tempres(k) = length(match);
    end

    RESULTS.SIMMED(pre, post) = median(Tempres) / sz1; % The median percentage of all simulation runs per cell pair is reported as the output.
end

%% Categorize by Cell Type (A/B/C) for Grouped Summary

AvRm(isinf(AvRm)) = 1000; % Filler data for handling infinite numbers which can occur from dead channels etc.

Aloc = find(CellID == 'A' & AvRm <= Rmthresh);
Bloc = find(CellID == 'B' & AvRm <= Rmthresh);
Cloc = find(CellID == 'C' & AvRm <= Rmthresh);

RESULTS.PSCFreq_ABC = NaN(max([length(Aloc), length(Bloc), length(Cloc)]), 3);
RESULTS.PSCFreq_ABC(1:length(Aloc),1) = RESULTS.EvtRate(Aloc)';
RESULTS.PSCFreq_ABC(1:length(Bloc),2) = RESULTS.EvtRate(Bloc)';
RESULTS.PSCFreq_ABC(1:length(Cloc),3) = RESULTS.EvtRate(Cloc)';

%% Multi-channel Coincidence Analysis
% The following section detects events that are correlated across more than
% two channels. This allows for detection of multi-cell 'network' events
% which are unlikely to result from monosynaptic input to pairs of neurons.

[~, ~, idxc] = unique(RESULTS.EvtList(:,1));
out = unique(RESULTS.EvtList(histcounts(idxc, numel(idxc)) > 1,1));

RESULTS.FracMulti = NaN(channels, 4);  % Channel | Frac1 | Frac2 | Frac3+
multilist = [];

for i = 1:channels
    idx = RESULTS.EvtList(:,3) == i;
    [freqs, chans] = groupcounts(RESULTS.EvtList(idx,1));
    if isempty(freqs)
        RESULTS.FracMulti(i,2:4) = [0 0 0];
    else
        RESULTS.FracMulti(i,2:4) = [sum(freqs == 1), sum(freqs == 2), sum(freqs >= 3)] ./ length(freqs);
    end
    multilist = [multilist; chans(freqs > 1)];
end
