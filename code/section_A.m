%%
% File         : section_A.m   
% Author       : Mauricio Villarroel
% Created on   : May 2024
% Last updated : $Id$
% ________________________________________________________________________
%
%
% Copyright (C) 2024 Mauricio Villarroel. 
% All rights reserved.
%
% SPDX-License-Identifer:  GPL-2.0-only
%
% ________________________________________________________________________
%
%
% DESCRIPTON
% ----------
%
% Launch script for the ECG R peak detection practice session
%
% ________________________________________________________________________


%% loading data file

load s0016lre;

% Define the sampling frequency
samp_freq=1000;

% Choose a sample ECG lead, Lead I in this case
data = s0016lre.data(:,2); 

%% QRS detector

[~, ~, ~, R_index, ~, ~] = rpeakdetect(data, samp_freq);

% Plot RR peaks

figure;

% Contruct time vector for data
t = (0:length(data) -1)/samp_freq;

plot( t,   data,   'b' ); hold on;
plot( t(R_index), data(R_index), 'r*');
legend('ECG', 'R peaks');
title('ECG R peak detection');
ylabel('ECG (AD units)');
xlabel('Time (sec)');


%% Computing HR



%% Compute HRV metrics

