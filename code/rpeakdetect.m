%%
% File         : rpeakdetect_deb.m         
% Author       : Gari Clifford, Mauricio Villarroel
% Created on   : 1995
% Last updated : $Id$
% ________________________________________________________________________
%
% Written by G. Clifford gari@ieee.org and made available under the
% GNU general public license. If you have not received a copy of this
% license, please download a copy from http://www.gnu.org/
%
% Please distribute (and modify) freely, commenting
% where you have added modifications.
% The author would appreciate correspondence regarding
% corrections, modifications, improvements etc.
%
% gari@ieee.org
%
% ________________________________________________________________________
%
% USAGE
% -----
%      [hrv, R_t, R_amp, R_index, S_t, S_amp, S_index]  = ...
%                   rpeakdetect(data, samp_freq, thresh, testmode)
%
% DESCRIPTON
% ----------
%
% WARNING: rpeakdetec version for DEBUGGING, contains errors
%
%     A batch QRS detector based upon that of Pan, Hamilton and Tompkins:
%
% * J. Pan & W. Tompkins - A real-time QRS detection algorithm
%   IEEE Transactions on Biomedical Engineering, vol. BME-32 NO. 3. 1985.
%
% * P. Hamilton & W. Tompkins. Quantitative Investigation of QRS
%   Detection  Rules Using the MIT/BIH Arrythmia Database.
%   IEEE Transactions on Biomedical Engineering, vol. BME-33, NO. 12.1986.
%
% Similar results reported by the authors above were achieved, without
% having to tune the thresholds on the MIT DB. An online version in C
% has also been written.
%
%
% INPUT
% -----
%
%     data          Original ECG waveform and possible time
%                   (assumed to be 1st column or row).
%
%     samp_freq     Sample frequency of "data" (samp_freq = 256Hz by
%                   default)
%
%     thresh        The 'triggering' threshold 'thresh' for the peaks
%                   in the 'integrated' waveform is 0.2 by default.
%
%     testmode      testmode = 0 (default) indicates no graphics
%                   diagnostics. Otherwise, you get to scan through
%                   each segment.
%
% OUTPUT
% ------
%
%     hrv
%
%     R_t == RR points in time
%
%     R_amp == amplitude of R peak in bpf data
%
%     S_amp == amplitude of following minmum.
%
% ________________________________________________________________________

function [hrv, R_t, R_amp, R_index, S_t, S_amp, S_index]  = ...
    rpeakdetect(data, samp_freq, thresh, testmode )

error(nargchk(2, inf, nargin, 'struct'));

%% Process function parameters with default values
%

% make threshold default 0.2 -> this was 0.15 on MIT data
if nargin < 4
    testmode = 0;
end
% make threshold default 0.2 -> this was 0.15 on MIT data
if nargin < 3
    thresh = 0.2;
end

%% Global configurable parameters
%

maxval=[];
minval=[];

% define a refactory period in seconds
refact = 0.2;

% and a minimum vetricular activation time (R peak to S wave)
min_vat = 0.01;

% Length of the data set
len = 0;

%% check format of data
%

[a b] = size(data);

% Check if data is in columns or rows
% and convert it into rows.

if(a>b)
 len =a;
end
if(b>a)
 len =b;
end

x = data(1, :);

% if there's no time axis - make one

if (a | b == 1);
    tt = 1/samp_freq:1/samp_freq:ceil(len/samp_freq);
    t = tt(1:len);
end

%% bandpass filter data
%

% remove mean
x = x-mean(x);

% FIR filtering stage

bpf=x; %Initialise

if( (samp_freq==128) && (exist('filterECG128Hz', 'file')~=0) )
    bpf = filterECG128Hz(x);
    
elseif ( (samp_freq==256) && (exist('filterECG256Hz', 'file')~=0) )
    bpf = filterECG256Hz(x);
    
elseif ( (samp_freq==1000) && (exist('filterECG1kHz', 'file')~=0) )
    if (testmode ~=0)
        fprintf('apply 1kHz zerophase IIR filter.\n');
    end
    bpf = filterECG1kHz(x,'bp');
end

%% Differentiate data
%

dff=[];
for i=1:length(bpf)-1
    dff=[dff; bpf(i+1)-bpf(i)];
end

%% square data
%

sqr = dff*dff;

len = len-1; % how long is the new vector now?

%% Integrate data over window 'd'
%

Nint=7;
sqr_padded=[zeros(Nint-1,1); sqr];
for i=7:length(sqr_padded)
    sqrint(i-6)=sum(sqr_padded(i-6:i));
end

% integrate
mdfint = medfilt1(sqrint,10)';

% remove filter delay for scanning back through ECG
delay = ceil(length(d)/2);
mdfint = mdfint(delay:length(mdfint));


%% Segment search area
%

% first find the highest bumps in the data
max_h = max (mdfint(round(len/4):round(3*len/4)));

% then build an array of segments to look in
% thresh = 0.2;
poss_reg = mdfint>(thresh*max_h);


%% Find peaks
%

% find indices into boudaries of each segment

left=[];
right=[];
for i=1:length(poss_reg)
    if (poss_reg(i)==0) && (poss_reg(i+1)==1)
        left=[left; i+1];
    elseif (poss_reg(i)==1) && (poss_reg(i+1)==0)
        right=[right; i+1];
    end
end


% intialise
maxloc=1;
minloc=1;
R_amp=0;
S_amp=0;

% loop through all possibilities
for(i=1:length(left))
    [maxval(i) maxloc(i)] = max( bpf(left(i):right(i)) );
    [minval(i) minloc(i)] = min( bpf(left(i):right(i)) );
    maxloc(i) = maxloc(i)-1+left(i); % add offset of present location
    minloc(i) = minloc(i)-1+left(i); % add offset of present location
end

R_index = maxloc;
R_t   = t(maxloc);
R_amp = maxval;

S_index = minloc;
S_t   = t(minloc);

% Assuming the S-wave is the lowest
% amp in the given window

S_amp = minval;   


%% Check for lead inversion
%

% i.e. do minima precede maxima?
% if (minloc(length(minloc))<maxloc(length(minloc)))
%    R_t   = t(minloc);
%    R_amp = minval;
%    S_t   = t(maxloc);
%    S_amp = maxval;
% end


% Remove artifacts

%%% due to too short a Ventricular Activation Time (VAT)
vat = abs(S_t-R_t);
min_vat=0.03;

good_qrs_index = find(vat>min_vat);

len_good_qrs_index = length(good_qrs_index);

if(len_good_qrs_index > 0)
    R_t = R_t(good_qrs_index);
    R_amp = R_amp(good_qrs_index);
    R_index = R_index(good_qrs_index);
    S_t   = S_t(good_qrs_index);
    S_amp = S_amp(good_qrs_index);
    S_index = S_index(good_qrs_index);
end

%%% due to too short a refactory time

% add a dummy R peak, to prevent "diff" removing the last beat.
R_t(length(R_t) + 1) = R_t(length(R_t)) + 1;
hrv  = diff(R_t);
sinus_beat_index = find(hrv>refact);

len_sinus_beat_index = length(sinus_beat_index);

if(len_sinus_beat_index > 0)
    R_t = R_t(sinus_beat_index);
    R_amp = R_amp(sinus_beat_index);
    R_index = R_index(sinus_beat_index);
    S_t   = S_t(sinus_beat_index);
    S_amp = S_amp(sinus_beat_index);
    S_index = S_index(sinus_beat_index);
end

% Question: What for?
hrv  = diff(R_t);
resp = R_amp-S_amp;

% Find the local max for every r peak point in the real ecg stream.

%window to look a maximum for every r peak
dt = ceil(samp_freq * 0.015);

for i=1:length(R_index)
    
    rpoint = R_index(i);
    
    if( (rpoint-dt) > 0)
        tmpSegR = data(rpoint-dt:rpoint+dt);
    else
        tmpSegR = data(1:rpoint+dt);
    end
    %h2 = figure;
    %plot(tmpSegR);hold on;
    
    if(rpoint >= 0)
        maxp=max(tmpSegR);
    else
        maxp=min(tmpSegR);
    end
    
    loc= find(tmpSegR == maxp);
    
    if(length(loc) > 1)
        loc=floor(mean(loc));
    end
    
    %plot(loc, tmpSegR(loc), 'r*');
    R_index(i)= rpoint - dt + loc - 1;
    %pause
    %close(h2);
end

%% Plot results
%

if (testmode~=0)

    figure;
    
    axh(1) = subplot(4,1,1);
    plot(t, x,   'b'); hold on;
    plot(t, bpf, 'r');
    legend('Raw ECG', 'Zero-pahse FIR filtered ECG (red)');
    title('Raw ECG');
    ylabel('ECG');
    
    axh(2) = subplot(4,1,2);
    plot( t(1:length(mdfint)), mdfint, 'b'); hold on;
    %plot(t(1:length(sqr)),sqr);hold on;
    plot( t,        max(mdfint)*bpf/(2*max(bpf)), 'r'  );
    plot( t(left),  mdfint(left),                 'og' );
    plot( t(right), mdfint(right),                'om' );
    title('Integrated data with scan boundaries over scaled ECG');
    ylabel('Int ECG')
    
    axh(3) = subplot(4,1,3);
    plot( t,   bpf,   'r' ); hold on;
    plot( R_t, R_amp, '+k');
    plot( S_t, S_amp, '+g');
    legend('ECG', 'R peaks', 'S peaks');
    title('ECG peak detection');
    ylabel('ECG+R+S');
    
    axh(4) = subplot(4,1,4);
    plot(R_t(1:length(hrv)), hrv, 'r+' );
    title('RR intervals')
    ylabel('RR (s)');
    xlabel('Time (sec)');
    
    linkaxes(axh, 'x');
    
end
