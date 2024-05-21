function [data]  = filterECG1kHz(data, options)

% [filtdata]  = filterECG1kHz(data, options);
% apply a zerophase IIR filter (coeffs designed with sptool)
% Sampling Frequency of Fs = 1000Hz.
%
% options can be 'lp', 'hp' or 'bp' for low, high or band pass.
% Default == 'bp'; band pass.
%
% Should work nicely for 16bit, 256Hz ECG. You may experience 
% You might like to try wavelets for a better filter 
%
% These routines are made available under the GNU general public license. 
% If you have not received a copy of this license, please download from 
% http://www.gnu.org/
%
% Please distribute (and modify) freely, commenting where you have 
% added modifications. The author would appreciate correspondence 
% regarding corrections, modifications, improvements etc.
%
% G. Clifford : gari@ieee.org

if nargin<2 
    options = 'bp';
end

if ( (strcmp(options,'lp')| strcmp(options,'hp')| strcmp(options,'bp') ) ==0)
 options = 'bp';
 fprintf('Using BP filter implementation\n')
end
    
% sampling frequency
Fs      = 1000;  % Hz    

% Data should be zeromeaned before calling this function
data = data-mean(data);

% set up transfer function
denominator_lp = [
   0.01000000000000
  -0.10031567744300
   0.45782014065187
  -1.25470427029339
   2.29432839868804
  -2.93911935390185
   2.69148006270018
  -1.76184961046748
   0.80791962849992
  -0.24716791918462
   0.04540208495209
  -0.00379348419655
];
denominator_lp = denominator_lp*100;

numerator_lp =[  
  0.00074008529810
  -0.00636191114604
   0.02384554037348
  -0.04970037731263
   0.05854575597887
  -0.02706909293097
  -0.02706909293097
   0.05854575597886
  -0.04970037731263
   0.02384554037348
  -0.00636191114604
   0.00074008529810
];


numerator_hp = [ 
     0.99077462632576  
    -0.99077462632576
];

denominator_hp = [
    1.00000000000000  
   -0.98154925265151
]; 


if ( (strcmp(options,'lp')| strcmp(options,'bp'))==1)
% low pass filter
if(exist('filtfilt'))
 data = filtfilt(numerator_lp,denominator_lp,data);
else 
 data = filter(numerator_lp,denominator_lp,data);
end
end



% Don't high pass filter it if you don't want to remove the baseline 
% fluctuations due to resp, BP? and electrode noise?
if ( (strcmp(options,'hp')| strcmp(options,'bp'))==1)
if(exist('filtfilt'))
 data = filtfilt(numerator_hp,denominator_hp,data);
else 
 data = filter(numerator_hp,denominator_hp,data);
end
end


% correct for amplitude distortion. !!!??????
%mean_s = mean(new_data);
%mean_a = mean(aff_hp);
%a = aff*(mean_s/mean_a);

