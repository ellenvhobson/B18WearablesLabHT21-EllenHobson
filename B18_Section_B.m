%% Section B: The Electroencephalogram
%% Execise 1: Investigating EEG Data
%% Q1. Load EEG Data
%% Q1. Load EEG Data & Plot entire EEG

LoadOpenEyes = load('/Users/ellenhobson/Desktop/DATA/B18_EEG_data/EEGeyesopen.mat');
LoadClosedEyes  = load('/Users/ellenhobson/Desktop/DATA/B18_EEG_data/EEGeyesclosed.mat');

% Trial 1
EyesOpen = LoadOpenEyes.eyesopen;
% Trial 2
EyesClosed = LoadClosedEyes.eyesclosed;

% Plot EEG samples, first 512 samples
sample = 1:512;
fs = 256;       % Sampling Freq Hz
t = sample./fs;

figure
subplot(2,1,1)
plot(t, EyesOpen(sample))
xlabel('Time, s')
ylabel('Amplitude, microV')
title('Open Eyes - EEG Signal')

subplot(2,1,2)
plot(t, EyesClosed(sample))
xlabel('Time, s')
ylabel('Amplitude, microV')
title('Closed Eyes - EEG Signal')

%% Q2. Pre-process the data by low-pass filtering parameters by reporting and analysing various parameter combinations

% Detrending & Normalising
EyesOpen_mean = mean(EyesOpen);
EyesClosed_mean = mean(EyesClosed);

Z_score_EyesOpen = normalize(EyesOpen - EyesOpen_mean);
Z_score_EyesClosed = normalize(EyesClosed  - EyesClosed_mean);

EyesOpen = Z_score_EyesOpen;
EyesClosed = Z_score_EyesClosed;


% Design a LPF
fc = 30; %cut−off frequcny in Hertz (Hz) 
norder = 4; %filter order
%The cutoff frequency Wn must be 0.0 < Wn < 1.0; 
Wn= fc/fs ;
%Design the filter
[B,A] = butter(norder, Wn);
%apply the LPF
EyesOpen_filtered = filter(B,A, EyesOpen);
EyesClosed_filtered = filter(B,A, EyesClosed);

% Plot first 512 samples after pre-processing
figure
subplot(2,1,1)
plot(t, EyesOpen_filtered(sample))
xlabel('Time, s')
ylabel('Amplitude, microV')
title('Open Eyes - Pre-Processed Filtered EEG Signal')

subplot(2,1,2)
plot(t, EyesClosed_filtered(sample))
xlabel('Time, s')
ylabel('Amplitude, microV')
title('Closed Eyes - Pre-Processed Filtered EEG Signal')

%% Q3. Plot modified periodogram PSD, upto 60Hz

% Hamming window function - EyesOpen:
window = hamming(length(EyesOpen_filtered)); 
% modified - periodogram:
[pxxEyesOpen_filtered,fEyesOpen_filtered] = periodogram(EyesOpen_filtered,window,[],fs,'psd');

% Hamming window function - EyesOpen:
window = hamming(length(EyesClosed_filtered)); 
[pxxEyesClosed_filtered,fEyesClosed_filtered] = periodogram(EyesClosed_filtered,window,[],fs,'psd');

figure
subplot(2,1,1)
plot(fEyesOpen_filtered, pxxEyesOpen_filtered)
xlim([0, 60])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (mV^{2}Hz^{-1})')
title('PSD - Eyes Open')


subplot(2,1,2)
plot(fEyesClosed_filtered, pxxEyesClosed_filtered)
xlim([0, 60])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (mV^{2}Hz^{-1})')
title('PSD - Eyes Closed')

%% Exercise 2: Quantify Sleep through EEG
%% Q1. Load Sleep Data
Sleep = load('/Users/ellenhobson/Desktop/DATA/B18_EEG_data/EEGSleepStateData.mat');
SleepData = Sleep.EEGSleepStateData;
%% Q2. Plot hnogram of 4-stage sleep categories
sleep_category_num = SleepData.sleep_category_num;

%determine a time vector: as sleep states are assigned 
%by a clinician in one−minute intervals :
time = 1:length(sleep_category_num);%plot the sleep hypnogram as a stairstep plot;

figure
stairs (time ,sleep_category_num)
%modify the ylimits as appropriate clearer visualisation;

ylim([0.5, 4.5]) 
xlabel('Time, mins') 
ylabel ('Sleep State')
title('Hynogram of the 4-stage sleep categories')

%% Q3. Plot modified periodogram PSD
% Detrend data
sleep_category_num = sleep_category_num - mean(sleep_category_num);

% Hamming window functionn:
window = hamming(length(sleep_category_num));
% cycles per hour, x60

[pxxsleep_category_num,fsleep_category_num] = periodogram(sleep_category_num,window,[],60,'psd');

% cycles per hour, x60
figure
subplot(2,1,1)
plot(fsleep_category_num, pxxsleep_category_num)

xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (mV^{2}Hz^{-1})')
title('PSD - 4-stage sleep categories')

