%% Section A: The Electrocardiogram
%% Exercise 1: Exploratory ECG Data Analysis
%% Load the PhysioNet data
LoadData = load('/Users/ellenhobson/Desktop/DATA/B18_ECG_data/PhysioNetData.mat');
%% Q1. Determine proportion of classes in the dataset from RecordingInfo
InfoValues = LoadData.RecordingInfo;
tabulate(InfoValues{:,2})

%% Q2. Plot ECG data
fs = 300;  %Sample frequency (Hz)
Normal = LoadData.Signals{16}; % Normal Rhythm

samples = 1:length(Normal);  % 'signal' is your vector of ECG values
% 'time' is a corresponding monotonic vector of time values
time = samples./fs;
range = [find(time==2):find(time==4.5)];
figure;
hold on 
plot(time(range), Normal(range))
xlabel('Time, s')
ylabel('Amplitude, mV')
title('ECG data, subject A04984')
% Label PQRST points
text(3.637,39,'P')
text(3.733,19,'Q')
text(3.763,1123,'R')
text(3.827,-59,'S')
text(4.01,164,'T')

%% Q3. Visualise 4 classes of ECG

Normal = LoadData.Signals{16};
AFib = LoadData.Signals{13};
Other = LoadData.Signals{12};
Noise = LoadData.Signals{2};
range = [find(time==2):find(time==6)];

% Plot above signal
figure;
hold on
subplot(4,1,1);
plot(time(range), Normal(range))
xlabel('Time, s')
ylabel('Amplitude, mV')
title('Normal - ECG data, subject A04984')

subplot(4,1,2);
plot(time(range), AFib(range))
xlabel('Time, s')
ylabel('Amplitude, mV')
title('Afib - ECG data, subject A02784')

subplot(4,1,3);
plot(time(range), Other(range))
xlabel('Time, s')
ylabel('Amplitude, mV')
title('Other - ECG data, subject A03519')

subplot(4,1,4);
plot(time(range), Noise(range))
xlabel('Time, s')
ylabel('Amplitude, mV')
title('Noisy - ECG data, subject A05233')

%% Exercise 2: Pre-Processing ECG
%% Q1. Detrend the data
Normal = LoadData.Signals{16};
AFib = LoadData.Signals{13};

Normal_mean = mean(Normal);
AFib_mean = mean(AFib);

Z_score_Normal = normalize(Normal - Normal_mean);
Z_score_AFib = normalize(AFib  - AFib_mean);

Normal = Z_score_Normal;
AFib = Z_score_AFib;

%% Q2. Apply band-pass filter

fs = 300;       % sampling freq Hz
f_low = 100;    % low-pass freq Hz
f_high = 0.5;   % high-pass freq Hz
norder = 4;     % filter order

% filter_ecg.m filters the original signal and outputs
% the band-pass filtered ECG signal

AFib_filtered_signal = filter_ecg(AFib, fs);
Normal_filtered_signal = filter_ecg(Normal, fs);

%% Q3.  Plot resulting recordings
samples_bpf = 1:length(AFib_filtered_signal);
time_bpf = samples_bpf./fs;

figure;
hold on
subplot(2,1,1);
plot(time_bpf, AFib_filtered_signal)
xlabel('Time, s')
ylabel('Amplitude, mV')
title('Afib - ECG data, subject A02784')


subplot(2,1,2);
plot(time_bpf, Normal_filtered_signal)
xlabel('Time, s')
ylabel('Amplitude, mV')
title('Normal - ECG data, subject A04984')

%% Exercise 3: Freq Represenation of ECG Signals

%% Q1. Plot modified periodogram PSD
%% To return the modified periodogram power spectral density (PSD) estimate, pxx, of the input signal, x, using a hamming window:

% Input ECG signal, measured in mV
    % x = Normal, AFib
% Note: ensure x is detrended for prior to input into the
% periodogram / FFT;

% f_low - the low-pass frequency used during filtering previously
% Hamming window function - AFib:
window = hamming(length(AFib)); 
% modified - periodogram:
[pxxAFib,fAFib] = periodogram(AFib,window,[],fs,'psd');

% Hamming window function - Normal:
window = hamming(length(Normal)); 
% modified -periodogram:
[pxxNormal,fNormal] = periodogram(Normal,window,[],fs,'psd');


% Plot the PSD
figure
subplot(2,1,1)
plot(fAFib, pxxAFib)
xlim([0, f_low])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (mV^{2}Hz^{-1})')
title('PSD - Afib, subject A02784')


subplot(2,1,2)
plot(fNormal, pxxNormal)
xlim([0, f_low])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (mV^{2}Hz^{-1})')
title('PSD - Normal, subject A04984')

% Change x boundary conditions to 30 as negligible value of y after this
figure
subplot(2,1,1)
plot(fAFib, pxxAFib)
xlim([0, 30])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (mV^{2}Hz^{-1})')
title('PSD - Afib, subject A02784')


subplot(2,1,2)
plot(fNormal, pxxNormal)
xlim([0, 30])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (mV^{2}Hz^{-1})')
title('PSD - Normal, subject A04984')

% hint: for decibels (dB) plot 10*log10(pxx) instead;


%% Q2. Plot Welch's PSD

% Welch's power spectral density estimate - AFib
% (pwelch.m automatically implements a Hamming window):
[pxx_pwelchAFib,f_pwelchAFib] = pwelch(AFib,[],[],[], fs,'psd');

% Welch's power spectral density estimate - Normal
[pxx_pwelchNormal,f_pwelchNormal] = pwelch(Normal,[],[],[], fs,'psd');

% plot Welch's PSD
figure
subplot(2,1,1)
plot(f_pwelchAFib,pxx_pwelchAFib)
xlim([0, 30])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (mV^{2}Hz^{-1})')
title('Welchs PSD - Afib, subject A02784')


subplot(2,1,2)
plot(f_pwelchNormal,pxx_pwelchNormal)
xlim([0, 30])
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (mV^{2}Hz^{-1})')
title('Welchs PSD - Normal, subject A04984')

%% Exercise 4: Determining Heart Rate and Heart Rate Variability
%% Part 1: Determining Heart Rate through R-Peak Detection
%% Q1. R-Peak detection using findpeaks.m


%Your chosen scalar threshold values, note MinPeakDistance is in samples not seconds
% AFib
MinPeakHeight = 2; MinPeakDistance = 200;
%findpeaks will output the R-peaks detected (R_pks) and the locations of the peaks (R_locs);
[R_pksAFib,R_locsAFib]=findpeaks(AFib, 'MinPeakHeight',MinPeakHeight,'MinPeakDistance', MinPeakDistance);

% Normal
MinPeakHeight = 1; MinPeakDistance = 200;
[R_pksNormal,R_locsNormal]=findpeaks(Normal, 'MinPeakHeight',MinPeakHeight,'MinPeakDistance', MinPeakDistance);

% Plot the R-peak annotations;
figure
subplot(2,1,1)
plot(time, AFib); 
hold on; 
plot(time(R_locsAFib), AFib(R_locsAFib), 'ro');
xlabel('Time, s')
ylabel('Amplitude, mV')
title('AFib R-Peaks, subject A02784')

subplot(2,1,2)
plot(time, Normal); 
hold on; 
plot(time(R_locsNormal), Normal(R_locsNormal), 'ro');
xlabel('Time, s')
ylabel('Amplitude, mV')
title('Normal R-Peaks, subject A04984')

%% Q2. Determine mean heart rate
% Mean HR, number of R-Peaks in 60 secs
% 60  / Average distance between peaks 
% Mean R-R peak difference
% AFib
Rpk_diffAFib = diff(R_locsAFib)./fs;
Rpk_diffAFib_mean = mean(Rpk_diffAFib);
HRAFib = 60/Rpk_diffAFib_mean;

% Normal
Rpk_diffNormal = diff(R_locsNormal)./fs;
Rpk_diffNormal_mean = mean(Rpk_diffNormal);
HRNormal = 60/Rpk_diffNormal_mean


%% Part 2: Determining Heart Rate Variability Measures
%% Q3. Plot R-R interval times over the duration of the recording

%from the first interval (i.e. the second R-peak detected);
time_AFib = time(R_locsAFib(2:end)); 
R_R_AFib = Rpk_diffAFib*1000     ;   % 1000 as ms
R_R_Normal = Rpk_diffNormal*1000;

figure
subplot(2,1,1)
plot(time(R_locsAFib(2:end)), (R_R_AFib));  
hold on
plot(time(R_locsAFib(2:end)), (R_R_AFib), 'ro');
xlabel('Time, s')
ylabel('R-R Interval, ms')
title('AFib R-R Intervals, subject A02784')

subplot(2,1,2)
plot(time(R_locsNormal(2:end)), (R_R_Normal));  %1000 as ms
hold on
plot(time(R_locsNormal(2:end)), (R_R_Normal), 'ro');
xlabel('Time, s')
ylabel('R-R Interval, ms')
title('Normal R-R Intervals, subject A04984')

%% Q4. Plot Histogram of the R-R intervals
%% Constructing a simple probabilistic histogram 
BinWidth=20; %the width of the bins

figure
subplot(2,1,1)
histogram(R_R_AFib, 'BinWidth', BinWidth)
xlabel('R-R Interval, ms') 
ylabel('Count')
title('Histogram of the R-R intervals - AFib')

subplot(2,1,2)
histogram(R_R_Normal, 'BinWidth', BinWidth)
xlabel('R-R Interval, ms') 
ylabel('Count')
title('Histogram of the R-R intervals - Normal')

%% Q5. Plot Lomb-Scargle periodogram
%Lomb-Scargle PSD RR-interval Spectrum;

R_R_AFib = R_R_AFib - mean(R_R_AFib); %detrend the RR intervals (best practice);

%Determine the time points in the ecg each R-peak occurred
%from the first interval (i.e. the second R-peak detected);
t_AFib=time(R_locsAFib(2:end));

%Denote the frequencies we are interest in evaluating:
%[0.001, 0.6] Hz;
f_interest =0.001:0.001:0.6;

% Compute the Lomb-Scargle periodogram;
[pxx_plombAFib,f_plombAFib] = plomb((R_R_AFib), (t_AFib),(f_interest), 'psd');


R_R_Normal = R_R_Normal- mean(R_R_Normal); %detrend the RR intervals (best practice);

%Determine the time points in the ecg each R-peak occurred
%from the first interval (i.e. the second R-peak detected);
t_Normal=time(R_locsNormal(2:end));


% Compute the Lomb-Scargle periodogram;
[pxx_plombNormal,f_plombNormal] = plomb(R_R_Normal, t_Normal, f_interest, 'psd');

% Plot the PSD
figure
subplot(2,1,1)
plot(f_plombAFib , pxx_plombAFib);
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (s^{2}Hz^{-1})')
title('AFib - Lomb-Scargle periodogram PSD')

subplot(2,1,2)
plot(f_plombNormal , pxx_plombNormal);
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density (s^{2}Hz^{-1})')
title('Normal - Lomb-Scargle periodogram PSD')