%% ADD CODEBASES TO PATH
addpath(genpath('~/code/pac'));
addpath(genpath('~/code/svndl/svndl2017'));

%% IMPORT DATA
clear all
close all

parentDir = '~/code/pac/Data'; % top-level directory with data and RCA figure subfolders
subj = 'DM'; % name of data folder

EEG = pop_fileio(sprintf('%s/%s/%s.vhdr',parentDir,subj,subj));
EEG.setname = 'EEG';
EEG = eeg_checkset(EEG);

%% define epochs
EEG = pop_epoch( EEG, {'S  2'}, [0 2], 'newname', 'EEG epochs', 'epochinfo', 'yes');
EEG = eeg_checkset(EEG);

data = EEG.data;
% what are the dimensions? 3 x 1000 x 100

%% TEST PLOTS

figure
for i = 1:100
    subplot(10,10,i)
    plot(data(1,:,i))
end

%% SAMPLE FFT
Fs = EEG.srate;
L = 1000;
sampdata = data(1,:,1);

spec = fft(sampdata);
% spec = fft(avgData(1,:));
P2 = abs(spec/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% avg over epochs?
avgData = squeeze(nanmean(data,3));