%% Project 2 Revisited
clear;clc;clear all;close all;

% load message
file = load('Message.mat');
y = file.data;
Fs = file.fs;
n = length(y);
nn = 1:n;

% DFT and frequency range
dft = fft(y);
f1 = (0:n-1)*(Fs/n); 
f1_range = find(f1>100 & f1<700);

% Standard lowpass filter

denoised_y = lowpass(y,0.75);

dft_mod = zeros(size(dft));
dft_mod(f1_range,:) = dft(f1_range,:);
dft_mod1 = dft_mod;
signal = real(ifft(dft_mod));

% Gaussian Filter
cconv = @(a,b) real(ifft(fft(a).*fft(b)));

normalize = @(h) h/sum(h(:)); 
t = [0:n/2-1, -n/2:-1]';
h = @(mu) normalize(exp(-(t.^2.25)/(2*mu^2)));
mu = 15;

sigma = 1;
denoise = @(x,mu) cconv(h(mu), x);
denoised = denoise(y,mu);

[Q,R] = qr(denoised,0);
[U,S,V] = svd(R);

yy1 = y*pinv(V)*pinv(U')*S;
yy1 = yy1./max(yy1);
yy1(:,1) = yy1(:,2);

% FFT Filter
% sigfilt = fftfilt(ones(1,10)/25,y); % signal
% sigfilt_1 = wdenoise(sigfilt,'DenoisingMethod','Minimax');
% dftfilt_2 = fft(sigfilt_1);
% dftfilt_2mod = zeros(size(dftfilt_2));
% dftfilt_2mod(f1_range,:) = dftfilt_2(f1_range,:);
% dftfilt_2mod1 = dftfilt_2mod(:,1);
% dftfilt_2mod2 = dftfilt_2mod(:,2);

% Signals 
% signal_1 = imag(ifft(dftfilt_2mod));
% signal_1mod = real(ifft(dftfilt_2modmod));
% signal1 = real(ifft(dft_mod1));
% signal2 = real(ifft(dft_mean));

% Last hunnid thou of signal
% yy_voice = y(722640:752640)';
% yy_noise = y(652640:722640)';
% dft_yy_voice = fft(yy_voice);
% dft_yy_noise = fft(yy_noise);
% dft_yy_voice_mod = zeros(size(dft_yy_voice));
% dft_yy_voice_mod(f1_range,:) = dft_yy_voice(f1_range,:);

% for ii = 1:1708
%     if dft_yy_voice_mod(ii) >= 500
%         dft_yy_voice_mod(ii) = ;
%     end
% end
% 
% for ii = 1706:17068
%     if dft_mod(ii) >= dft_yy1_noise_mod(ii)
%         dft_mod(ii) = 0.25*dft_mod(ii);
%     end
% end

% Convolution of the Signal

% sig_conv = conv2(y,y,'same');
% dft_conv = fft(sig_conv);
% 
% max_conv = max(sig_conv);
% max_sig1 = max(signal1);
% factor = max_sig1./max_conv;
% scale_sig_conv = factor.*sig_conv;
