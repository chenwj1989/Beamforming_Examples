%DAS_BF
clear;
close all;
load('Computed_RIRs.mat');

nChannels = size(m_pos,1);
dMic = m_pos(2, 2) - m_pos(1, 2);
c = 340;
delta = pdist2(s_pos(1,:), m_pos(2,:)) - pdist2(s_pos(1,:), m_pos(1,:));
DOA_est = acosd(delta/dMic);

%==================Generate array signals==================================%
speechfilename = {'6319-275224-0008.flac', '6319-275224-0011.flac'};
noisefilename = {'noise1.wav', 'noise2.wav'};

[source1, fs] = audioread(speechfilename{1});
[source2, ~] = audioread(noisefilename{1});
[noise, fs_n] = audioread(noisefilename{1});
noise = resample(noise,fs,fs_n);

n_f = fs * 10; %2 seoconds
source1 = source1(1:n_f);
source2 = source2(1:n_f) ./ 3;
noise = noise(1:n_f);

rir = RIR_sources(:,:,1);
speech1 = fftfilt(rir, source1).*30;

rir = RIR_sources(:,:,2);
speech2 = fftfilt(rir,source2).*30;

arraySignal = speech1 + repmat(noise, 1, 5);
%==========================================================================%

%=============================plot Beamformer==============================%
%h_mvdr = mvdr(d, Phi_u);
theta = linspace(0, 2*pi, 200); % incidence angle range
%plot
figure(1);
f = [10000];
hw = dsb(nChannels, dMic, DOA_est, f);        % delay-and-sum filter
p = plotBeamformer(nChannels, f, dMic, hw, theta);

% f = linspace(0, fs, 256);% freqeuncy vector
% hw = dsb(nChannels, dMic, DOA_est, f);        % delay-and-sum filter
% p = plotBeamformer(nChannels, f, dMic, hw, theta);
%==========================================================================%

%=============================apply Beamformer=============================%
x_das = applyBeamforming(arraySignal.', dMic, fs,  DOA_est, 'DSB', false);
audiowrite('output/x_dsb.wav', x_das, fs);

x_mvdr = applyBeamforming(arraySignal.', dMic, fs,  DOA_est, 'MVDR', false);
audiowrite('output/x_mvdr.wav', x_mvdr, fs);
%==========================================================================%

%====================Quality evaluation====================================%
%Compare
%[y,fs] = audioread('wav/speech1.wav'); 
%audiowrite('wav/speech1_16k.wav',y,16000);
%[y,fs] = audioread('output/x_dsb.wav'); 
%audiowrite('output/x_out_16k.wav',y,16000);
%%pesq = pesq('wav/speech1_16k.wav', 'wav/x_out_16k.wav')
%==========================================================================%