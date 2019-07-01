%DAS_BF
clear;
close all;
load('Computed_RIRs.mat');

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

%=============================STFT==============================%

frame_length = 1024;
frame_shift = 160;
fft_len = 1024;
[nFrames, ffts] = arrayStft(arraySignal, frame_length, frame_shift, fft_len);
%[frames, ffts] = multi_fft(arraySignal, frame_length, frame_shift, fft_len);
%nFrames = size(frames, 2)
figure(1);
subplot(3,1,1);
plotSpectrogram((1:nFrames)*frame_length/fs, (1:fft_len/2)*fs/fft_len, squeeze(ffts(1,:,:)));

[nFrames, spec_s] = arrayStft(speech1(:, 1), frame_length, frame_shift, fft_len);
subplot(3,1,2);
plotSpectrogram((1:nFrames)*frame_length/fs, (1:fft_len/2)*fs/fft_len, squeeze(spec_s(1,:,:)));
%==========================================================================%

%=============================apply CGMM=============================%
lambda_y = zeros(size(ffts, 2), size(ffts, 3));
lambda_n = zeros(size(ffts, 2), size(ffts, 3));

mini_batch = size(ffts,2); %estimate cgmm per  mini_batch 
iters = 10;
for t = 1:mini_batch:size(ffts, 2)
    idx =  t:min(t+mini_batch-1, size(ffts,2));
    [lambda_y_one, lambda_n_one,  Ry, Rn, Q] = cgmm_em(ffts(:, idx, :), iters);
    lambda_y(idx, :) = lambda_y_one;
    lambda_n(idx, :) = lambda_n_one;
end
%==========================================================================%
figure(2);
subplot(2,1,1);
plotMask(lambda_y.');
subplot(2,1,2);
plotMask(lambda_n.');
Rx = Ry -Rn;      

[M, T, F]  = size(ffts); %fft bins number
d = zeros(M, F);         %steering vectors
hw = d;                   %mvdr beamforming weight 
output = zeros(T, F);    %beamforming outputs

%========== Steering vectors esimation  ===========%
for f= 1:F
    d(:, f) = rtf(squeeze(Ry(:, :, f)), squeeze(Rn(:, :, f)), 'EVD');
    hw(:, f) = mvdr(d, Rn(:, :, f));
    output(:, f) =   hw(:, f) .'* squeeze(ffts(:, :, f));
end
figure(1);
subplot(3,1,3)
plotSpectrogram((1:nFrames)*frame_length/fs, (1:fft_len/2)*fs/fft_len, squeeze(output(1,:,:)));

%====================Quality evaluation====================================%
%output = [output, fliplr(conj(output(:, 2:end-1)))];
%rec_frames = real(ifft(output, fft_len, 2));
%rec_frames = rec_frames(:,1:frame_length);
%sig = overlapadd(rec_frames, hamming(frame_length, 'periodic'), frame_shift);
sig = invStft(output, frame_length, frame_shift);
audiowrite('output/output_cgmm.wav', sig./max(abs(sig)), fs);
save('output/testCGMM.mat')
%==========================================================================%