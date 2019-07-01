function [nFrames, stft] = arrayStft(arraySignal, frame_size, frame_shift, nFFT)

x = arraySignal.'; 
[nChannels, len] = size(x); % input signal, matrix  [nChannels * time]

nFrames = floor((len - frame_size) / frame_shift);
if rem(nFFT, 2)==1, 
    nFFT = nFFT + 1; 
end;

win = hamming(frame_size, 'periodic')';  % define window
win2D = repmat(win, nChannels, 1);
stft = zeros(nChannels, nFrames,  nFFT/2+1);

k = 1;
for t = 1 : nFrames
    xk = win2D .* x(:, k:k+frame_size-1);
    tmp = fft(xk, nFFT, 2);
    stft(:, t, :) = tmp(:, 1:nFFT/2+1);
    k = k + frame_shift; 
end