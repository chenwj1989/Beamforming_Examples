function [istft] = invStft(ffts, frame_size, frame_shift)

x = [ffts, fliplr(conj(ffts(:, 2:end-1)))];
[nFrames, nFFT] = size(x);

win = hamming(frame_size, 'periodic')';  % define window
istft = zeros(1, (nFrames-1)*frame_shift+frame_size);

k = 1;
for t = 1 : nFrames-1
    tmp = real(ifft(x(t, :), nFFT));
    istft(:,  k:k+frame_size-1) = istft(:,  k:k+frame_size-1)  + tmp .* win;
    k = k + frame_shift; 
end