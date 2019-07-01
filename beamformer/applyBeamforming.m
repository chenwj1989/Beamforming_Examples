
function [x_out] = applyBeamforming(arraySignal, dMic, fs, DOA, typeBF, flagGSC)

nChannels = size(arraySignal,1);
x = arraySignal;  % input signal, matrix  [nChannels * time]

%=========initialize parameters for spectral processing ===================%
len=floor(20 * fs / 1000); % Frame size in samples
if rem(len,2)==1, 
    len=len+1; 
end;
overlap = 0.50; % window overlap in percent of frame size
len1 = floor(len * overlap);
len2 = len - len1;

nFrames = floor(size(x, 2)/len2) - floor(len/len2);
xfinal = zeros(1, nFrames*len2+len1);

win = hamming(len)';  % define window
win2D = repmat(win, nChannels, 1);

nFFT = len;
f = linspace(0, fs, nFFT);% freqeuncy vector
%==========================================================================%

%============================Initialize covariance matirx==================%
noise_mean=zeros(nChannels, nFFT);
noiseConv = zeros(nChannels, nChannels, nFFT);
j=1;
for k=1:10   %use fisrt 10 segments to  estimate noise spectrum
    noise_spec = fft(win2D .* x(:, j:j+len-1), nFFT, 2);
    noise_mean = noise_mean + noise_spec;
    for freqidx = 1:nFFT   
        noiseConv(:,:,freqidx) = noiseConv(:,:,freqidx) + (noise_spec(:, freqidx)' * noise_spec(:, freqidx));
    end
    j = j + len2;
end
noiseConv = noiseConv / 10;
%==========================================================================%

%================= generate beamformer ====================================%
hw = zeros(nChannels, nFFT);
if strcmp(typeBF, 'DSB')
    hw = dsb(nChannels, dMic, DOA, f);        % delay-and-sum filter
elseif  strcmp(typeBF, 'MVDR')
    hw = mvdrFF(nChannels, dMic, DOA, f, noiseConv);
end
%==========================================================================%

%================= speech enhancement per frame============================%
k = 1;
for n = 1 : nFrames
  
    xk = win2D .* x(:, k:k+len-1);
    spec = fft(xk, nFFT, 2);
    
    %vad or speech presence probability
    %power covariance matrix
    %updatebeamformer
    %rtf = rtf(, noiseConv)
    %hw = mvdr(rtf, noiseConv);
    
    if flagGSC == false
        spec_bf = sum(hw .* spec, 1);
    else
        [spec_bf, ak] = lmsGSC(X(:, 1:(L/2)+1), hw, 0, B, ak, mu);
    end
        
    xi_w = ifft( spec_bf, nFFT);
    xi_w = real(xi_w);
    
   % Overlap and add 
    xfinal( k:k+len-1) =  xfinal(k:k+len-1) + win.* xi_w;
    
    k = k + len2;
    
end
%==========================================================================%

%===============================Output=====================================%
x_out = xfinal;
%==========================================================================%

end