%estimate_noise_PSD: 
%--------------------
% Input parameters
%--------------------
%spectrum:   The steering vector, or relative transfer function [numMics x num_frequencies ]
%old_nPSD: The old PSD matrix of the interfering signal [num_mics x num_mics x num_freq ]
%mask: The vad or speech presence probability [a real number between 0~1 ]
%alpha: The adaptation parameter to update noise PSD [a real number between 0~1 ]
%--------------------
% Output parameters
%--------------------
%nPSD: The updated PSD matrix of the interfering signal [num_mics x num_mics x num_freq ]
function nPSD = estimate_noise_PSD(spectrum, old_nPSD, mask, alpha)

% Memory allocation
num_freq = size(spectrum,2);

% initialization
nPSD = old_nPSD;

for freqidx = 1:num_freq

    % select the signal for the current time-frequency bin
    currentSig = spectrum(:,freqidx);

    % compute the average parameter based on the mask
    currentAlpha = alpha + (1 - mask(freqidx)) * (1 - alpha);

    % Perform the recursive step for the current TF-bin
    nPSD(:,:,freqidx) = currentAlpha * old_nPSD( :, :, freqidx) + (1 - currentAlpha) * (currentSig * currentSig');

end

end
    

