%MVDR_FF: Minimum Variance Distortionless Response Beamformer using
%Free-field Steering Vector
%--------------------
% Input parameters
%--------------------
%nChannels :   The steering vector, or relative transfer function [numMics x num_frequencies ]
%dMic : The distance between two microphones [a real positive value]
%DOA: The direction of arrival in degree [ 0~360 ]
%freqVec: The frequencies vector [ num_frequencies]%Phi_u: The PSD matrix of the interfering signal [num_mics x num_mics x num_freq ]
%--------------------
% Output parameters
%--------------------
%H_mvdr: An mvdr beamformer [ num_frequencies x num_freq ] 
function H_mvdr = mvdrFF(nChannels, dMic, DOA, freqVec, Phi_u)

c = 340;
theta0 = DOA * pi / 180;                     
d = exp(-1j*[0:nChannels-1]'*2*pi*freqVec*cos(theta0)*dMic/c) ./ nChannels;  % steering vector

H_mvdr = zeros(nChannels, length(freqVec));
%Phi_u_inv = zeros(nChannels, nChannels);
for freqidx = 1:length(freqVec)  
   Phi_u_inv = pinv(Phi_u(:,:,freqidx));
   H_mvdr(:, freqdix) = (Phi_u_inv * d) / (d.' * Phi_u_inv * d);
end

end