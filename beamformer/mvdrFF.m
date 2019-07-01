%MVDR_FF: Minimum Variance Distortionless Response Beamformer using
%Free-field Steering Vector
%--------------------
% Input parameters
%--------------------
%nChannels :   The steering vector, or relative transfer function [numChannels x numFreqs ]
%dMic : The distance between two microphones [a real positive value]
%DOA: The direction of arrival in degree [ 0~360 ]
%freqVec: The frequencies vector [ numFreqs]
%Phi_u: The PSD matrix of the interfering signal [numChannels x numChannels x numFreqs ]
%--------------------
% Output parameters
%--------------------
%H_mvdr: An mvdr beamformer [ numChannels x numFreqs ] 
function H_mvdr = mvdrFF(nChannels, dMic, DOA, freqVec, Phi_u)

c = 340;
theta0 = DOA * pi / 180;                     
d = exp(-1j*[0:nChannels-1]'*2*pi*freqVec*cos(theta0)*dMic/c) ./ nChannels;  % steering vector

H_mvdr = zeros(nChannels, length(freqVec));
%Phi_u_inv = zeros(nChannels, nChannels);
for f = 1:length(freqVec)  
   Phi_u_inv = pinv(Phi_u(:,:,f));
   df = d(:, f);
   H_mvdr(:, f) = (Phi_u_inv * df) / (df' * Phi_u_inv * df);
end

end