%DSB: Delay and Sum Beamformer
%--------------------
% Input parameters
%--------------------
%nChannels :   The steering vector, or relative transfer function [numChannels x numFreqs ]
%dMic : The distance between two microphones [a real positive value]
%DOA: The direction of arrival in degree [ 0~360 ]
%freqVec: The frequencies vector [ numFreqs]
%--------------------
% Output parameters
%--------------------
%H_dsb: A delay-and-sum beamformer [ numChannels x numFreqs ] 
function H_dsb = dsb(nChannels, dMic, DOA, freqVec)

c = 340;
theta0 = DOA * pi / 180;                     
d = exp(-1j*2*pi*freqVec'*[0:nChannels-1]*cos(theta0)*dMic/c) ./ nChannels;  % steering vector
H_dsb = d';

end