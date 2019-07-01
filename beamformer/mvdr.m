%MVDR: Minimum Variance Distortionless Response Beamformer
%--------------------
% Input parameters
%--------------------
%d :   The steering vector, or relative transfer function [numChannels x numFreqs ]
%Phi_u: The PSD matrix of the interfering signal [numChannels x numChannels x numFreqs ]
%--------------------
% Output parameters
%--------------------
%H_mvdr: An mvdr beamformer [ numChannels x numFreqs ] 
function H_mvdr = mvdr(d, Phi_u)

[numChannels, ~, numFreqs] = size(Phi_u);
H_mvdr = ones(numChannels, numFreqs);

for f = 1:numFreqs
    Phi_u_inv = pinv(Phi_u(:, :, f));
    df = d(:, f);
    H_mvdr(:, f) = (Phi_u_inv * df) / (df' * Phi_u_inv * df);
end

end

