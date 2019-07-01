%GEV: Generalized Eigenvalue Beamformer
%--------------------
% Input parameters
%--------------------
%Phi_u: The PSD matrix of the interfering signal [numChannels x numChannels x numFreqs ]
%Phi_y: The PSD matrix of the recorded signal [numChannels x numChannels x numFreqs ]
%--------------------
% Output parameters
%--------------------
%H_mvdr: An mvdr beamformer [ numChannels x numFreqs ] 
%--------------------
%Reference: 
%--------------------
%Warsitz, E., & Haeb-Umbach, R. (2007). Blind acoustic beamforming based on generalized eigenvalue decomposition. 
%IEEE Transactions on audio, speech, and language processing, 15(5), 1529-1539.

function H_gev = gev(Phi_u, Phi_y)

[numChannels, ~, numFreqs] = size(Phi_u);
H_gev = ones(numChannels, numFreqs);

for f = 1:numFreqs
    Phi_u_inv = pinv(Phi_u(:, :, f));
    [Vmax, ~] = eigs( Phi_u_inv * Phi_y(:, :, f), 1 ); % the eigen vector with biggest eigen value
    H_gev(:, f) = Vmax;
end

end

