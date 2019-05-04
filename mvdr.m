%MVDR: Minimum Variance Distortionless Response Beamformer
%--------------------
% Input parameters
%--------------------
%d :   The steering vector, or relative transfer function [numMics x num_frequencies ]
%Phi_u: The PSD matrix of the interfering signal [num_mics x num_mics x num_freq ]
%--------------------
% Output parameters
%--------------------
%H_mvdr: An mvdr beamformer [ num_frequencies x num_freq ] 
function H_mvdr = mvdr(d, Phi_u)

Phi_u_inv = pinv(Phi_u);
H_mvdr = (Phi_u_inv * d) / (d' * Phi_u_inv * d);

end

