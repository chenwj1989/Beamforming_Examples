%GSC: Generalized Sidelobe Cancellation (GSC) beamformer
%-----------------
% Input parameters
%-----------------
% X:  input spectrogram from each microphone, a matrix of size [numChannels x num_frequencies]
% H:  main beamformer, a matrix of size [numChannels x num_frequencies]
% phi_delta:  an extra delay on main beamformer, [a real number]
% B:  blocking matrix, a matrix of size [numChannels - 1 x numChannels]
% W:  current sidelobe canceller coefficients, a matrix of size [numChannels-1 x num_frequencies]
% mu: lms update rate, [a real number]
%-----------------
% Output parameters
%-----------------
% Xout:  output signal after GSC beamformer , an array of lenght [ num_frequencies]
% Wnext: updated  sidelobe canceller coefficients, a matrix of size [numChannels x num_frequencies]
%-----------------
function [Xout, Wnext] = lmsGSC(X, H, phi_delta, B, W, mu)

%main beamformer
HT = H.';
Xtemp = HT .* X;
Xmain = sum(Xtemp, 1) * exp(-1i * phi_delta);

%sidelobe 
Xb = B * X;
Xs = sum(W .* Xb, 1); 

%lms
N = size(Xb, 1);
Xout = Xmain - Xs;
err = repmat(Xout, N, 1) .* Xb;
Wnext = W + mu * err ./ Xb.^2;   
  
end
