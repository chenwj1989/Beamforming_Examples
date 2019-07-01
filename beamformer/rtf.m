%RTF: Relative Transfer Function
%--------------------
% Input parameters
%--------------------
%Ry :   The PSD matrix of input signal,  [numChannels x numChannels x numFreqs ]
%Rn :   The PSD matrix of estimated noise,  [numChannels x numChannels x numFreqs ]
%type : The type of method used to estimate RTF, [string] 
%--------------------
% Output parameters
%--------------------
%h: The estimated reative transfer function [ numChannels x numFreqs ] 
%--------------------
% Reference
%--------------------
%Gannot, S., Vincent, E., Markovich-Golan, S., Ozerov, A., Gannot, S., Vincent, E., ... & Ozerov, A. (2017). A consolidated perspective on multimicrophone speech enhancement and source separation. 
%IEEE/ACM Transactions on Audio, Speech and Language Processing (TASLP), 25(4), 692-730.

function [h] = rtf(Ry, Rn, type)
[numChannels, ~, numFreqs] = size(Ry);
h = ones(numChannels, numFreqs);
e =  [1; zeros(numChannels-1, 1)];
    
if strcmp(type, 'CS') %covariance subtraction
    Rx = Ry - Rn;
    for f = 1:numFreqs
        h(:, f) = (Rx(:, :, f) * e) ./ (e.' * Rx(:, :, f) * e);
    end
    
elseif  strcmp(type, 'CW') %covariance whitening
    for f = 1:numFreqs
        invRn = pinv(Rn(:, :, f));
        [q, ~] = eigs(Ry(:, :, f), 1);
        h(:, f) = (e'*invRn*q) \ (invRn*q);
    end
    
elseif  strcmp(type, 'EVD') %eigen vector decomposition
    Rx = Ry - Rn;
    for f = 1:numFreqs
         [q, ~] = eigs(Rx(:, :, f), 1);
        h(:, f) = q;
    end

end

