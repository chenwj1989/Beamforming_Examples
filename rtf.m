%RTF: Relative Transfer Function using Covariace Subtraction
%--------------------
% Input parameters
%--------------------
%Ry :   The PSD matrix of input signal,  [numMics x numMics x numFreqs ]
%Rn :   The PSD matrix of estimated noise,  [numMics x numMics x numFreqs ]
%--------------------
% Output parameters
%--------------------
%h: The estimated reative transfer function [ numMics x numFreqs ] 
function [h] = rtf(Ry, Rn)
    numMics =  size(Ry, 1);
    numFreqs = size(Ry, 3);
    Rx = Ry - Rn;
    
    e =  [1; zeros(numMics-1, 1)];
    h = zeros(numMics, numFreqs);
    for i = 1:numFreqs
        h(:, i) = (Rx(:, :, i) * e) ./ (e.' * Rx(:, :, i) * e);
    end
    


