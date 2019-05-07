%RTF: Relative Transfer Function using Covariace Subtraction
%--------------------
% Input parameters
%--------------------
%Ry :   The PSD matrix of input signal,  [numMics x num_mics x num_freq ]
%Rn :   The PSD matrix of estimated noise,  [numMics x num_mics x num_freq ]
%--------------------
% Output parameters
%--------------------
%h: The estimated reative transfer function [ numMics x num_freq ] 
function [h] = rtf(Ry, Rn)
    N =  size(Ry, 1);
    Rx = Ry - Rn;
    e =  [1; zeros(N-1, 1)];
    h = (Rx * e) ./ (e.' * Rx * e);


