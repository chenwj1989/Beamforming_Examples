%CGMM_EM: Complex Gaussian Mixture Model - EM Algorithm
%--------------------
%Input parameters
%--------------------
%stft:   The multi-channel time-frequency spectrum for one mini-batch [numMics x numFrames x numFreqs ]
%--------------------
%Output parameters
%--------------------
%nPSD: The updated PSD matrix of the interfering signal [numMics x numFrames x numFreqs ]
%--------------------
%Reference: 
%--------------------
%Higuchi, Takuya, et al. "Robust MVDR beamforming using time-frequency masks for online/offline ASR in noise." 
%2016 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP). IEEE, 2016.

function [lambda_xn, lambda_n, R_xn, R_n, Q] = cgmm_em(y_stft, iters)

%Get input spectrum dimention
[numChannels, numFrames, numFreqs] = size(y_stft);

%Calculate covariance matrix of input signal 
yyh = zeros(numChannels, numChannels, numFrames, numFreqs);
yyh_norm = zeros(numChannels, numChannels, numFrames, numFreqs);


%initilization
lambda_xn = ones(numFrames, numFreqs);
lambda_n = zeros(numFrames, numFreqs);
p_xn = ones(numFrames, numFreqs);
p_n = ones(numFrames, numFreqs);
phi_xn = ones(numFrames, numFreqs);
phi_n = ones(numFrames, numFreqs);

for t = 1:numFrames
    for f = 1:numFreqs
     yyh(:, :, t, f) = y_stft(:, t, f) * y_stft(:, t ,f)';  
     yyh_norm(:, :, t, f) = yyh(:, :, t, f) /  yyh(1, 1, t, f);  
    end
end
R_xn = squeeze(mean(yyh, 3));
R_n = repmat(eye(numChannels), [1, 1, numFreqs]);
R_nn = repmat(eye(numChannels), [1, 1, numFreqs]);

%R_xn = zeros(numChannels, numChannels, numFreqs);
%R_n = zeros(numChannels, numChannels, numFreqs);


%
theta = 0.0001;
d = 0.0001;
Q = zeros(iters);

for i = 1 : iters
    for f = 1:numFreqs
        
        acc_n = zeros(numChannels, numChannels);
        R_xn_sum = zeros(numChannels, numChannels);
        R_n_sum = zeros(numChannels, numChannels);
        
        Rf_xn = R_xn(:, :, f);
        if rcond(Rf_xn) < theta
            Rf_xn = Rf_xn + rand(numChannels)*d;
        end            
        Rf_n = R_n(:, :, f); 
        if rcond(Rf_n) < theta
            Rf_n = Rf_n + rand(numChannels)*d;
        end  
        invRf_xn = inv(Rf_xn);
        invRf_n = inv(Rf_n);

        for t = 1:numFrames       
            y = y_stft(:, t, f);
            yyh_one = yyh(:, :, t, f);
            
            % EM
            phi_xn(t, f) = trace(yyh_one * invRf_xn) / numChannels; %eq(11)
            phi_n(t, f) = trace(yyh_one * invRf_n) / numChannels; %eq(11)
            
            kxn = y' * (1 / phi_xn(t, f)) * invRf_xn * y;
            p_xn(t, f) = real( exp(-kxn) / (pi)^numChannels / det(phi_xn(t, f)*Rf_xn));  %eq(8)
            kn = y' * (1 / phi_n(t, f)) * invRf_n * y;
            p_n(t, f) = real( exp(-kn) / (pi)^numChannels  / det(phi_n(t, f)*Rf_n));

            lambda_xn(t, f) = p_xn(t, f) / (p_xn(t, f) + p_n(t, f)); %eq(10)
            lambda_n(t, f) = p_n(t, f) / (p_xn(t, f) + p_n(t, f));

            R_xn_sum = R_xn_sum  + lambda_xn(t, f) / phi_xn(t, f) * yyh_one; %eq(12)
            R_n_sum =  R_n_sum  + lambda_n(t, f) / phi_n(t, f) * yyh_one; %eq(12)
            acc_n = acc_n + lambda_n(t, f)*yyh_one; %for eq(4)
        end
        R_nn(:, :, f) = 1/sum(lambda_n(:, f)) * acc_n; %eq(4)
        
        R_xn(:, :, f) = R_xn_sum / sum(lambda_xn(:, f)); %eq(12)
        R_n(:, :, f) = R_n_sum / sum(lambda_n(:, f)); %eq(12)
        
        tmp_Rxn_f = 1/sum(lambda_xn(:, f)) * R_xn_sum;
        tmp_Rn_f = 1/sum(lambda_n(:, f)) * R_n_sum;
       [V1, ~] = eig(squeeze(tmp_Rxn_f));
       [V2, ~] = eig(squeeze(tmp_Rn_f));

        entropy1 = -diag(V1, 0)'/sum(diag(V1, 0)) * log(diag(V1, 0)/sum(diag(V1, 0)));
        entropy2 = -diag(V2, 0)'/sum(diag(V2, 0)) * log(diag(V2, 0)/sum(diag(V2, 0)));
        if entropy1 > entropy2
            R_xn(:, :, f) = tmp_Rn_f;
            R_n(:, :, f) = tmp_Rxn_f;
        else
            R_xn(:, :, f) = tmp_Rxn_f;
            R_n(:, :, f) = tmp_Rn_f;
        end
        
    end %numFreqs
    
    Qn = sum(sum(lambda_n .* log(p_n + 0.00001))) / (numFrames * numFreqs);
    Qx = sum(sum(lambda_xn .* log(p_xn + 0.00001))) / (numFrames * numFreqs);
    Q(i) = Qx + Qn;
    fprintf('--- iter = %2d, Q = %.4f + %.4f = %.4f\n', i, Qn, Qx, Q(i));

end %iter

%R_n = R_nn;
end
    





