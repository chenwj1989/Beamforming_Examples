function [phi] = tdoa(arraySignals, fs, d, c)

x1 = arraySignals(:,1);
x2 = arraySignals(:,2);

[acor,lag] = xcorr(x2,x1);
[~,I] = max(abs(acor));
t = lag(I);   

phi = acos(t/fs*c/d) * 180 / pi;

end

function [phi] = gcc_phat(arraySpectrum, fs, d, c)

X1 = arraySpectrum(:,1);
X2 = arraySpectrum(:,2);

G = X1 .* conj(X2);
corr = ifft(G ./ abs(G));

[~,I] = max(abs(corr));
tdoa = lag(I);   

phi = acos(tdoa/fs*c/d) * 180 / pi;

end