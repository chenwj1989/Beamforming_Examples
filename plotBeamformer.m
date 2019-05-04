function [res] = plotBeamformer(M, f, d, bfFilter, theta)
    c = 340;
    w = bfFilter;       % beamformer
    p = zeros(length(theta), length(f));
    for  j=1:length(theta)                    %scan angles                 
        a=exp(-1j*[0:M-1]'*2*pi*f*cos(theta(j))*d/c);
        p(j, :) = sum(w.*a);                                       
    end
   % p=p';
    res = p;
    
    if length(f) == 1
        figure;
        subplot(1, 2, 1);
        plot(theta/pi*180, abs(p.')); grid on;
        xlabel('Degree')
        subplot(1, 2, 2);
        polar(theta,abs(p.'))
    else
        figure;
        subplot(1, 2, 1);
        mesh(f, theta/pi*180, abs(p))
        xlabel('Frequency(Hz)')
        ylabel('Degree')
        title('Beamformer Directivity')
        subplot(1, 2, 2);
        imagesc(abs(p).');
    end

end

