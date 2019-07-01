function plotSpectrogram(frameTime, freq, y)


%set(gcf,'Position',[20 100 600 500]);            
%axes('Position',[0.1 0.1 0.85 0.5]);  
imagesc(frameTime, freq, abs(y)');  
axis xy; ylabel('Freq/Hz');xlabel('Time/s');
title('Spectrogram');
m = 64;
LightYellow = [0.6 0.6 0.6];
MidRed = [0 0 0];
Black = [0.5 0.7 1];
Colors = [LightYellow; MidRed; Black];
%colormap(SpecColorMap(m,Colors));
end

