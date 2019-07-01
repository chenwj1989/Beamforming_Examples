
% To Generate a room impulse response,
% an open-souce program "Room Impulse Response Generator (Version 2.1.20141124) by
% Emanuel Habets" is used. The source codes locate at:
% https://github.com/ehabets/RIR-Generator.git
% 
% The program is already compiled and put at
% contribute/rir_generator.mexw64, ready for matlab to call by api:
% rir_generator(c, fs, r, s, L, beta, n).


m_pos = [4 1.96 2 ; 
         4 1.98 2 ; 
         4 2.00 2 ;
         4 2.02 2 ;
         4 2.04 2 ; ];      % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)   
s_pos = [2 0 2 ;
         2 4 2 ; ];         % Source position [x y z] (m)

c = 340;                    % Sound velocity (m/s)
fs_RIR = 44100;             % Sample frequency (samples/s)
room_dim = [5 4 6];         % Room dimensions [x y z] (m)
rev_time = 0;             % Reverberation time (s)
n = 22050;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

RIR_sources = zeros(n, size(m_pos, 1), size(s_pos, 1));
for i = 1:size(s_pos, 1)
h = rir_generator(c, fs_RIR, m_pos, s_pos, room_dim, rev_time, n, mtype, order, dim, orientation, hp_filter);
RIR_sources(:, :, i) = h.';
end

save('Computed_RIRs','RIR_sources','fs_RIR','m_pos','s_pos','room_dim','rev_time')
