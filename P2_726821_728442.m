clear all;
%****** Generating Sine waves********%
Fs = 96000; % Sampling Fs 
Ts = 1/Fs;
F_1 = 5000;  % SIN_1 frequency
F_2 = 10000; % SIN_2 frequency
T = 10;      % Duration
t = 0:Ts:(T-Ts);    % Vector time
Sin_5k  = sin(2*pi*F_1*t); 
Sin_10k = sin(2*pi*F_2*t);

%%
%***** Transmiting signals ********%
soundsc(Sin_5k,Fs);
soundsc(Sin_10k,Fs);
%%
%****** Receiving Sine waves********%


Fs=96000; Nb=16;Chs=1; 
recObj = audiorecorder(Fs, Nb, Chs); 
get(recObj); 
disp('Start speaking.') 
recordblocking(recObj, 5); 
disp('End of Recording.'); 
play(recObj); % Play back the recording. 
% Store data in double-precision array. 
myRecording = getaudiodata(recObj); 
% Plot the waveform. 
plot(myRecording); 
title('Recording plot');
% Power Spectrum Densitiy: 
figure();
pwelch(myRecording,500,300,500,'one-side','power',Fs) 
%%
%sendign a pulse 
soundsc( [zeros(1,Fs) 1 zeros(1,Fs)], Fs );

%%
%Sending AWGN 
xa = ones(1,Fs*5);
xa_noised  = awgn(xa,10);
plot(xa_noised);

soundsc(xa_noised,Fs);


%%
%Generate a Power Signal: Chirp Signal for example
start_T=0;
end_T=5;
t = start_T:1/Fs:end_T;
fo = 20;
f1 = 20e3;
y = chirp(t,fo,end_T,f1,'linear');
figure(1)
pwelch(y,500,300,500,'one-side','power',Fs)
plot(y);

%sending signal
soundsc(y,Fs);