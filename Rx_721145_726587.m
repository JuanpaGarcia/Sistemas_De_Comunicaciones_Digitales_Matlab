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
%EL RECEPTOR TAMBIEN CONOCE EL PULSO FORMADOR
%modificar el valor de B en 4000 y12000 para los diferentes ejercicios
Fs      =   96e3;              % Samples per second  
Ts      =   1/Fs;              % Sampling period 
beta    =   0.25;              % Roll-off factor 
B       =   4000;              % Bandwidth available 
Rb      =   2*B/(1+beta);      % Bit rate = Baud rate 
mp      =   ceil(Fs/Rb)        % samples per pulse 
Rb      =   Fs/mp;             % Recompute bit rate 
Tp      =   1/Rb;              % Symbol period 
B       =   (Rb*(1+beta)/2)    % Bandwidth consumed 
D       =   10;                % Time duration in terms of Tp 
type    =   'srrc';            % Shape pulse: Square Root Rise Cosine 
E       =   Tp;                % Energy 
[pbase ~] = rcpulse(beta, D, Tp, Ts, type, E);    % Pulse Generation 


%%
%forma1: se puede utilizar un archivo WAV o escuchar la señal directamante y guardarla

%filename= 'B4K.wav'; 
%info=audioinfo(filename); 
%b = info.BitsPerSample; 
%Fs_audio = info.SampleRate; 
%sec = 5; start = 1; 
%samples = [start*Fs_audio,(sec+start)*Fs_audio-1]; 
%[Rx_signal,Fs_audio] = audioread(filename,samples); 

%forma2:

Fs=96e3; sec=5;  % Time duration of the whole communication including the silence 
recObj = audiorecorder(Fs,16,1); 
disp('Start speaking'); 
recordblocking(recObj, sec); 
disp('End of Recording'); 
Rx_signal = getaudiodata(recObj); 


threshold = 0.1;                            % Detecting the channel energization 
start = find(abs(Rx_signal)> threshold,3,'first'); % Initial 
stop  = find(abs(Rx_signal)> threshold,1,'last');  % End 
Rx_signal = Rx_signal (start:stop); 

match = fliplr(pbase);%el pulso formador tambien lo conoce el receptor
Rx_signal_match = conv(match, Rx_signal);
eyediagram(Rx_signal_match,mp*3);
figure();
pwelch(Rx_signal_match,500,300,500,'one-side','power',Fs)

figure();
plot(Rx_signal_match(1:mp*128))

%cambiar el inicio del preambulo dependiendo de cual B se esta utilizando
start_preamble = 94;
preamble = sign(Rx_signal_match(start_preamble:mp:56*mp+(start_preamble-1)));
preamble = (preamble + 1)/2
%preamble'

start_sfd = 56*mp + start_preamble;
sfd = sign(Rx_signal_match(start_sfd:mp:56*mp+(start_preamble-1)+8*mp));
sfd = (sfd + 1)/2
%sfd'

start_sda = start_sfd +8*mp;
dsa = sign(Rx_signal_match(start_sda:mp:56*mp+(start_preamble-1)+8*mp+280*mp));
dsa = (dsa + 1)/2
%sda'
end_preamble = 56*mp+(start_preamble-1)+8*mp+280*mp;

Rx_signal_match = Rx_signal_match(1:numel(bits2Tx));
Rx_signal_match = Rx_signal_match';

bitsM1= Rx_signal_match(start_preamble:end_preamble);
bitsM1 = reshape(bitsM1,16,2);
bitsM1 = bitsM1';
decval1 = bi2de(bitsM1,'left-msb');
lenaRS1 = reshape(decval1, size(bits2Tx));
imshow(uint8(lenaRS1))
