%******                     RECEPTION                     ********%
%************************** PARTE I ******************************%
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
%************************** PARTE II ******************************%
%**************************   4KHz   ******************************%
% RECEIVING SIGNAL
Fs = 96e3; 
file = 'AP2.wav'; % SAVED AUDIO IN THE TRANSMITION 
FILE_INFO = audioinfo(file);
[Rx_signal,Fs] = audioread(file); 
sec = FILE_INFO.Duration;
threshold = 0.1;                           
start = find(abs(Rx_signal)> threshold,3,'first'); % Initial 
stop  = find(abs(Rx_signal)> threshold,1,'last');  % End 
Rx_signal = Rx_signal (start:stop); 
% CATCHING WITH MATCH FILTER 
MF = fliplr(pbase);
Mfil = conv(Rx_signal,MF);
plot(Mfil(1:6000)); % To view preamble, SFD, DSA and header bits
eyediagram(Mfil(1:5000), mp*3);
figure;
pwelch(Mfil,[],[],[],Fs,'power'); % PSD
title('PSD Non wait time PNRZ recovered signal');
%%
% SAMPLING

start =  96; % Starting medition point 
sampled_signal = Mfil(start:mp:end); % Sample every mp
bits_Rxp = zeros(1,numel(sampled_signal));
bits_Rxp(sampled_signal >= 0) = 1;
bits_Rxp(sampled_signal < 0) = 0;
plot(bits_Rxp(1:1000)); 
figure;
bits_Rxp = bits_Rxp(1:numel(bits2Tx));
bits_Rxp = bits_Rxp';

% DETECTING ONLY SIGNAL BITS

% Creación del objeto
header_detect = comm.PreambleDetector(header_bits,'Input','Bit');
preamble_detect = comm.PreambleDetector(trash_bits,'Input','Bit');   
% Deteccion de preámbulo. El índice indica dónde termina la trama
idx = preamble_detect(bits_Rxp);  
idx_h = header_detect(bits_Rxp); 
% Una vez que encuentra el índice, se descartan los “bits basura”
% Una forma de hacerlo es la siguiente:
bits_Rxp_h = bits_Rxp(idx_h+1:idx_h+32);
bits_Rxp_h = reshape(bits_Rxp_h,16,2);
bits_Rxp_h = bits_Rxp_h';
bits_Rxp_h = bi2de(bits_Rxp_h,'left-msb');
bits_Rxp= bits_Rxp(idx+1:end);
bits_Rxp(:);
error =  sum(xor(bits2Tx(idx+1:end),bits_Rxp));
BER = (error/numel(bits2Tx)) * 100; % Calculating Bit Error Rate 

% REBUILDING IMAGE WITH RECEIVED BITS 

recuperado = zeros(bits_Rxp_h(1),bits_Rxp_h(2),'uint8'); % allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : bits_Rxp_h(2)
    for j = 1: bits_Rxp_h(1)
        recuperado(j,i) = bi2de(bits_Rxp(counter-7:counter)','left-msb');
        counter = counter +8;
    end
end

imshow(uint8 (recuperado));
%%
%**************************   12KHz   ******************************%
% RECEIVING SIGNAL
Fs = 96e3;  
file = 'exa2.wav'; % SAVED AUDIO IN THE TRANSMITION 
FILE_INFO = audioinfo(file);
[Rx_signal,Fs] = audioread(file); 
sec = FILE_INFO.Duration;
threshold = 0.1;                            % Detecting the channel energization 
start = find(abs(Rx_signal)> threshold,3,'first'); % Initial 
stop  = find(abs(Rx_signal)> threshold,1,'last');  % End 
Rx_signal = Rx_signal (start:stop); 
MF = fliplr(pbase);
Mfil = conv(Rx_signal,MF);
plot(Mfil(1:2000)); % To view preamble, SFD, DSA and header bits
eyediagram(Mfil(1:3000), mp*3);
figure;
pwelch(Mfil,[],[],[],Fs,'power'); % PSD
title('PSD Non wait time PNRZ recovered signal');
%%
% SAMPLING
start =  32; % Starting medition point by eye of good cuber
sampled_signal = Mfil(start:mp:end); % Sample every mp
 
bits_Rxp = zeros(1,numel(sampled_signal));
bits_Rxp(sampled_signal >= 0) = 1;
bits_Rxp(sampled_signal < 0) = 0;
bits_Rxp = bits_Rxp(1:numel(bits2Tx));
bits_Rxp = bits_Rxp';

% DETECTING ONLY SIGNAL BITS

% Creación del objeto
header_detect = comm.PreambleDetector(header_bits,'Input','Bit');
preamble_detect = comm.PreambleDetector(trash_bits,'Input','Bit');   
% Deteccion de preámbulo. El índice indica dónde termina la trama
idx = preamble_detect(bits_Rxp);  
idx_h = header_detect(bits_Rxp); 
% Una vez que encuentra el índice, se descartan los “bits basura”
% Una forma de hacerlo es la siguiente:
bits_Rxp_h = bits_Rxp(idx_h+1:idx_h+32);
bits_Rxp_h = reshape(bits_Rxp_h,16,2);
bits_Rxp_h = bits_Rxp_h';
bits_Rxp_h = bi2de(bits_Rxp_h,'left-msb');
bits_Rxp= bits_Rxp(idx+1:end)
bits_Rxp(:);
error =  sum(xor(bits2Tx(idx+1:end),bits_Rxp));
BER = (error/numel(bits2Tx)) * 100;
% REBUILDING IMAGE WITH RECEIVED BITS 

recuperado = zeros(bits_Rxp_h(1),bits_Rxp_h(2),'uint8'); % allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : bits_Rxp_h(2)
    for j = 1: bits_Rxp_h(1)
        recuperado(j,i) = bi2de(bits_Rxp(counter-7:counter)','left-msb');
        counter = counter +8;
    end
end

imshow(uint8 (recuperado));
%%
%**************************   20KHz   ******************************%
% RECEIVING SIGNAL
Fs = 96e3;  % Time duration of the whole communication including the silence 
file = '20K.wav'; % SAVED AUDIO IN THE TRANSMITION 
FILE_INFO = audioinfo(file);
[Rx_signal,Fs] = audioread(file); 
sec = FILE_INFO.Duration;
threshold = 0.1;                            % Detecting the channel energization 
start = find(abs(Rx_signal)> threshold,3,'first'); % Initial 
stop  = find(abs(Rx_signal)> threshold,1,'last');  % End 
Rx_signal = Rx_signal (start:stop); 
MF = fliplr(pbase);
Mfil = conv(Rx_signal,MF);
eyediagram(Mfil(1:3000), mp*3);
figure;
pwelch(Mfil,[],[],[],Fs,'power'); % PSD
title('PSD Non wait time PNRZ recovered signal');
%%
start =  20; % Starting medition point by eye of good cuber
sampled_signal = Mfil(start:mp:end); % Sample every mp
scatterplot(sampled_signal);
figure;
bits_Rxp = zeros(1,numel(sampled_signal));
bits_Rxp(sampled_signal >= 0) = 1;
bits_Rxp(sampled_signal < 0) = 0;
bits_Rxp = bits_Rxp(1:numel(bits2Tx));
bits_Rxp = bits_Rxp';

% DETECTING ONLY SIGNAL BITS

% Creación del objeto
header_detect = comm.PreambleDetector(header_bits,'Input','Bit');
preamble_detect = comm.PreambleDetector(trash_bits,'Input','Bit');   
% Deteccion de preámbulo. El índice indica dónde termina la trama
idx = preamble_detect(bits_Rxp);  
idx_h = header_detect(bits_Rxp); 
% Una vez que encuentra el índice, se descartan los “bits basura”
% Una forma de hacerlo es la siguiente:
bits_Rxp_h = bits_Rxp(idx_h+1:idx_h+32);
bits_Rxp_h = reshape(bits_Rxp_h,16,2);
bits_Rxp_h = bits_Rxp_h';
bits_Rxp_h = bi2de(bits_Rxp_h,'left-msb');
bits_Rxp= bits_Rxp(idx+1:end)
bits_Rxp(:);
error =  sum(xor(bits2Tx(idx+1:end),bits_Rxp));
BER = (error/numel(bits2Tx)) * 100;
% REBUILDING IMAGE WITH RECEIVED BITS 

recuperado = zeros(bits_Rxp_h(1),bits_Rxp_h(2),'uint8'); % allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : bits_Rxp_h(2)
    for j = 1: bits_Rxp_h(1)
        recuperado(j,i) = bi2de(bits_Rxp(counter-7:counter)','left-msb');
        counter = counter +8;
    end
end

imshow(uint8 (recuperado));