clear all;
%******                             TRANSMITION                   ********%
%********************************** FASE I *******************************%
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
%%
%******************************** Fase II ********************************%
% PREAMBLE CREATION AND LENA LOAD 
clear all;
preamble= [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0]' ; % poner 56
SFD = [1 0 1 0 1 0 1 1]' ;
DSA = de2bi(uint8('Practica II FASE II: 726821 y 728442'),8,'left-msb'); 
DSA = reshape(DSA',numel(DSA),1); 
load lena512.mat; img = uint8(lena512); 
img = img(248:247+41,245:244+42,1); % Image size= 4(UDE1) x 4(UDE2) pixels 
imshow(img);  % Where UDE =Último Dígito Expediente 
size_img = de2bi(size(img),16,'left-msb');
header= [size_img(1,:) size_img(2,:)]'; 
b = de2bi(img,8,'left-msb');
b = b'; 
bits = b(:);  
header_bits = [preamble; SFD; DSA];
trash_bits = [preamble; SFD; DSA; header];
bits2Tx = [preamble; SFD; DSA; header; bits];
%%
% PULSE DESIGN 
%************************** 4KHz ***********************************%
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

% PNRZ LINE CODING
s1 = int8(bits2Tx);
s1(s1==0) = -1;  % Convert “0” to “-1”
s = zeros(1,numel(s1)*mp);
s(1:mp:end) = s1;

%PULSE TRAIN 
xPNRZ = conv(pbase,s); % Line code 
powe_polar = sum(xPNRZ.^2)/numel(xPNRZ);
norm_polar = sqrt(1/powe_polar)*xPNRZ; % 1 Watt Polar NRZ

pulse_time_transmition = mp/Fs;
pulseTrain_time_transmition = numel(bits2Tx)/Rb;

eyediagram(norm_polar, mp*2);
%%
% TRANSMITION 
silence = zeros(1,(0.5*Fs));
TxPulseTrain = [silence, norm_polar];

soundsc( TxPulseTrain, Fs );

%%
% PULSE DESIGN 
%************************** 12KHz ***********************************%
Fs      =   96e3;              % Samples per second  
Ts      =   1/Fs;              % Sampling period 
beta    =   0.25;              % Roll-off factor 
B       =   12000;              % Bandwidth available 
Rb      =   2*B/(1+beta);      % Bit rate = Baud rate 
mp      =   ceil(Fs/Rb)        % samples per pulse 
Rb      =   Fs/mp;             % Recompute bit rate 
Tp      =   1/Rb;              % Symbol period 
B       =   (Rb*(1+beta)/2)    % Bandwidth consumed 
D       =   10;                % Time duration in terms of Tp 
type    =   'srrc';            % Shape pulse: Square Root Rise Cosine 
E       =   Tp;                % Energy 
[pbase ~] = rcpulse(beta, D, Tp, Ts, type, E);    % Pulse Generation 

% PNRZ LINE CODING
s1 = int8(bits2Tx);
s1(s1==0) = -1;  % Convert “0” to “-1”
s = zeros(1,numel(s1)*mp);
s(1:mp:end) = s1;

%PULSE TRAIN 
xPNRZ = conv(pbase,s); % Line code 
powe_polar = sum(xPNRZ.^2)/numel(xPNRZ);
norm_polar = sqrt(1/powe_polar)*xPNRZ; % 1 Watt Polar NRZ

pulse_time_transmition = mp/Fs;
pulseTrain_time_transmition = numel(bits2Tx)/Rb;

eyediagram(norm_polar, mp*2);
%%
% TRANSMITION
silence = zeros(1,(0.5*Fs));
TxPulseTrain = [silence, norm_polar];

soundsc( TxPulseTrain, Fs );
%%
% PULSE DESIGN 
%************************** 20KHz ***********************************%
Fs      =   96e3;              % Samples per second  
Ts      =   1/Fs;              % Sampling period 
beta    =   0.25;              % Roll-off factor 
Rb      =   19200*2;           % Doubling Rb from the 12KHz BW design 
mp      =   ceil(Fs/Rb)        % samples per pulse 
Rb      =   Fs/mp;             % Recompute bit rate 
Tp      =   1/Rb;              % Symbol period 
B       =   (Rb*(1+beta)/2)    % Bandwidth consumed 
D       =   10;                % Time duration in terms of Tp 
type    =   'srrc';            % Shape pulse: Square Root Rise Cosine 
E       =   Tp;                % Energy 
[pbase ~] = rcpulse(beta, D, Tp, Ts, type, E);    % Pulse Generation 

% PNRZ LINE CODING
s1 = int8(bits2Tx);
s1(s1==0) = -1;  % Convert “0” to “-1”
s = zeros(1,numel(s1)*mp);
s(1:mp:end) = s1;

%PULSE TRAIN 
xPNRZ = conv(pbase,s);
powe_polar = sum(xPNRZ.^2)/numel(xPNRZ);
norm_polar = sqrt(1/powe_polar)*xPNRZ; % 1 Watt Polar NRZ

pulse_time_transmition = mp/Fs;
pulseTrain_time_transmition = numel(bits2Tx)/Rb;

eyediagram(norm_polar(1:3000), mp*2);
%%
% TRANSMITION
silence = zeros(1,(0.5*Fs));
TxPulseTrain = [silence, norm_polar];

soundsc( TxPulseTrain, Fs );