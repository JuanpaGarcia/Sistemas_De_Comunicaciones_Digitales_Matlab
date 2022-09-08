% Energy Signal
Fs=10;            % Sampling Frequency
Ts=1/Fs;          % Sampling Time
t = [0:Ts:2*pi];  % Times at which to sample the sine function
Wo= 1;            % 1 rad
sigE = sin(Wo*t); % Original signal, a sine wave
Energy = Ts*sum(sigE.^2); % Otra forma (sigE*sigE')*Ts;
Rxx= Ts*xcorr(sigE,sigE); % Autocorrelation of Energy Signal
% Ts*conv(sigE,fliplr(sigE))
plot(Rxx)
%% Power Signal
sigP= randn(1,1e6);
Power = (sigP*sigP')/numel(sigP);
Rxx = xcorr(sigP,sigP);
plot(Rxx)
histogram(sigP,'Normalization','pdf');

%%
%ejemplo como pa tarea jijij

%Grafique el espectro de un pulso cuadrado
Fs = 100;
Ts = 1/Fs;
t = 0:Ts:1;
x = zeros(1,length(t));
x( (t>=0.4) & (t<0.6) ) = 1;
title('Pulso cuadrado, duración = 0.2s');
% Analice el Ancho de Banda (según criterios)

wvtool(x)%windows view tool
% caida 3 db en freq 0.03125*Fs/2

%%
%SEÑALES DE POTENCIA 
% Power Signal
sigP= randn(1,1e6);
Power = (sigP*sigP')/numel(sigP);
Rxx = xcorr(sigP,sigP);
plot(Rxx)
histogram(sigP,'Normalization','pdf'); 
%Eestimar la densidad de potencia con
pwelch(x,500,300,500,Fs,'power')


