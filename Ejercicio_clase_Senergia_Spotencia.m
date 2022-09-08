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