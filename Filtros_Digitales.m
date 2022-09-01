Fs=2000; Ts = 1/Fs; t = 0:Ts:1; % intervalo de muestreo y vector de tiempo
x = cos(2*pi*100*t)+cos(2*pi*500*t)+cos(2*pi*800*t);
f = [0 1/3 1/3 1];
m = [1 1 0 0]; o = 99;
b = fir2(o,f,m);
fvtool(b);  % Espectro en frecuencia del filtro.
y = conv(x,b); % filtrado de x con filtro b o usando:
%y = filter(b,1,x]); % filtrado de x con filtro b, %función filter
%y = fftfilt(b,x);   % A more efficient FIR filtering for large operands