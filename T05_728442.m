clear all;
close all;
clc 
%%
mp = 11; % muestras por pulso 
pbase = triang(mp); % ejemplo con pulso base triangular de mp muestras
b = [1 0 1 0 1 1 0 0 1 0 1 0 1 1 0 0]; 
s = zeros(1,numel(b)*mp);
s(1:mp:end) = b *2-1;
stem(s)
%%
x = conv(s,pbase);  % Verify dimentions
plot(x) 
hold on 
stem(x,'LineStyle','none') 

%%
%Ejercicio1
%load lena512.mat
%tuve que cargar la imagen manual porque no es soportada en mi version de
%matlab
lena512 = imread('lena.tif');
imshow(uint8(lena512)) 
lenarec=lena512(243:284,309:350); 
imshow(uint8(lenarec)) 

%%
b=de2bi(lenarec,8); 
b=b'; 
bits=b(:);   % Bits vector 
%%
%Image recovery

recuperado = zeros(42,42,'uint32'); %allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : 42
    for j = 1: 42
        recuperado(j,i) = bi2de(b(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end
imshow(uint8(recuperado)) %recovered image

%%
%Pulse generation

mp = 20;

V_16bit = b(1:16);

%Pulses
%unipolar NRZ
UPNRZ = ones(1,mp);
stem(UPNRZ)
wvtool(UPNRZ) 
%%
%Polar RZ
PRZ= zeros(1,mp); 
for i = 1:mp
    if i <= mp/2
        PRZ(i) = 1;  
    else
        PRZ(i) = 0; 
    end
end
stem(PRZ)
wvtool(PRZ) 

%%
%Manchester
Manchester = zeros(1,mp); 
for i = 1:mp
    if i <= mp/2
        Manchester(i) = 1;  %load 1
    else
        Manchester(i) = -1; %load -1
    end
end
stem(Manchester)
wvtool(Manchester) 




%%
%Signal Generation

V_16bit_polar = zeros(1,numel(V_16bit)*mp);
%V_16bit_polar(1:mp:end) = V_16bit*2-1; 
counter = 0;
for i= 0 : numel(V_16bit)-1
    if V_16bit(i+1) == 0
        value = -1;
    else
        value = 1;
    end
    V_16bit_polar(counter*i+1) = value;
    counter =20;
end

V_16bit_Unipolar = zeros(1,numel(V_16bit)*mp);
V_16bit_Unipolar(1:mp:end) = V_16bit;
%%
Fs = 9600;
%Unipolar NRZ LineCode
Unipolar_NRZ_Sig = conv(UPNRZ ,V_16bit_Unipolar);
plot(Unipolar_NRZ_Sig)
title('Signal on Unipolar NRZ line code');
figure;

pwelch(Unipolar_NRZ_Sig,[],[],[],Fs,'power');  % PSD of Unipolar NRZ
title('Spectral power density');
%%
% Polar NRZ LineCode
Polar_NRZ_Sig = conv(UPNRZ ,V_16bit_polar);
plot(Polar_NRZ_Sig)
title('Signal on Polar NRZ line code');

pwelch(Polar_NRZ_Sig,[],[],[],Fs,'power');  % PSD of Unipolar NRZ
title('Power spectral density');
%%
% Polar RZ LineCode
Polar_RZ_sig = conv(PRZ, V_16bit_polar);
plot(Polar_RZ_sig)
title('Signal on Polar RZ line code');

pwelch(Polar_RZ_sig,[],[],[],Fs,'power');  % PSD of Unipolar NRZ
title('Power spectral density');
%%
% Bipolar NRZ LineCode
Bipolar_NRZ = conv(UPNRZ, V_16bit_polar);
plot(Bipolar_NRZ);
title('Signal on Bipolar NRZ line code');

pwelch(Bipolar_NRZ,[],[],[],Fs,'power');  % PSD of Unipolar NRZ
title('Power spectral density');

%%
% Manchester LineCode
Manchester = conv(Manchester, V_16bit_polar);
plot(Manchester);
title('Signal on manchester line code');

pwelch(Bipolar_NRZ,[],[],[],Fs,'power');  % PSD of Unipolar NRZ
title('Power spectral density');

%%
%Canal comunicacion
fv=   [0 0.2 0.2 1];  % Vector de Frecuencias 
mv= [1 1 0 0];        % Vector de Magnitudes 
fordv=100;             % Orden del Filtro 
f1v = fir2(fordv,fv,mv);  % Coeficientes del Filtro usando FIR2( ) 
fvtool(f1v)               % Herramienta para analizar el Filtro 
%%
%Filter design
%Filtro

%F 0.052
f=[0 0.052 0.052 1];
m=[1 1 0 0];
ford=100;
filter_1 = fir2(ford,f,m);
fvtool(filter_1);

%%
%Frecuencia de corte 0.22
f=[0 0.22 0.22 1];
filter_2=fir2(ford,f,m);
fvtool(filter_2);

%%
%Frecuencia de corte 0.42
f=[0 0.42 0.42 1];
filter_3=fir2(ford,f,m);
fvtool(filter_3);

%%
%Line code filtered

%Unipolar 
Unipolar_f1 = conv(Unipolar_NRZ_Sig,filter_1);
Unipolar_f2 = conv(Unipolar_NRZ_Sig,filter_2);
Unipolar_f3 = conv(Unipolar_NRZ_Sig,filter_3);

plot(Unipolar_f3);
title('Filtered signal Unipolar');

%%
%Line code filtered

%polar 
polar_f1 = conv(Polar_NRZ_Sig,filter_1);
polar_f2 = conv(Polar_NRZ_Sig,filter_2);
polar_f3 = conv(Polar_NRZ_Sig,filter_3);

plot(polar_f3);
title('Filtered signal Polar');
%%
%Line code filtered

%polar RZ

polarz_f1 = conv(Polar_RZ_sig,filter_1);
polarz_f2 = conv(Polar_RZ_sig,filter_2);
polarz_f3 = conv(Polar_RZ_sig,filter_3);
plot(polarz_f3);
title('Filtered signal Polar RZ');
%%
%Line code filtered
%bipolar nRZ

bpolanrz_f1 = conv(Bipolar_NRZ,filter_1);
bpolanrz_f2 = conv(Bipolar_NRZ,filter_2);
bpolanrz_f3 = conv(Bipolar_NRZ,filter_3);
plot(bpolanrz_f3);
title('Filtered signal Biolar NRZ');

%%
%Line code filtered
%Manchester
M_f1 = conv(Manchester,filter_1);
M_f2 = conv(Manchester,filter_2);
M_f3 = conv(Manchester,filter_3);
plot(M_f3);
title('Filtered signal Machester');

%%
%Spectrums
%Unipolar 
tiledlayout(2,2);
nexttile;
pwelch(Unipolar_f1,[],[],[],Fs,'power');
title('Unipolar_f1')
nexttile;
pwelch(Unipolar_f1,[],[],[],Fs,'power');
title('Unipolar_f2')
nexttile([1 2])
pwelch(Unipolar_f1,[],[],[],Fs,'power');
title('Unipolar_f3')

%%
%Spectrums
%polar 
tiledlayout(2,2);
nexttile;
pwelch(polar_f1,[],[],[],Fs,'power');
title('polar_f1')
nexttile;
pwelch(polar_f2,[],[],[],Fs,'power');
title('polar_f2')
nexttile([1 2])
pwelch(polar_f3,[],[],[],Fs,'power');
title('polar_f3')

%%
%Spectrums
%polar RZ
tiledlayout(2,2);
nexttile;
pwelch(polarz_f1,[],[],[],Fs,'power');
title('polarz_f1')
nexttile;
pwelch(polarz_f2,[],[],[],Fs,'power');
title('polarz_f2')
nexttile([1 2])
pwelch(polarz_f3,[],[],[],Fs,'power');
title('polarz_f3')

%%
%Spectrums
%bipolar NRZ
tiledlayout(2,2);
nexttile;
pwelch(bpolanrz_f1,[],[],[],Fs,'power');
title('bpolanrz_f1')
nexttile;
pwelch(bpolanrz_f2,[],[],[],Fs,'power');
title('bpolanrz_f2')
nexttile([1 2])
pwelch(bpolanrz_f3,[],[],[],Fs,'power');
title('bpolanrz_f3')

%%
%Spectrums
%Manchester
tiledlayout(2,2);
nexttile;
pwelch(M_f1,[],[],[],Fs,'power');
title('M_f1')
nexttile;
pwelch(M_f2,[],[],[],Fs,'power');
title('M_f2')
nexttile([1 2])
pwelch(M_f3,[],[],[],Fs,'power');
title('M_f3')

%%
%discrete time graphs
%Unipolar NRZ
delay = 48;
tiledlayout(2,2);
nexttile;
stem(Unipolar_NRZ_Sig(1:(16*mp)));
title('Señal Original');
nexttile;
stem(Unipolar_f1(1:(16*mp + delay)));
title('Señal filtrada 1');
nexttile;
stem(Unipolar_f2(1:(16*mp + delay)));
title('Señal filtrada 2');

nexttile;
stem(Unipolar_f3(1:(16*mp + delay)));
title('Señal filtrada 3');
%%
%Polar NRZ
delay = 48;
tiledlayout(2,2);
nexttile;
stem(Polar_NRZ_Sig(1:(16*mp)));
title('Señal Original');
nexttile;
stem(polar_f1(1:(16*mp + delay)));
title('Señal filtrada 1');
nexttile;
stem(polar_f2(1:(16*mp + delay)));
title('Señal filtrada 2');

nexttile;
stem(polar_f3(1:(16*mp + delay)));
title('Señal filtrada 3');

%%
%Biolar NRZ
delay = 48;
tiledlayout(2,2);
nexttile;
stem(Bipolar_NRZ(1:(16*mp)));
title('Señal Original');
nexttile;
stem(bpolanrz_f1(1:(16*mp + delay)));
title('Señal filtrada 1');
nexttile;
stem(bpolanrz_f2(1:(16*mp + delay)));
title('Señal filtrada 2');

nexttile;
stem(bpolanrz_f3(1:(16*mp + delay)));
title('Señal filtrada 3');


%%
%Polar RZ
delay = 48;
tiledlayout(2,2);
nexttile;
stem(Polar_RZ_sig(1:(16*mp)));
title('Señal Original');
nexttile;
stem(polarz_f1(1:(16*mp + delay)));
title('Señal filtrada 1');
nexttile;
stem(polarz_f1(1:(16*mp + delay)));
title('Señal filtrada 2');

nexttile;
stem(polarz_f1(1:(16*mp + delay)));
title('Señal filtrada 3');
%%
%Manchester
delay = 48;
tiledlayout(2,2);
nexttile;
stem(Manchester(1:(16*mp)));
title('Señal Original');
nexttile;
stem(M_f1(1:(16*mp + delay)));
title('Señal filtrada 1');
nexttile;
stem(M_f2(1:(16*mp + delay)));
title('Señal filtrada 2');
nexttile;
stem(M_f3(1:(16*mp + delay)));
title('Señal filtrada 3');

