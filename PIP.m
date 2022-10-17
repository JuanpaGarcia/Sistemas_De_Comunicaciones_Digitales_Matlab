
%Lena
numexpediente=42;
plots = 0;
sizematrix=numexpediente^2;
load lena512.mat;
lenarec=lena512(252:298,318:364);
b=de2bi(lenarec,8,'left-msb'); 
b=b'; 
bits=b(:);   % Vector de bits concatenado

%%

% Pulsos 
Fs=96000;
Ts=1/Fs;
mp=10;
baudrate=Fs/mp;
potencia_deseada=sqrt(1);

%Polar NRZ
pbasePNRZ = rectwin(mp);
s1 = bits;
s1(s1==0) = -1; 
s = zeros(1,numel(s1)*mp);
s(1:mp:end) = s1;
xPNRZ = conv(pbasePNRZ,s);
p_xPNRZ= (sum(xPNRZ.^2))/(numel(xPNRZ));
xPNRZ = (xPNRZ/sqrt(p_xPNRZ))*potencia_deseada;
p_xPNRZ= (sum(xPNRZ.^2))/(numel(xPNRZ));


%%
%Generación de señal half-sine

n = 0:mp-1; 
w0=pi/(mp);
Half_sine= sin(w0*n);

xPNRZ_HS = conv(Half_sine,s);
p_xPNRZ_HS= (sum(xPNRZ_HS.^2))/(numel(xPNRZ_HS));
xPNRZ_HS = (xPNRZ_HS/sqrt(p_xPNRZ_HS))*potencia_deseada;
p_xPNRZ_HS= (sum(xPNRZ_HS.^2))/(numel(xPNRZ_HS));

%%
%Filtro

%Frecuencia de corte 0.4
f=[0 0.6 0.6 1];
m=[1 1 0 0];
ford=60;
filter_delay=ford/2;
f1=fir2(ford,f,m);
if plots ==1
    fvtool(f1);
end


%%
%Canal
fxPNRZ=conv(xPNRZ,f1);
fxPNRZ_HS=conv(xPNRZ_HS,f1);

%%
%Polar NRZ

%PNRZ
p_fxPNRZ= (sum(fxPNRZ.^2))/(numel(fxPNRZ));
fxPNRZ = (fxPNRZ/sqrt(p_fxPNRZ))*potencia_deseada;
p_fxPNRZ= (sum(fxPNRZ.^2))/(numel(fxPNRZ));

%PNRZ HS
p_fxPNRZ_HS= (sum(fxPNRZ_HS.^2))/(numel(fxPNRZ_HS));
fxPNRZ_HS = (fxPNRZ_HS/sqrt(p_fxPNRZ_HS))*potencia_deseada;
p_fxPNRZ_HS= (sum(fxPNRZ_HS.^2))/(numel(fxPNRZ_HS));

%%

%Ruido

PNoise=1;

%PNoise=0:3:30;

%Ruido
Noise_PNRZ=sqrt(PNoise)*randn(1,numel(fxPNRZ));
Noise_PNRZ_HS=sqrt(PNoise)*randn(1,numel(fxPNRZ_HS));

%Potencia del ruido
PNoise_PNRZ=var(Noise_PNRZ);
PNoise_PNRZ_HS=var(Noise_PNRZ_HS);

%SNR
SNR_dB_PNRZ=10*log10(p_fxPNRZ/PNoise);
SNR_dB_PNRZ_HS=10*log10(p_fxPNRZ_HS/PNoise);

%Añadir ruido
fxPNRZ_AWGN=fxPNRZ+Noise_PNRZ;
fxPNRZ_HS_AWGN=fxPNRZ_HS+Noise_PNRZ_HS;

%%
%Match Filter

filter_recovery_delay = mp/2;

%PNRZ
pbasePNRZ_receptor=fliplr(pbasePNRZ);

recover_PNRZ=conv(fxPNRZ_AWGN,pbasePNRZ_receptor);

%PNRZ_HS
pbasePNRZ_HS_receptor=fliplr(Half_sine);

recover_PNRZ_HS=conv(fxPNRZ_HS_AWGN,pbasePNRZ_HS_receptor);

%%
%Estimador espectral de potencia
if plots ==1
    figure;
    pwelch(recover_PNRZ,500,300,500,Fs,'power');
    figure;
    pwelch(recover_PNRZ_HS,500,300,500,Fs,'power');
end


%%
%Diagrama de ojo
if plots ==1
    eyediagram(recover_PNRZ,2*mp);
    eyediagram(recover_PNRZ_HS,2*mp);
end


%%
%Muestreo

%PNRZ
start=filter_delay+(mp/2)+filter_recovery_delay;

MfxPNRZ=recover_PNRZ(start:mp:end);

if plots ==1
    scatterplot(MfxPNRZ);
end


umbral_PolarNRZ=0;

bits_Rx_PNRZ=zeros(1,numel(MfxPNRZ));

bits_Rx_PNRZ(MfxPNRZ>=umbral_PolarNRZ)=1;

bits_Rx_PNRZ(MfxPNRZ<umbral_PolarNRZ)=0;

bits_Rx_PNRZ=bits_Rx_PNRZ(1:numel(bits));

bits_Rx_PNRZ=bits_Rx_PNRZ';

bits_Rx_PNRZ=bits_Rx_PNRZ(:);

bits_error=sum(xor(bits,bits_Rx_PNRZ(1:numel(bits))));

BER_PNRZ=(bits_error/numel(bits))*100;

%Recuperación
bits_reshape=reshape(bits_Rx_PNRZ, 8, sizematrix);

bits_reshape=bits_reshape';

decVal=bi2de(bits_reshape,'left-msb');

lena_reshape=reshape(decVal, size(lenarec));

if plots ==1
    figure;
    imshow(uint8(lena_reshape));
end


%%
%PNRZ_HS

start=filter_delay+(mp/2)+filter_recovery_delay;

MfxPNRZ_HS=recover_PNRZ_HS(start:mp:end);

if plots ==1
    scatterplot(MfxPNRZ_HS);
end

umbral_PolarNRZ_HS=0;

bits_Rx_PNRZ_HS=zeros(1,numel(MfxPNRZ_HS));

bits_Rx_PNRZ_HS(MfxPNRZ_HS>=umbral_PolarNRZ_HS)=1;

bits_Rx_PNRZ_HS(MfxPNRZ_HS<umbral_PolarNRZ_HS)=0;

bits_Rx_PNRZ_HS=bits_Rx_PNRZ_HS(1:numel(bits));

bits_Rx_PNRZ_HS=bits_Rx_PNRZ_HS';

bits_Rx_PNRZ_HS=bits_Rx_PNRZ_HS(:);

bits_error=sum(xor(bits,bits_Rx_PNRZ_HS(1:numel(bits))));

BER_PNRZ_HS=(bits_error/numel(bits))*100;

%Recuperación
bits_reshape=reshape(bits_Rx_PNRZ_HS, 8, sizematrix);

bits_reshape=bits_reshape';

decVal=bi2de(bits_reshape,'left-msb');

lena_reshape=reshape(decVal, size(lenarec));

if plots ==1
    figure;
    imshow(uint8(lena_reshape));
end




























































































%%
%Cambiando PNoise

n=0:3:30;
PNoise=10.^(n./10);

%%
%Ruido
i=11;
Noise_PNRZ=sqrt(PNoise(i)).*randn(1,numel(fxPNRZ));
Noise_PNRZ_HS=sqrt(PNoise(i)).*randn(1,numel(fxPNRZ_HS));

%Potencia del ruido
PNoise_PNRZ=var(Noise_PNRZ);
PNoise_PNRZ_HS=var(Noise_PNRZ_HS);

%SNR
SNR_dB_PNRZ=10*log10(p_fxPNRZ/PNoise(i));
SNR_dB_PNRZ_HS=10*log10(p_fxPNRZ_HS/PNoise(i));

%Añadir ruido
fxPNRZ_AWGN=fxPNRZ+Noise_PNRZ;
fxPNRZ_HS_AWGN=fxPNRZ_HS+Noise_PNRZ_HS;

%%
%Match Filter

filter_recovery_delay = mp/2;

%PNRZ
pbasePNRZ_receptor=fliplr(pbasePNRZ);

recover_PNRZ=conv(fxPNRZ_AWGN,pbasePNRZ_receptor);

%PNRZ_HS
pbasePNRZ_HS_receptor=fliplr(Half_sine);

recover_PNRZ_HS=conv(fxPNRZ_HS_AWGN,pbasePNRZ_HS_receptor);

%%
%Estimador espectral de potencia

figure;
pwelch(recover_PNRZ,500,300,500,Fs,'power');
figure;
pwelch(recover_PNRZ_HS,500,300,500,Fs,'power');

%%
%Diagrama de ojo

eyediagram(recover_PNRZ,2*mp);
eyediagram(recover_PNRZ_HS,2*mp);

%%
%Muestreo

%PNRZ
start=filter_delay+(mp/2)+filter_recovery_delay;

MfxPNRZ=recover_PNRZ(start:mp:end);

scatterplot(MfxPNRZ);

umbral_PolarNRZ=0;

bits_Rx_PNRZ=zeros(1,numel(MfxPNRZ));

bits_Rx_PNRZ(MfxPNRZ>=umbral_PolarNRZ)=1;

bits_Rx_PNRZ(MfxPNRZ<umbral_PolarNRZ)=0;

bits_Rx_PNRZ=bits_Rx_PNRZ(1:numel(bits));

bits_Rx_PNRZ=bits_Rx_PNRZ';

bits_Rx_PNRZ=bits_Rx_PNRZ(:);

bits_error=sum(xor(bits,bits_Rx_PNRZ(1:numel(bits))));

BER_PNRZ=(bits_error/numel(bits))*100;

%Recuperación
bits_reshape=reshape(bits_Rx_PNRZ, 8, sizematrix);

bits_reshape=bits_reshape';

decVal=bi2de(bits_reshape,'left-msb');

lena_reshape=reshape(decVal, size(lenarec));

figure;

imshow(uint8(lena_reshape));

%%
%PNRZ_HS

start=filter_delay+(mp/2)+filter_recovery_delay;

MfxPNRZ_HS=recover_PNRZ_HS(start:mp:end);

scatterplot(MfxPNRZ_HS);

umbral_PolarNRZ_HS=0;

bits_Rx_PNRZ_HS=zeros(1,numel(MfxPNRZ_HS));

bits_Rx_PNRZ_HS(MfxPNRZ_HS>=umbral_PolarNRZ_HS)=1;

bits_Rx_PNRZ_HS(MfxPNRZ_HS<umbral_PolarNRZ_HS)=0;

bits_Rx_PNRZ_HS=bits_Rx_PNRZ_HS(1:numel(bits));

bits_Rx_PNRZ_HS=bits_Rx_PNRZ_HS';

bits_Rx_PNRZ_HS=bits_Rx_PNRZ_HS(:);

bits_error=sum(xor(bits,bits_Rx_PNRZ_HS(1:numel(bits))));

BER_PNRZ_HS=(bits_error/numel(bits))*100;

%Recuperación
bits_reshape=reshape(bits_Rx_PNRZ_HS, 8, sizematrix);

bits_reshape=bits_reshape';

decVal=bi2de(bits_reshape,'left-msb');

lena_reshape=reshape(decVal, size(lenarec));

figure;

imshow(uint8(lena_reshape));

