%Miranda Elizabeth Dávila Velarde

%Tarea 5

%Ejercicio 1

numexpediente=47;
sizematrix=numexpediente^2;
load lena512.mat;
figure;
imshow(uint8(lena512));
lenarec=lena512(252:298,318:364);
figure;
imshow(uint8(lenarec));
b=de2bi(lenarec,8,'left-msb'); 
b=b'; 
bits=b(:);   % Vector de bits concatenado

%%
%Reconstrucción
bits_reshape=reshape(bits, 8, sizematrix);
bits_reshape=bits_reshape';
decVal=bi2de(bits_reshape,'left-msb');
lena_reshape=reshape(decVal, size(lenarec));
figure;
imshow(uint8(lena_reshape));

%%

% Pulsos 
Fs=96000;
Ts=1/Fs;
mp=20;

%Unipolar NRZ
pbase = rectwin(mp); 
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = bits;   % Impulse Train
xUNRZ = conv(pbase,s); % Pulse Train

%%
%Polar NRZ
pbase = rectwin(mp);
s1 = bits;
s1(s1==0) = -1; 
s = zeros(1,numel(s1)*mp);
s(1:mp:end) = s1;
xPNRZ = conv(pbase,s);

%%
%Polar RZ
pbase2 = rectwin(mp/2);
s2=bits;
s2(s2==0) = -1; 
s = zeros(1,numel(s2)*mp);
s(1:mp:end) = s2;
xPRZ = conv(pbase2,s);

%%
%Bipolar NRZ
pbase = rectwin(mp);
s2=bits;
lastbit=1;
for i=1:length(bits)
  if bits(i)==1
    s2(i)=lastbit;
    lastbit=-lastbit;
  end
end
s = zeros(1,numel(s2)*mp);
s(1:mp:end) = s2;
xBNRZ = conv(pbase2,s);

%%
%Manchester
pbase3 = rectwin(mp);
pbase3 (mp/2:mp) = -1; 
s3=bits;
s3(s3==1) = -1; 
s3(s3==0) = 1; 
s = zeros(1,numel(s3)*mp);
s(1:mp:end) = s3;
xM = conv(pbase3,s);

%%

%Filtro

%Frecuencia de corte 0.057
f=[0 0.057 0.057 1];
m=[1 1 0 0];
ford=100;
f1=fir2(ford,f,m);
fvtool(f1);

%%
%Frecuencia de corte 0.27
f=[0 0.27 0.27 1];
f2=fir2(ford,f,m);
fvtool(f2);

%%
%Frecuencia de corte 0.47
f=[0 0.47 0.47 1];
f3=fir2(ford,f,m);
fvtool(f3);

%%
%Filtrar los códigos de línea 

%Unipolar NRZ
fxUNRZ1=conv(xUNRZ,f1);
fxUNRZ2=conv(xUNRZ,f2);
fxUNRZ3=conv(xUNRZ,f3);
figure; 
plot(fxUNRZ3(1:mp*16));
title("Convolución del Unipolar NRZ");

%Polar NRZ
fxPNRZ1=conv(xPNRZ,f1);
fxPNRZ2=conv(xPNRZ,f2);
fxPNRZ3=conv(xPNRZ,f3);
figure; 
plot(fxPNRZ3(1:mp*16));
title("Convolución del Polar NRZ");

%Polar RZ
fxPRZ1=conv(xPRZ,f1);
fxPRZ2=conv(xPRZ,f2);
fxPRZ3=conv(xPRZ,f3);
figure; 
plot(fxPRZ3(1:mp*16));
title("Convolución del Polar RZ");


%Bipolar NRZ
fxBNRZ1=conv(xBNRZ,f1);
fxBNRZ2=conv(xBNRZ,f2);
fxBNRZ3=conv(xBNRZ,f3);
figure; 
plot(fxBNRZ3(1:mp*16));
title("Convolución del Bipolar NRZ");


%Manchester
fxM1=conv(xM,f1);
fxM2=conv(xM,f2);
fxM3=conv(xM,f3);
figure; 
plot(fxM3(1:mp*16));
title("Convolución del Manchester");

%%
%Estimador espectral de potencia

figure;
pwelch(fxUNRZ1,500,300,500,Fs,'power');
figure;
pwelch(fxUNRZ2,500,300,500,Fs,'power');
figure;
pwelch(fxUNRZ3,500,300,500,Fs,'power');

%%
figure;
pwelch(fxPNRZ1,500,300,500,Fs,'power');
figure;
pwelch(fxPNRZ2,500,300,500,Fs,'power');
figure;
pwelch(fxPNRZ3,500,300,500,Fs,'power');

%%
figure;
pwelch(fxPRZ1,500,300,500,Fs,'power');
figure;
pwelch(fxPRZ2,500,300,500,Fs,'power');
figure;
pwelch(fxPRZ3,500,300,500,Fs,'power');

%%
figure;
pwelch(fxBNRZ1,500,300,500,Fs,'power');
figure;
pwelch(fxBNRZ2,500,300,500,Fs,'power');
figure;
pwelch(fxBNRZ3,500,300,500,Fs,'power');

%%
figure;
pwelch(fxM1,500,300,500,Fs,'power');
figure;
pwelch(fxM2,500,300,500,Fs,'power');
figure;
pwelch(fxM3,500,300,500,Fs,'power');

%%
%Gráficas stem
%Colores 

figure; 
subplot(2,2,1)
stem(xUNRZ(1:mp*16), 'r');
title("Señal original");
subplot(2,2,2)
stem(fxUNRZ1(48:mp*16+48), 'r');
title("Señal filtrada 0.057");
subplot(2,2,3)
stem(fxUNRZ2(48:mp*16+48), 'r');
title("Señal filtrada 0.27");
subplot(2,2,4)
stem(fxUNRZ3(48:mp*16+48), 'r');
title("Señal filtrada 0.47");
sgtitle("Unipolar NRZ");

%%
figure; 
subplot(2,2,1)
stem(xPNRZ(1:mp*16), 'c');
title("Señal original");
subplot(2,2,2)
stem(fxPNRZ1(48:mp*16+48), 'c');
title("Señal filtrada 0.057");
subplot(2,2,3)
stem(fxPNRZ2(48:mp*16+48), 'c');
title("Señal filtrada 0.27");
subplot(2,2,4)
stem(fxPNRZ3(48:mp*16+48), 'c');
title("Señal filtrada 0.47");
sgtitle("Polar NRZ");

%%
figure; 
subplot(2,2,1)
stem(xPRZ(1:mp*16), 'm');
title("Señal original");
subplot(2,2,2)
stem(fxPRZ1(48:mp*16+48), 'm');
title("Señal filtrada 0.057");
subplot(2,2,3)
stem(fxPRZ2(48:mp*16+48), 'm');
title("Señal filtrada 0.27");
subplot(2,2,4)
stem(fxPRZ3(48:mp*16+48), 'm');
title("Señal filtrada 0.47");
sgtitle("Polar RZ");

%%
figure; 
subplot(2,2,1)
stem(xBNRZ(1:mp*16), 'g');
title("Señal original");
subplot(2,2,2)
stem(fxBNRZ1(48:mp*16+48), 'g');
title("Señal filtrada 0.057");
subplot(2,2,3)
stem(fxBNRZ2(48:mp*16+48), 'g');
title("Señal filtrada 0.27");
subplot(2,2,4)
stem(fxBNRZ3(48:mp*16+48), 'g');
title("Señal filtrada 0.47");
sgtitle("Bipolar RZ");

%%
figure; 
subplot(2,2,1)
stem(xM(1:mp*16), 'k');
title("Señal original");
subplot(2,2,2)
stem(fxM1(48:mp*16+48), 'k');
title("Señal filtrada 0.057");
subplot(2,2,3)
stem(fxM2(48:mp*16+48), 'k');
title("Señal filtrada 0.27");
subplot(2,2,4)
stem(fxM3(48:mp*16+48), 'k');
title("Señal filtrada 0.47");
sgtitle("Manchester");