clear all;
close all;
clc 

%Tarea 6
%%
%Low pass filter generation

%F 0.3042
f=[0 0.3042 0.3042 1];
m=[1 1 0 0];
ford=100;
filter_1 = fir2(ford,f,m);
fvtool(filter_1);

%F 0.1542
f=[0  0.1542  0.1542 1];
m=[1 1 0 0];
ford=100;
filter_2 = fir2(ford,f,m);
fvtool(filter_2);

%F 0.042
f=[0  0.042  0.042 1];
m=[1 1 0 0];
ford=100;
filter_3 = fir2(ford,f,m);
fvtool(filter_3);

%%
%load lena512.mat
%tuve que cargar la imagen manual porque no es soportada en mi version de
%matlab
lena512 = imread('lena.tif');
imshow(uint8(lena512)) 
lenarec=lena512(243:284,309:350); 
imshow(uint8(lenarec)) 

b=de2bi(lenarec,8); 
b=b'; 
bits=b(:);   % Bits vector

%%
%Pulse generation

mp = 10;

V_16bit = b(1:42*42*8);
%V_16bit = b(1:16);
%Pulses
%unipolar NRZ
UPNRZ = ones(1,mp);
stem(UPNRZ);

%Polar RZ
PRZ= zeros(1,mp); 
for i = 1:mp
    if i <= mp/2
        PRZ(i) = 1;  
    else
        PRZ(i) = 0; 
    end
end

%Manchester
Manchester = zeros(1,mp); 
for i = 1:mp
    if i <= mp/2
        Manchester(i) = -1;  %load 1
    else
        Manchester(i) = 1; %load -1
    end
end

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
    counter = mp;
end

V_16bit_Unipolar = zeros(1,numel(V_16bit)*mp);
V_16bit_Unipolar(1:mp:end) = V_16bit;


%%
%Signal Filter
Fs = 96000;
%Unipolar NRZ LineCode
Unipolar_NRZ_Sig = conv(UPNRZ ,V_16bit_Unipolar);

% Polar NRZ LineCode
Polar_NRZ_Sig = conv(UPNRZ ,V_16bit_polar);

% Polar RZ LineCode
Polar_RZ_sig = conv(PRZ, V_16bit_polar);

%Bipolar NRZ

help_vetor = V_16bit_Unipolar;
lastbit = 1;
for i=1:length(V_16bit_Unipolar)
  if V_16bit_Unipolar(i)==1
    help_vetor(i) = lastbit;
    lastbit=-lastbit;
  end
end
BP = zeros(1,numel(help_vetor)*mp);
BP(1:mp:end) = help_vetor;
Bipolar_NRZ = conv(UPNRZ,BP);


% Manchester LineCode
Manchester = conv(Manchester, V_16bit_polar);

%Power
pow = (sum(Unipolar_NRZ_Sig.^2))/numel(Unipolar_NRZ_Sig);

%Since we want power equal to 1 watt

% pow_wanted = sqrt(wanted_pow/actualpow)*signal
Unipolar_NRZ_Sig = sqrt(1/pow).*Unipolar_NRZ_Sig;

%check that we obtained the desired power
pow = (sum(Unipolar_NRZ_Sig.^2))/numel(Unipolar_NRZ_Sig);

%change for polar
Polar_NRZ_Sig = sqrt(1/((sum(Polar_NRZ_Sig.^2))/numel(Polar_NRZ_Sig))).*Polar_NRZ_Sig;

%change for polar
Polar_RZ_sig = sqrt(1/((sum(Polar_RZ_sig.^2))/numel(Polar_RZ_sig))).*Polar_RZ_sig;

%change for polar
Bipolar_NRZ = sqrt(1/((sum(Bipolar_NRZ.^2))/numel(Bipolar_NRZ))).*Bipolar_NRZ;

%change for polar
Manchester = sqrt(1/((sum(Manchester.^2))/numel(Manchester))).*Manchester;

%%
%Filter the signals

Signal_filtered_Unipolar_NRZ = conv(Unipolar_NRZ_Sig, filter_1);
tiledlayout(2,2);

nexttile;
plot(Signal_filtered_Unipolar_NRZ)
title('16 bits Filtered signal Unipolar NRZ');

nexttile;
Signal_filtered_Polar_NRZ = conv(Polar_NRZ_Sig, filter_1);
plot(Signal_filtered_Polar_NRZ)
title('16 bits Filtered signal Polar NRZ');

nexttile;
Signal_filtered_Bipolar_NRZ = conv(Bipolar_NRZ, filter_1);
plot(Signal_filtered_Bipolar_NRZ)

Signal_filtered_Manchester = conv(Manchester, filter_3);
plot(Signal_filtered_Manchester)
title('16 bits Filtered signal Manchester');

nexttile;
Signal_filtered_Polar_RZ = conv(Polar_RZ_sig, filter_1);
plot(Signal_filtered_Polar_RZ)
title('16 bits Filtered signal AMI RZ');
%%
%Power spectral density



%%
%bit recovery 
%UPNRZ
delay_signal = ford/2;
start_recovery_count = delay_signal + mp/2;
UPNRZ_recovery = Signal_filtered_Unipolar_NRZ(start_recovery_count:mp:end);
Decition_treshold_UPNRZ = 0.7;
%view the samples in a complex plane.
scatterplot(UPNRZ_recovery);
%%
%PNRZ
delay_signal = ford/2;
start_recovery_count = delay_signal + mp/2;
PNRZ_recovery = Signal_filtered_Polar_NRZ(start_recovery_count:mp:end);
Decition_treshold_PNRZ = 0.5;
%view the samples in a complex plane.
scatterplot(PNRZ_recovery);

%%
%PRZ
delay_signal = ford/2;
start_recovery_count = round(delay_signal + mp/4);
PRZ_recovery = Signal_filtered_Polar_RZ(start_recovery_count:mp:end);
Decition_treshold_PRZ = 0.5;
%view the samples in a complex plane.
scatterplot(PRZ_recovery);

%%
%Manchester
delay_signal = ford/2;
start_recovery_count = round(delay_signal + mp/4);
Manchester_recovery_y1 = Signal_filtered_Manchester(start_recovery_count:mp:end);
start_recovery_count = round(delay_signal + 3*mp/4);
Manchester_recovery_y2 = Signal_filtered_Manchester(start_recovery_count:mp:end);
Decition_treshold_PRZ = 0.5;

%view the samples in a complex plane.
%scatterplot(Manchester_recovery);

%%
%convert simbols to bits
UPNRZ_recovery_bits = zeros(1,numel(V_16bit)); 
UPNRZ_recovery_bits((UPNRZ_recovery >= Decition_treshold_UPNRZ)) = 1; 

PNRZ_recovery_bits = zeros(1,numel(V_16bit));
PNRZ_recovery_bits((PNRZ_recovery > Decition_treshold_UPNRZ)) = 1; 


PRZ_recovery_bits = zeros(1,numel(V_16bit));
PRZ_recovery_bits((PRZ_recovery > Decition_treshold_UPNRZ)) = 1; 
PRZ_recovery_bits((PNRZ_recovery <= -Decition_treshold_UPNRZ)) = 0;


Manchester_recovery_bits = zeros(1,numel(V_16bit));
Manchester_recovery_bits = ( sign(Manchester_recovery_y2 - Manchester_recovery_y1) +1 )/2;
 


%%
%checar bits recibidos con error
%UPNRZ
bits_error_UPNRZ = xor(V_16bit, UPNRZ_recovery_bits(1: numel(V_16bit)) );
error_UPNRZ = sum(bits_error_UPNRZ);
Bit_error_rate_UPNRZ = (error_UPNRZ/numel(V_16bit)) * 100;

%PNRZ
bits_error_PNRZ = xor(V_16bit, PNRZ_recovery_bits(1: numel(V_16bit)) );
error_PNRZ = sum(bits_error_PNRZ);
Bit_error_rate_PNRZ = (error_PNRZ/numel(V_16bit)) * 100;

%PRZ
bits_error_PRZ = xor(V_16bit, PRZ_recovery_bits(1: numel(V_16bit)) );
error_PRZ = sum(bits_error_PRZ);
Bit_error_rate_PRZ = (error_PRZ/numel(V_16bit)) * 100;

%Manchester
bits_error_Manchester = xor(V_16bit, Manchester_recovery_bits(1: numel(V_16bit)) );
error_Manchester = sum(bits_error_PRZ);
Bit_error_rate_Manchester = (error_Manchester/numel(V_16bit)) * 100;


%%
%Recuperar imagen
%UPNRZ
pixels = 42;
recuperado = zeros(pixels,pixels,'uint32'); %allocate memory

%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : pixels
    for j = 1: pixels
        recuperado(j,i) = bi2de(UPNRZ_recovery_bits(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end
imshow(uint8(recuperado)) %recovered image

%%
%PNRZ
recuperado = zeros(42,42,'uint32'); %allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : pixels
    for j = 1: pixels
        recuperado(j,i) = bi2de(PNRZ_recovery_bits(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end
imshow(uint8(recuperado)) %recovered image

%%
%PRZ
recuperado = zeros(42,42,'uint32'); %allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : pixels
    for j = 1: pixels
        recuperado(j,i) = bi2de(PRZ_recovery_bits(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end
imshow(uint8(recuperado)) %recovered image

%%
%Manchester
recuperado = zeros(42,42,'uint32'); %allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : pixels
    for j = 1: pixels
        recuperado(j,i) = bi2de(Manchester_recovery_bits(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end
imshow(uint8(recuperado)) %recovered image

%%
%Bipolar
recuperado = zeros(42,42,'uint32'); %allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : pixels
    for j = 1: pixels
        recuperado(j,i) = bi2de(Manchester_recovery_bits(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end
imshow(uint8(recuperado)) %recovered image

