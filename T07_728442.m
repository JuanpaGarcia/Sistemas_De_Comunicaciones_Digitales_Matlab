clear all;
close all;
clc 

%Tarea 7

%%
mp = 10;
Fs = 96000;

lena512 = imread('lena.tif');
imshow(uint8(lena512)) 
lenarec=lena512(243:284,309:350); 
imshow(uint8(lenarec)) 

b=de2bi(lenarec,8); 
b=b'; 
bits=b(:);   % Bits vector

pixels = 42;
V_bit = b(1:pixels*pixels*8); %8 because 8bit pixel

V_bit_polar = zeros(1,numel(V_bit)*mp);

counter = 0;
for i= 0 : numel(V_bit)-1
    if V_bit(i+1) == 0
        value = -1;
    else
        value = 1;
    end
    V_bit_polar(counter*i+1) = value;
    counter = mp;
end

%%
n = 0:mp-1;
w = pi/mp;
hs = sin(w*n);
pr = ones(1,mp);

% Polar NRZ LineCode rectangular
Polar_NRZ_sig_rec = conv(pr ,V_bit_polar);

% Polar NRZ LineCode halfsine
Polar_NRZ_sig_hs = conv(hs ,V_bit_polar);

%change power to 1
Polar_NRZ_sig_rec = sqrt(1/((sum(Polar_NRZ_sig_rec.^2))/numel(Polar_NRZ_sig_rec))).*Polar_NRZ_sig_rec;

%change power to 1
Polar_NRZ_sig_hs = sqrt(1/((sum(Polar_NRZ_sig_hs.^2))/numel(Polar_NRZ_sig_hs))).*Polar_NRZ_sig_hs;
%%
plot(Polar_NRZ_sig_rec);
title('Señal PNRZ rectangular');

figure();
plot(Polar_NRZ_sig_hs);
title('Señal PNRZ half sine');

eyediagram(Polar_NRZ_sig_rec,2*mp);
title('Eyediagrma rectangular');

eyediagram(Polar_NRZ_sig_hs,2*mp);
title('Eyediagrma half sine');

%%
%Ejercicio 3
%Low pass filter generation, communications channel
%F 0.3042
f=[0 0.4 0.4 1];
m=[1 1 0 0];
ford=60;
filter_1 = fir2(ford,f,m);
%%
%transmit the signal in the communication channel
Signal_filtered_Polar_NRZ_rec = conv(Polar_NRZ_sig_rec, filter_1);
plot(Signal_filtered_Polar_NRZ_rec)
title('Filtered signal Polar NRZ rectangular');


Signal_filtered_Polar_NRZ_hs = conv(Polar_NRZ_sig_hs, filter_1);
figure();
plot(Signal_filtered_Polar_NRZ_hs)
title('Filtered signal Polar NRZ half sine');

eyediagram(Signal_filtered_Polar_NRZ_rec,2*mp);
title('Eyediagrma rectangular');

eyediagram(Signal_filtered_Polar_NRZ_hs,2*mp);
title('Eyediagrma half sine');

%%
%obtener ancho de banda 
% fs/mp
pwelch(Polar_NRZ_sig_hs);
title('Line coding bandwidth');

%%
%change power to 1
Signal_filtered_Polar_NRZ_hs = sqrt(1/((sum(Signal_filtered_Polar_NRZ_hs.^2))/numel(Signal_filtered_Polar_NRZ_hs))).*Signal_filtered_Polar_NRZ_hs;
plot(Signal_filtered_Polar_NRZ_hs)
title('Normalized Filtered signal Polar NRZ half sine');

%check power pow = (sum(Unipolar_NRZ_Sig.^2))/numel(Unipolar_NRZ_Sig);

%change power to 1
Signal_filtered_Polar_NRZ_rec = sqrt(1/((sum(Signal_filtered_Polar_NRZ_rec.^2))/numel(Signal_filtered_Polar_NRZ_rec))).*Signal_filtered_Polar_NRZ_rec;
figure();
plot(Signal_filtered_Polar_NRZ_rec)
title('Normalized Filtered signal Polar NRZ rectangular');

%%
%Ejercicio 5
%apply match filter

match_filtered_PNRZ_rec = conv(Signal_filtered_Polar_NRZ_rec, fliplr(pr));

match_filtered_PNRZ_hs = conv(Signal_filtered_Polar_NRZ_hs, fliplr(hs));

eyediagram(match_filtered_PNRZ_rec,2*mp);
title('Eyediagrma rectangular');

eyediagram(match_filtered_PNRZ_hs,2*mp);
title('Eyediagrma half sine');

%%
%Exercise 6
%Resample the signal and apply a treshold and start sampling

delay_signal = ford/2 + numel(hs)/2; %calculate the signal delay, channel delay + match filter delay

start_recovery_count = delay_signal + mp/2;

%plotting the signal we can see that the maximum value resides between 12
%units por V+ and V- so a good treshold would be either plus 0 or 6

PNRZ_recovery_rec = match_filtered_PNRZ_rec(start_recovery_count:mp:end);
PNRZ_recovery_hs = match_filtered_PNRZ_hs(start_recovery_count:mp:end);
Decition_treshold_PNRZ = 0;

plot(PNRZ_recovery_rec);
title('Rectangular pulse sampled');

figure();
plot(PNRZ_recovery_hs);
title('Rectangular pulse sampled');

%%
%Exercise 7
%Graph the constelation

scatterplot(match_filtered_PNRZ_rec);
title('Rectangular pulse constelation');

scatterplot(match_filtered_PNRZ_hs);
title('Half-sine pulse constelation');


%%
%Exercise 8
%Recover the bits, calculate BER and recover the image

%start symbols to bit convertions

PNRZ_recovery_bits_rec = zeros(1,numel(V_bit));
PNRZ_recovery_bits_rec((PNRZ_recovery_rec > Decition_treshold_PNRZ)) = 1; 

PNRZ_recovery_bits_hs = zeros(1,numel(V_bit));
PNRZ_recovery_bits_hs((PNRZ_recovery_hs > Decition_treshold_PNRZ)) = 1; 

plot(PNRZ_recovery_bits_rec);
title('Recovered bits in rectangular pulses');

figure();
plot(PNRZ_recovery_bits_hs);
title('Recovered bits in half sine pulses');
%%
%calculate BER

%Half sine
bits_error_PNRZ_hs = xor(V_bit, PNRZ_recovery_bits_hs(1:numel(V_bit)));
error_PNRZ_hs = sum(bits_error_PNRZ_hs);
Bit_error_rate_PNRZ_hs = (error_PNRZ_hs/numel(V_bit)) * 100

%Rectangular
bits_error_PNRZ_rec = xor(V_bit, PNRZ_recovery_bits_rec(1:numel(V_bit)));
error_PNRZ_rec = sum(bits_error_PNRZ_rec);
Bit_error_rate_PNRZ_rec = (error_PNRZ_rec/numel(V_bit)) * 100

%%
%Recover the image

%Half sine
recuperado = zeros(pixels,pixels,'uint32'); %allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : pixels
    for j = 1: pixels
        recuperado(j,i) = bi2de(PNRZ_recovery_bits_hs(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end

nexttile;
imshow(uint8(recuperado)) %recovered image
title('Image recovered half-sine pulse');

%%
%Recover the image rectangular

%Bipolar
recuperado = zeros(42,42,'uint32'); %allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : pixels
    for j = 1: pixels
        recuperado(j,i) = bi2de(PNRZ_recovery_bits_rec(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end

nexttile;
imshow(uint8(recuperado)) %recovered image
title('Image recovered rectangular pulse');




















