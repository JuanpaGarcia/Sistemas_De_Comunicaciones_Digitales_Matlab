clear all;
close all;
clc 

%Tarea 8

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

plots = 1;

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

% Polar NRZ LineCode halfsine
Polar_NRZ_sig_hs = conv(hs ,V_bit_polar);

%change power to 1
Polar_NRZ_sig_hs = sqrt(1/((sum(Polar_NRZ_sig_hs.^2))/numel(Polar_NRZ_sig_hs))).*Polar_NRZ_sig_hs;

if plots == 1
    plot(Polar_NRZ_sig_hs);
    title('Señal PNRZ half sine');

    eyediagram(Polar_NRZ_sig_hs,2*mp);
    title('Eyediagrma half sine');
end

%%
%Low pass filter generation, communications channel
%F 0.3042
f=[0 0.6 0.6 1];
m=[1 1 0 0];
ford=60;
filter_1 = fir2(ford,f,m);

Signal_filtered_Polar_NRZ_hs = conv(Polar_NRZ_sig_hs, filter_1);

if plots == 1
    plot(Signal_filtered_Polar_NRZ_hs)
    title('Filtered signal Polar NRZ half sine');

    eyediagram(Signal_filtered_Polar_NRZ_hs,2*mp);
    title('Eyediagrma half sine');
end


%change power to 1
Signal_filtered_Polar_NRZ_hs = sqrt(1/((sum(Signal_filtered_Polar_NRZ_hs.^2))/numel(Signal_filtered_Polar_NRZ_hs))).*Signal_filtered_Polar_NRZ_hs;
if plots == 1
    plot(Signal_filtered_Polar_NRZ_hs)
    title('Normalized Filtered signal Polar NRZ half sine');
end

pow_PNRZ_hs= (sum(Signal_filtered_Polar_NRZ_hs.^2))/(numel(Signal_filtered_Polar_NRZ_hs));

%%
%NOISE 

PNoise= 0.015625*mp;

%noise signal
PNoise_PNRZ_hs = sqrt(PNoise)*randn(1,numel(Signal_filtered_Polar_NRZ_hs));

%Power
PNoise_PNRZ_pow_hs=var(PNoise_PNRZ_hs)

%SNR
SNR_dB_PNRZ_HS=10*log10(pow_PNRZ_hs/PNoise)

%ADD NOISE TO SIGNAL
noisy_PNRZ_hs= Signal_filtered_Polar_NRZ_hs + PNoise_PNRZ_hs;

%%

delay_mf = mp/2;

match_filter = fliplr(hs);
recovered_signal_PNRZ_hs = conv(match_filter, noisy_PNRZ_hs);

%%
%Estimador espectral de potencia
if plots ==1
    figure;
    pwelch(recovered_signal_PNRZ_hs,500,300,500,Fs,'power');
end

%%
%Diagrama de ojo
if plots ==1
    eyediagram(recovered_signal_PNRZ_hs,2*mp);
end

%%
%Muestreo
%PNRZ_HS

start= ford/2+(mp/2)+ delay_mf;

sampled_noisy_PNRZ_hs = recovered_signal_PNRZ_hs(start:mp:end);

if plots ==1
    scatterplot(sampled_noisy_PNRZ_hs);
end

treshold_PolarNRZ_HS=0;

PNRZ_recovery_bits_hs = zeros(1,numel(V_bit));
PNRZ_recovery_bits_hs((sampled_noisy_PNRZ_hs > treshold_PolarNRZ_HS)) = 1; 

if plots == 1
    figure();
    plot(PNRZ_recovery_bits_hs);
    title('Recovered bits in half sine pulses');
end


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

if plots == 1
    figure();
    imshow(uint8(recuperado)) %recovered image
    title('Image recovered half-sine pulse');
end


%%
%calculate BER

%Half sine
bits_error_PNRZ_hs = xor(V_bit, PNRZ_recovery_bits_hs(1:numel(V_bit)));
error_PNRZ_hs = sum(bits_error_PNRZ_hs);
Bit_error_rate_PNRZ_hs = (error_PNRZ_hs/numel(V_bit)) * 100

%%































