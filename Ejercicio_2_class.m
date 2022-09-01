% Test Bench Configuration
frameLength = 1024;  
% Audio object for reading an audio file
fileReader = dsp.AudioFileReader( ...
'Counting-16-44p1-mono-15secs.wav', ...
'SamplesPerFrame',frameLength);
% Audio object for writing the audio device
deviceWriter = audioDeviceWriter( ...
'SampleRate',fileReader.SampleRate);
% Object for Visualization
scope = timescope( ...
'SampleRate',fileReader.SampleRate, ...
'TimeSpanSource', 'property', ...
'TimeSpan', 2, ...
'BufferLength',fileReader.SampleRate*2*2, ...
'YLimits',[-1,1], ...
'TimeSpanOverrunAction',"Scroll");
% Quantization process (Q)
k=16;                              % Cantidad de bits
swing = (2^k-1)/2;         %   Asumir señal en el rango [-1 1] y usar                    
% un rango simétrico
% Audio Stream Loop
while ~isDone(fileReader)                   %
signal = fileReader();                        %
xq_int = round(signal*swing+swing); % Convertir a enteros en el 
% rango [0 2^k-1] usando k bits 
xq = (xq_int-swing)/swing;     % Proceso inverso
scope([signal,xq])                  % Visualizacion
deviceWriter(xq);                    % Enviar la señal al driver de audio
end
release(fileReader)                         
release(deviceWriter)                      
release(scope) 