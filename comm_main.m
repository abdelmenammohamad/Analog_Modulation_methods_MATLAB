%Recorder;  % Run Recorder.m to get audiodata in Workspace
%N_samples=record_object.TotalSamples;

%------------------

[audioData, sampling_freq] = audioread('main_audio.wav');
N_samples = length(audioData);

%----------------------------------------------------------
time_vector=(0:(N_samples)-1);% for zero indexing
half_N = floor(N_samples/2);

%----------------Calculatings------------
% Frequency axis

ft = fft(audioData);%discrete fourier transfrom
ft_shifted= fftshift(ft);%%discrete fourier transfrom(shifted) zero freq in the middle

freq_vector = linspace(-sampling_freq,sampling_freq,N_samples);%length(freq_vector)=length(N_samples)

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%-------------Ploting------------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

fig_1=figure;
subplot(2,1,1);%message time domain
plot(time_vector, audioData,'m');
title('Time Domain Signal-continous time','FontSize',18);
xlabel('Time (s)');
ylabel('Amplitude (v)');
legend('message');
grid on;

fig_2=figure;%frequency domain positive and negative
semilogy(freq_vector,abs(ft_shifted),'m');  
title('Frequency Spectrum','FontSize',18);
xlabel('Frequency (Hz)');
ylabel('Magnitude (log scale)');
legend('message fourier transform');
grid on;

figure(fig_1);
subplot(2,1,2);%filtered message
wn=3400/(sampling_freq / 2);%filtered ratio
[num_coef,den_coef]=butter(4, wn, 'low');% Design the Butterworth low-pass filter 4th order
filtered_audio = filter(num_coef,den_coef, audioData);
plot(time_vector,filtered_audio,'m');
title('Filtered message','FontSize',18);
xlabel('Time (s)');
ylabel('voltages');
grid on;

fig_3=figure;
ft_filtered = fft(filtered_audio);%foureri transform after filtering
ft_shifted_filtered= fftshift(ft_filtered);
semilogy(freq_vector,abs(ft_shifted_filtered),'m');  
title('Filtered Frequency Spectrum','FontSize',18);
xlabel('Frequency (Hz)');
ylabel('Magnitude (log scale)');
legend('message-filtered fourier transform');
grid on;

cutoff_freq = 350;               % Start at 7 kHz
min_cutoff = 10;                 % End at 500 Hz
step = -400;                      % Step down by 500 Hz
filter_order = 4;

while cutoff_freq >= min_cutoff
    % Normalize cutoff frequency
    Wn = cutoff_freq / (sampling_freq / 2);

    % Design low-pass Butterworth filter
    [b, a] = butter(filter_order, Wn, 'low');

    % Apply the filter
    filtered_audio = filter(b, a, audioData);

    % Play the filtered audio
    fprintf('Playing with cutoff = %d Hz\n', cutoff_freq);
    sound(filtered_audio, sampling_freq);
    
    pause(length(filtered_audio) / sampling_freq + 0.5); % Wait for audio to finish + pause

    % Decrease cutoff
    cutoff_freq = cutoff_freq + step;
end

