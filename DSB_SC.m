%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------setup-------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%{
we need to change the sampling speed because we now have the carrier plus
the message freq(3.4khz after filtering).
maxmimum freq we will have is 3.4+48=51.4
Nyquist sampling theory; we will choose fc to be bigger then 2*(51.4)
%}
sampling_freq=150000;
carrier_freq=48000;
[audio, samp_freq_message] = audioread('main_audio.wav');
audioData=resample(audio,sampling_freq,samp_freq_message);

%-------------------- Low-pass Filtering --------------------
wn=3400/(sampling_freq / 2);%filtered ratio
[num_coef,den_coef]=butter(4, wn, 'low');% Design the Butterworth low-pass filter 4th order
filtered_audio = filter(num_coef,den_coef, audioData);

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------General variabels-------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

N_samples = length(audioData);
time_vector=(0:(N_samples)-1);% for zero indexing
half_N = floor(N_samples/2);

t_seconds = (0:N_samples-1) / sampling_freq;  %to get time vector in seconds
freq_vector = linspace(-sampling_freq/2, sampling_freq/2, N_samples);

%--------------------------------------------------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------Modulation-------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%modulated_message=modulate(filtered_audio,carrier_freq,sampling_freq ,'amdsb-sc');
%demodulated_message=demod(modulated_message,carrier_freq,sampling_freq ,'amdsb-sc');
%{
x = y.*cos(2*pi*fc*t);
[b,a] = butter(5,fc*2/fs);
x = filtfilt(b,a,x);
filtfilt applies zero-phase filtering by running the filter forward and then backward,
eliminating phase distortion and doubling the filter order.
%}

carrier = cos(2 * pi * carrier_freq * t_seconds);
modulated_message = filtered_audio .* carrier';

received_carrier = cos(2 * pi * carrier_freq * t_seconds);
demodulated_fullterms = modulated_message .* received_carrier';%full terms means the terms rigth before LBF
demodulated_message = filter(num_coef, den_coef, demodulated_fullterms);
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%-----------------------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

fig_1=figure;
subplot(3,1,1);
plot(t_seconds, modulated_message);
title('DSB-SC Modulated Signal','FontSize',18);
xlabel('Time (s)');
ylabel('Amplitude (volt)');
grid on;

ft_modulated = fft(modulated_message);
ft_shifted_filtered_modulated = fftshift(ft_modulated);

fig_2=figure;
plot(freq_vector, abs(ft_shifted_filtered_modulated));
title('Spectrum of DSB-SC Modulated Signal','FontSize',18);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;    

%-------------------------------------------
sound(demodulated_message,sampling_freq);
pause((N_samples / sampling_freq)+1);
%--------------------------------------------

figure(fig_1);
subplot(3,1,2);
plot(t_seconds, demodulated_message);
title('DSB-SC Demodulated Signal','FontSize',18);
xlabel('Time (s)');
ylabel('Amplitude (volt)');
grid on;

ft_demod = fft(demodulated_message);
ft_shifted_filtered_demod = fftshift(ft_demod);

fig_3=figure;
plot(freq_vector, abs(ft_shifted_filtered_demod));
title('Spectrum of DSB-SC Demodulated Signal','FontSize',18);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------sound-----------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
Rx_LO_original_freq=carrier_freq;
offsets = [0,-8000,-5000, 5000, 8000];
for i = 1:length(offsets)
    freq_offset = offsets(i);
    fprintf('Demodulating with offset = %d Hz\n', freq_offset);
    % Demodulate with offset frequency
    demod_message_offseted = demod(modulated_message, Rx_LO_original_freq + freq_offset, sampling_freq, 'amdsb-sc');
    % Normalize audio to prevent clipping
    demod_message_offseted = demod_message_offseted / max(abs(demod_message_offseted));
    % LPF to remove high-frequency noise (sharp sounds)
    demod_message_offseted_filtered = filtfilt(num_coef,den_coef, demod_message_offseted);
    sound(demod_message_offseted_filtered, sampling_freq);
    pause((N_samples / sampling_freq)+1);%to get the full audio time in seconds+1second
end
