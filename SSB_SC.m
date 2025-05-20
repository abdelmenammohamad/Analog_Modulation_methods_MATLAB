sampling_freq = 150000;
carrier_freq = 48000;

[audio, samp_freq_message] = audioread('main_audio.wav');
audioData = resample(audio, sampling_freq, samp_freq_message);

%-------------------- Low-pass Filtering --------------------

wn = 3400 / (sampling_freq / 2); % Normalized cutoff frequency
[num_coef, den_coef] = butter(4, wn, 'low'); % 4th order Butterworth LPF
filtered_audio = filter(num_coef, den_coef, audioData);

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------General variabels--------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

N_samples = length(filtered_audio);
t_seconds = (0:N_samples - 1) / sampling_freq;
freq_vector = linspace(-sampling_freq / 2, sampling_freq / 2, N_samples);

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------Modulation---------------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

fc = carrier_freq;
% Lower Sideband SSB modulation 
modulated_ssb_sc = filtered_audio .* cos(2 * pi * fc * t_seconds)' + ...
                   imag(hilbert(filtered_audio)) .* sin(2 * pi * fc * t_seconds)';

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------Demodulation-------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%suppose follwing equation is Coherent detection
demodulated_fullterms = modulated_ssb_sc .* cos(2 * pi * carrier_freq * t_seconds') * 2;

% Low-pass filter
[b_demod, a_demod] = butter(4, wn, 'low');
demodulated_message = filter(b_demod, a_demod, demodulated_fullterms);

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%-----------------Ploting------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

fig_1 = figure;
subplot(2,1,1);
plot(t_seconds, modulated_ssb_sc);
title('SSB-SC Modulated Signal (Time Domain)', 'FontSize', 18);
xlabel('Time (s)');
ylabel('Amplitude (volt)');
grid on;

ft_ssb = fft(modulated_ssb_sc);
ft_ssb_shifted = fftshift(ft_ssb);

fig_2 = figure;
plot(freq_vector, abs(ft_ssb_shifted));
title('SSB-SC Modulated Signal Spectrum (Frequency Domain)', 'FontSize', 18);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

figure(fig_1);
subplot(2,1,2);
plot(t_seconds, demodulated_message);
title('SSB-SC Demodulated Signal (Time Domain)', 'FontSize', 18);
xlabel('Time (s)');
ylabel('Amplitude (volt)');
grid on;

ft_demod = fft(demodulated_message);
ft_demod_shifted = fftshift(ft_demod);

fig_3 = figure;
plot(freq_vector, abs(ft_demod_shifted));
title('SSB-SC Demodulated Signal Spectrum (Frequency Domain)', 'FontSize', 18);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------sound-----------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

sound(demodulated_message, sampling_freq);
pause((N_samples / sampling_freq) + 1); % Full audio length + 1 second

offsets = [0, -8000, -5000, 5000, 8000]; 

% Loop over different frequency offsets
for i = 1:length(offsets)
    freq_offset = offsets(i);
    fprintf('Demodulating with offset = %d Hz\n', freq_offset);
    
    % Demodulate with frequency offset
    demodulated_message_offseted = demod(modulated_ssb_sc, carrier_freq + freq_offset, sampling_freq, 'amssb');
    
    % Normalize to prevent clipping
    demodulated_message_offseted = demodulated_message_offseted / max(abs(demodulated_message_offseted));
    
    % Play the demodulated message
    sound(demodulated_message_offseted, sampling_freq);
    pause((N_samples / sampling_freq) + 1); % Full audio time + 1 second
end
