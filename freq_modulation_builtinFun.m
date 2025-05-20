%max deviation ratio in this file is 10 and 
%since bandwidth of the message is 3400Hz due to filtering.
%then max freq_devation we will have is 10*3400 
%thenby using carson formula: 2*(∆f+B), the max BW is 2*(34000+3400)=74800
%so the maximum freq that will in the S(t) is 48000+34000=82000
%then the we will chose a sampling frequnecy =200000 due to Nyquist sampling frequency

sampling_freq = 200000;
carrier_freq = 48000;


[audio, samp_freq_message] = audioread('recorded_audio.wav');
audioData = resample(audio, sampling_freq, samp_freq_message);

%-------------------- Low-pass Filtering --------------------
max_msg_freq = 3400;  % Max frequency in bandlimited signal
wn = 3400 / (sampling_freq / 2); % Normalized cutoff frequency
[num_coef, den_coef] = butter(4, wn, 'low'); % 4th order Butterworth LPF
filtered_audio = filter(num_coef, den_coef, audioData);

%----------------General variables-------------------
N_samples = length(filtered_audio);
t_seconds = (0:N_samples - 1) / sampling_freq;
freq_vector = linspace(-sampling_freq / 2, sampling_freq / 2, N_samples);

%----------------Modulation and Plotting-------------------
beta_values = [3, 5];

figure_time = figure;  % One figure for all time-domain plots
for i = 1:length(beta_values)
    beta = beta_values(i);
    freq_dev = beta * max_msg_freq;  % Frequency deviation

    % FM Modulation using built-in function
    fm_signal = modulate(filtered_audio, carrier_freq, sampling_freq, 'fm', freq_dev);

    % Time-Domain Plot in Subplot
    figure(figure_time);
    subplot(2, 1, i);
    plot(t_seconds, fm_signal);
    title(['FM Modulated Signal (Time Domain) - \beta = ', num2str(beta)]);
    xlabel('Time (s)');
    ylabel('Amplitude (v)');
    ylim([-2,2]);%FM does not change amplitude so I used 2 so the signal does not fill up the whole garph
    grid on;

    % Frequency-Domain Plot
    ft_fm_shifted = fftshift(fft(fm_signal));
    figure;
    plot(freq_vector, abs(ft_fm_shifted));
    title(['Spectrum of FM Signal - \beta = ', num2str(beta)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;

    %================= DEMODULATION using demod() =====================%
    recovered_message = demod(fm_signal, carrier_freq, sampling_freq, 'fm', freq_dev);
    %y = cos(2*pi*fc*t + opt*cumsum(x)) 

    % Normalize for playback to avoid clippings
    recovered_play = recovered_message / max(abs(recovered_message));

    % Playback
    fprintf('Playing demodulated message (β = %d)\n', beta);
    sound(recovered_play, sampling_freq);
    pause(N_samples / sampling_freq + 1);

    % Plot recovered signal
    figure;
    plot(t_seconds, recovered_message);
    title(['Demodulated Message using demod() - \beta = ', num2str(beta)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

%============= FM MODULATION OF A 3 kHz SINUSOID =============
tone_freq = 3000;  % 3 kHz
tone_signal = sin(2 * pi * tone_freq * t_seconds);  % Message signal
beta_list = [0.5, 1, 3, 5, 10];

figure_tone_time = figure;  % One figure for all tone time-domain plots
set(figure_tone_time, 'Position', [100, 100, 1200, 700]);%make figure bigger
for i = 1:length(beta_list)
    beta = beta_list(i);
    freq_dev = beta * tone_freq;

    % FM modulation of tone using built-in function
    fm_tone = modulate(tone_signal', carrier_freq, sampling_freq, 'fm', freq_dev);

    % Time-domain subplot
    figure(figure_tone_time);
    subplot(3, 2, i);
    plot(t_seconds(1:1000), fm_tone(1:1000));
    title(['FM Modulated 3 kHz Tone (Time Domain) - \beta = ', num2str(beta)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    ylim([-2,2]);%FM does not change amplitude so I used 2 so the signal does not fill up the whole garph
    grid on;

    % Frequency-domain
    ft_tone = fftshift(fft(fm_tone));
    figure;
    plot(freq_vector, abs(ft_tone));
    title(['Spectrum of FM 3 kHz Tone - \beta = ', num2str(beta)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
end
