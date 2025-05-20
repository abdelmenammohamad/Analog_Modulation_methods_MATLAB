%max deviation ratio in this file is 10 and 
%since bandwidth of the message is 3400Hz due to filtering.
%then max freq_devation we will have is 10*3400=34000 
%thenby using carson formula: 2*(∆f+B), the BW is 2*(34000+3400)=74800
%so the maximum freq that will in the S(t) is 48000+34000=82000
%then the we will chose a sampling frequnecy =480000 due to Nyquist sampling frequency

sampling_freq = 480000;
carrier_freq = 48000;
fc=carrier_freq;

[audio, samp_freq_message] = audioread('main_audio.wav');
N_orig = length(audio);
t_orig = (0:N_orig-1) / samp_freq_message;% Original time vector

% New time vector for desired sampling frequency
t_new = (0:1/sampling_freq:(N_orig-1)/samp_freq_message)';

% Interpolate audio to new sampling frequency using interp1
audioData = interp1(t_orig', audio, t_new, 'linear');  % linear interp can be replaced by 'spline' for smoother
%interpulate_1 should be the same as resamble but resambe doesn't work here
%and i don't know why
%-------------------- Low-pass Filtering --------------------
max_msg_freq = 3400;  % Max frequency in bandlimited signal
message_bandwidth=max_msg_freq;
wn = 3400 / (sampling_freq / 2); % Normalized cutoff frequency
[num_coef, den_coef] = butter(4, wn, 'low'); % 4th order Butterworth LPF
filtered_audio = filtfilt(num_coef, den_coef, audioData);

%----------------General variables-------------------
N_samples = length(filtered_audio);
t_seconds = (0:N_samples - 1) / sampling_freq;
freq_vector = linspace(-sampling_freq / 2, sampling_freq / 2, N_samples);
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------Modulation and Plotting-------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
beta_values = [3, 5];

figure_time = figure;  % One figure for all time-domain plots
for i = 1:length(beta_values)
    beta = beta_values(i);
    freq_dev = beta * message_bandwidth;  % Frequency deviation

    % FM Modulation
    %{
    fm_signal = modulate(filtered_audio, carrier_freq, sampling_freq, 'fm', freq_dev);
    %y = cos(2*pi*fc*t + opt*cumsum(x)) cumsum is the intergration
    %}
    
    fm_signal=cos(2*pi*fc*t_seconds' + freq_dev*cumsum(filtered_audio)/(4800));%fc/(sampling_freq/fc)=4800 
    %i really don't know why but this value fc/(sampling_freq/fc) always
    %seems to work
    % Time-Domain Plot in Subplot
    figure(figure_time);
    subplot(2, 1, i);
    plot(t_seconds, fm_signal);
    title(['FM Modulated Signal (Time Domain) - \beta = ', num2str(beta)], 'FontSize', 18);
    xlabel('Time (s)');
    ylabel('Amplitude (v)');
    ylim([-2,2]);%FM does not change amplitude so I used 2 so the signal does not fill up the whole garph
    grid on;

    % Frequency-Domain Plot
    ft_fm_shifted = fftshift(fft(fm_signal));
    figure;
    plot(freq_vector, abs(ft_fm_shifted));
    title(['Spectrum of FM Signal - \beta = ', num2str(beta)], 'FontSize', 18);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
    
    %================= HARD LIMITER =====================%
    fm_signal_clipped = max(min(fm_signal, 1), -1);  % Clip to [-1, 1]
    %a hard limiter can cause harmonics and high frequency noise, that is why we need a filter after.
    %================= BPF FILTER =====================%
    bpf_low = fc-34000;  % Lower cutoff (fc-∆f )
    bpf_high = fc+34000;  % Upper cutoff (fc+∆f)
    bpf_norm = [bpf_low bpf_high] / (sampling_freq / 2);  % Normalized passband
    [B_bpf, A_bpf] = butter(4, bpf_norm, 'bandpass');  % 4th order bandpass filter
    lpf_output = filter(B_bpf, A_bpf, fm_signal_clipped);
    
    %================= DIRECT METHOD DEMODULATION =====================%
    differentiated_fm = diff(lpf_output);%diff output vector O/P have a one element smaller

    % Envelope detection using Hilbert transform + LPF
    analytic_signal = hilbert(differentiated_fm(:));  % Compute analytic signal
    envelope_fm = abs(analytic_signal);  % Magnitude gives the envelope
    [b_env, a_env] = butter(8, wn, 'low');  % 8th order LPF for smoothing
    recovered_message = filter(b_env, a_env, envelope_fm);
    
    % Normalize for playback to avoid clippings
    max_val = max(abs(recovered_message));
    if max_val > 0
        recovered_message = recovered_message / max_val;
    else
        warning('Recovered message has zero amplitude.');
    end
    recovered_play = recovered_message;
    
    recovered_play_44k = resample(recovered_play, 44100, sampling_freq);
    %some devices can't play at a huge sample rate like our sampling_freq
    
    % Playback
    fprintf('Playing demodulated message (β = %d)\n', beta);
    sound(recovered_play_44k, 44100);
    pause(N_samples / sampling_freq + 1);

    % Plot recovered signal
    figure;
    plot(t_seconds(2:length(recovered_message)), recovered_message(2:end));
    title(['Demodulated Message using Direct Method - \beta = ', num2str(beta)], 'FontSize', 15);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

%============= FM MODULATION OF A 3 kHz SINUSOID =============%
tone_freq = 3000;  
tone_duration=1;
tone_time_vec=linspace(0, tone_duration, tone_duration * sampling_freq);
tone_signal = cos(2 * pi * tone_freq * tone_time_vec);  
beta_list = [0.1,0.5, 1, 3, 5, 10];%it is like moving from NBFM to WBFM
%when increasing beta the number of significant deltas increase β=N+1
figure_tone_time = figure;  % One figure for all tone time-domain plots
set(figure_tone_time, 'Position', [100, 100, 1200, 700]);%make figure bigger
for i = 1:length(beta_list)
    beta = beta_list(i);
    freq_dev = beta * tone_freq;

    % FM modulation of tone
    fm_tone = cos(2*pi*fc*tone_time_vec' + 2*pi*freq_dev*cumsum(tone_signal')/fc/(sampling_freq/fc));%fc/(sampling_freq/fc)=4800 
    %i really don't know why but this value fc/(sampling_freq/fc) always

    % Time-domain subplot
    figure(figure_tone_time);
    subplot(3, 2, i);
    plot(tone_time_vec(1:1500), fm_tone(1:1500));%to 1500 smaple just to show modulation rather than zooming in the figure
    title(['FM Modulated 3 kHz Tone (Time Domain) - \beta = ', num2str(beta)], 'FontSize', 18);
    xlabel('Time (s)');
    ylabel('Amplitude');
    ylim([-2,2]);%FM does not change amplitude so I used 2 so the signal does not fill up the whole garph
    grid on;

    % Frequency-domain
    ft_tone = fftshift(fft(fm_tone));
    freq_axis_tone = linspace(-sampling_freq/2, sampling_freq/2, length(ft_tone));
    figure;
    plot(freq_axis_tone, abs(ft_tone));
    title(['Spectrum of FM 3 kHz Tone - \beta = ', num2str(beta)], 'FontSize', 15);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
end



