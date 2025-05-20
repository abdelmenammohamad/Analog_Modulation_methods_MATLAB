%max deviation ratio in this file is 10 and 
%since bandwidth of the message is 3400Hz due to filtering.
%then max freq_devation we will have is 10*3400=34000 
%then by using carson formula: 2*(∆f+B), the BW is 2*(34000+3400)=74800
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
    

    % Frequency-Domain Plot
    ft_fm_shifted = fftshift(fft(fm_signal));
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
    [b_env, a_env] = butter(4, wn, 'low');  % 8th order LPF for smoothing
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
    %fprintf('Playing demodulated message (β = %d)\n', beta);
    %sound(recovered_play_44k, 44100);
    %pause(N_samples / sampling_freq + 1);
    
end

   


noise_power_level = [0.0000001, 0.01, 0.1];  % Example noise power levels
for i = 1:length(beta_values)
    beta = beta_values(i);
    freq_dev = beta * message_bandwidth;

    % FM modulation (again, to isolate from previous runs)
    fm_audio_clean = cos(2*pi*fc*t_seconds' + 2*pi*freq_dev*cumsum(filtered_audio)/(sampling_freq));

    for j = 1:length(noise_power_level)
        noise_power = noise_power_level(j);
        random=randn(size(fm_audio_clean));
        noise = sqrt(noise_power) *random;
        fm_audio_noisy = fm_audio_clean + noise;
        
        signal_power_input = mean(fm_audio_clean.^2);
        noise_power_input = mean(noise.^2);
        input_snr = 10 * log10(signal_power_input / (noise_power_input + eps));%eps is to ensure no zero values
        
        %====================================================%
        %================= Demodulation =====================%
        %====================================================%
        
        %================= HARD LIMITER =====================%
        limited_tone = max(min(fm_audio_noisy, 1), -1);

        %================= BPF =====================%
        [B_bpf, A_bpf] = butter(8, [fc-34000 fc+34000]/(sampling_freq/2), 'bandpass');
        filtered_tone = filter(B_bpf, A_bpf, limited_tone);

        %================= DIFFERENTIATION + HILBERT =====================%
        diff_tone = diff(filtered_tone);
        analytic = hilbert(diff_tone(:));
        envelope = abs(analytic);
        envelope=envelope - mean(envelope);
        [b_env, a_env] = butter(8, 3400/(sampling_freq/2), 'low');
        recovered = filter(b_env, a_env, envelope);

        % Normalize and resample
        recovered = recovered / max(abs(recovered));  
        recovered_resampled = resample(recovered, 44100, sampling_freq);
        
        % Output SNR calculation
        signal_power_output = mean(recovered.^2);
        snr_output = 10 * log10(signal_power_output / (noise_power + eps));

        % Print both input and output SNRs
        fprintf('β = %g | Noise Power = %.6f | SNR_in = %.2f dB | SNR_out = %.2f dB\n', ...
        beta, noise_power, input_snr, snr_output);
        
        % Playback (optional)
        fprintf('Playing β = %g with noise power = %.4f\n', beta, noise_power);
        sound(recovered_resampled, 44100);
        pause(N_samples/sampling_freq + 0.5);
        
        
    end
end


%================= Threshold Effect Test (Manual beta) =================%
fprintf('\n================ Threshold Effect Test ================\n');
noise_Level_fixed = 0.00001;  % fixed noise power as you requested
random=randn(size(fm_audio_clean));
noise_fixed = sqrt(noise_Level_fixed) *random ;
noise_power_fixed = mean(noise.^2);

beta =input('Enter the beta value: ');
freq_dev = beta * message_bandwidth;

% FM Modulation
fm_mod = cos(2*pi*fc*t_seconds' + 2*pi*freq_dev*cumsum(filtered_audio)/(sampling_freq));

% Add Noise
noise = sqrt(noise_power_fixed) * randn(size(fm_mod));
fm_noisy = fm_mod + noise;

% Hard limiter
limited = max(min(fm_noisy, 1), -1);

% BPF
[B_bpf, A_bpf] = butter(8, [fc-34000 fc+34000]/(sampling_freq/2), 'bandpass');
filtered = filter(B_bpf, A_bpf, limited);

% Differentiation + Hilbert
diff_sig = diff(filtered);
analytic = hilbert(diff_sig(:));
envelope = abs(analytic);
envelope = envelope - mean(envelope);
[b_env, a_env] = butter(8, 3400/(sampling_freq/2), 'low');
recovered = filter(b_env, a_env, envelope);

% Output SNR
signal_power_output = mean(recovered.^2);
snr_out = 10 * log10(signal_power_output / (noise_power_fixed + eps));

fprintf('β = %.1f | Output SNR = %.2f dB\n', beta, snr_out);

