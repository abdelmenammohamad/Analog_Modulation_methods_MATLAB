% Path to the folder containing the WAV files
folder_path = fullfile(pwd, 'Letters');  % Gets current folder + 'Letters'

% List of files to process
file_list = {'F.wav', 'B.wav', 'M.wav', 'N.wav', 'D.wav', 'S.wav'};

for k = 1:length(file_list)
    filename = file_list{k};
    full_filename = fullfile(folder_path, filename);
    fprintf('Processing file: %s\n', full_filename);

    % Load audio file
    [audioData, sampling_freq] = audioread(full_filename);
    N_samples = length(audioData);
    time_vector = (0:(N_samples-1)) / sampling_freq;  % In seconds
    freq_vector = linspace(-sampling_freq, sampling_freq, N_samples);

    % Frequency domain (original)
    ft = fft(audioData);
    ft_shifted = fftshift(ft);

    % Plot time domain (original)
    fig_1 = figure('Name', ['Time Domain - ', filename]);
    subplot(2,1,1);
    plot(time_vector, audioData, 'm');
    title(['Time Domain Signal - ', filename], 'FontSize', 18);
    xlabel('Time (s)');
    ylabel('Amplitude (v)');
    legend('Original Message');
    grid on;

    % Plot frequency domain (original) with log scale
    fig_2 = figure('Name', ['Frequency Spectrum - ', filename]);
    semilogy(freq_vector, abs(ft_shifted), 'm');
    title(['Frequency Spectrum - ', filename], 'FontSize', 18);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (log scale)');
    legend('Original FFT');
    grid on;

    % Filtering for 3.4 kHz display
    wn = 3400 / (sampling_freq / 2);
    [num_coef, den_coef] = butter(4, wn, 'low');
    filtered_audio = filter(num_coef, den_coef, audioData);

    % Plot time domain (filtered)
    figure(fig_1);
    subplot(2,1,2);
    plot(time_vector, filtered_audio, 'r');
    title(['Filtered Signal - ', filename], 'FontSize', 18);
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Filtered Message');
    grid on;

    % Plot frequency domain (filtered)
    ft_filtered = fft(filtered_audio);
    ft_shifted_filtered = fftshift(ft_filtered);
    fig_3 = figure('Name', ['Filtered Spectrum - ', filename]);
    semilogy(freq_vector, abs(ft_shifted_filtered), 'r');
    title(['Filtered Frequency Spectrum - ', filename], 'FontSize', 18);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (log scale)');
    legend('Filtered FFT');
    grid on;

    % -------- Cutoff Sweep Loop with Skip Option ----------
    cutoff_freq = 800;   % Start
    min_cutoff = 100;    % End
    step = -50;

    while cutoff_freq >= min_cutoff
        Wn = cutoff_freq / (sampling_freq / 2);
        [b, a] = butter(4, Wn, 'low');
        filtered_audio = filter(b, a, audioData);

        % Ask user what to do
        fprintf('\nReady to play %s at cutoff = %d Hz\n', filename, cutoff_freq);
        user_input = '';
        while ~strcmp(user_input, 'c') && ~strcmp(user_input, 's')
            user_input = lower(input('Type "c" to continue, "s" to skip this letter: ', 's'));
        end

        if strcmp(user_input, 's')
            fprintf('Skipping remaining cutoffs for %s...\n', filename);
            break;  % Exit the cutoff loop, go to next letter
        end

        % If 'c' → play the filtered sound
        sound(filtered_audio, sampling_freq);
        pause(length(filtered_audio)/sampling_freq + 0.5);

        cutoff_freq = cutoff_freq + step;
    end
end
