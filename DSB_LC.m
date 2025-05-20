%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------setup-------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%{
we need to change the sampling speed because we now have the carrier plus
the message freq(3.4khz after filtering).
maxmimum freq we will have is 3.4+48=51.4
Nyquist sampling theory; we will choose fc to be bigger then 2*(51.4)
%}
sampling_freq=128500;
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


mod_index = 0.8;
Ac = max(abs(filtered_audio)) / mod_index;

t_seconds = (0:N_samples-1) / sampling_freq;  %to get time vector in seconds
%--------------------------------------------------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------Modulation-------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
fc=carrier_freq;
carrier = cos(2 * pi * fc * t_seconds);
modulated_message = ( Ac+filtered_audio) .* carrier';
%{
opt=-Ac;
modulated_message = modulate(filtered_audio, carrier_freq, sampling_freq, 'amdsb-tc',opt);
y = (x-opt).*cos(2*pi*fc*t)
%}
ft = fft(modulated_message);
ft_shifted_filtered = fftshift(ft);
freq_vector = linspace(-sampling_freq/2, sampling_freq/2, N_samples);

%demodulation is later in the file

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%-----------------Ploting------------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


%-------------time domain of DSB-LC--------------
fig_1=figure;
subplot(3,1,1);
plot(t_seconds, modulated_message);
title('DSB-LC Modulated Signal','FontSize',18);
xlabel('Time (s)');
ylabel('Amplitude (volt)');
grid on;


%-------------frequency domain of the DSB-lC--------------



fig_2=figure;
plot(freq_vector, abs(ft_shifted_filtered));
title('Spectrum of DSB-LC Modulated Signal','FontSize',18);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


%-------------Envlope in time domian--------------
%///////envelope method one//////////////////
%[yupper,ylower] = envelope(modulated_message);
%message_demod= yupper -mean(yupper);
% OR you can do this message_demod= yupper - abs(min(yupper)/mod_index_estimation);
%problem is that envelope method one is almost ideal in values
%///////method two:square-law detection//////////////////
squared_signal = modulated_message .^ 2;
%Ac^2+2*Ac*m(t)+m(t)^2*cos(2*pi*fc*t)^2
% Low-pass filter extract envelope 
wn_sq = 3400 / (sampling_freq / 2);  
[b_LBF, a_LBF] = butter(4, wn_sq, 'low');
envelope_squared = filter(b_LBF, a_LBF, squared_signal);

% Remove DC component which will be AC^2
message_demod = envelope_squared - mean(envelope_squared);
max_amp = max(abs(message_demod));

figure(fig_1);
subplot(3,1,2);

plot(t_seconds,message_demod);
ylim([-max_amp, max_amp]);
title('demodulated message','FontSize',18);
xlabel('Time (s)');
ylabel('amplitude (v)');
grid on;





energy_filtered = sum(filtered_audio .^ 2);
energy_demod = sum(message_demod .^ 2);
scaling_factor = sqrt(energy_filtered / energy_demod);
message_demod_scaled= message_demod* scaling_factor;
figure;
subplot(2,1,1);
plot(t_seconds, filtered_audio, 'b');
xlabel('Time (s)');
ylabel('Amplitude (volt)');
grid on;
title('Original Demodulated Signal after filtering','FontSize',18);
subplot(2,1,2);
plot(t_seconds, message_demod_scaled);
xlabel('Time (s)');
ylabel('Amplitude (volt)');
title('Scaled Demodulated Signal','FontSize',18);
grid on;

figure;
ft_demod = fft(message_demod_scaled);
ft_demod_shifted = fftshift(ft_demod);
% Plot spectrum
plot(freq_vector, abs(ft_demod_shifted));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Demodulated Signal Frequency Domain','FontSize',18);
grid on;



%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%----------------sound-----------------
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

sound(message_demod_scaled, sampling_freq);




