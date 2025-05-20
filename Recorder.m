sampling_freq = 48000;             % Sampling frequency 
%Nyquist Theorem, explained in your NOTES.txt
duration = 1.5;           % Durationin in seconds
record_object = audiorecorder(sampling_freq, 16, 1);  % 16-bit, mono channel
%audiorecorder() return audiorecorder object which has multiple other
%properites about the audio
disp('start speaking');
recordblocking(record_object, duration);%start recording
disp('recording ended');


audioData = getaudiodata(record_object);%just to get the vector that has all the arrays.
audiowrite('D.wav', audioData, sampling_freq);