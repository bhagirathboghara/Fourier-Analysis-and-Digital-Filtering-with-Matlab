% Define the sampling parameters
fs = 1000; % Sample rate in Hz
t = 0:1/fs:1-1/fs; % Time vector for 1 second

% Generate a composite signal with multiple frequency components
signal = sin(2*pi*50*t) + 0.5*sin(2*pi*120*t) + 0.3*sin(2*pi*250*t);

% Design a low-pass filter with a cutoff frequency of 100 Hz
fc = 100; % Cut-off frequency
[n, Wn] = buttord(fc/(fs/2), (fc+10)/(fs/2), 3, 40); % Calculate the minimum filter order
[b, a] = butter(n, Wn, 'low'); % Generate filter coefficients

% Apply a window to the original signal to reduce spectral leakage
windowed_signal = signal .* hann(length(signal))';  % Using Hann window

% Apply the filter to the windowed signal
filtered_signal = filter(b, a, windowed_signal);

% Compute the frequency spectrum of the filtered signal
L = length(filtered_signal); % Length of the signal
Y = fft(filtered_signal); % Compute the Fourier Transform
P2 = abs(Y/L); % Two-sided spectrum P2
P1 = P2(1:L/2+1); % Single-sided spectrum P1
P1(2:end-1) = 2*P1(2:end-1); % Since we dropped half the points we multiply by 2 for the actual amplitude
f = fs*(0:(L/2))/L; % Define the frequency domain f

% Plot the original and filtered signal
figure;
subplot(2,1,1); % Subplot for the original signal
plot(t, signal);
title('Original Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');

% Compute dB magnitude
P1_dB = 20*log10(P1);

% Plot the frequency spectrum in dB
figure;
plot(f, P1_dB);
title('Frequency Spectrum of the Filtered Signal (dB)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');


% Plot the frequency spectrum of the filtered signal
%subplot(2,1,2); % Subplot for the frequency spectrum
%plot(f, P1);
%title('Frequency Spectrum of the Filtered Signal');
%xlabel('Frequency (Hz)');
%ylabel('Magnitude');
