% Insert Description

micro = 10^(-6);
mega = 10^(6);

ti = 0; % Begin of the computational window
tf = 200*micro; % End of the computational window

fs = 50*mega; % Sampling frequency
Ts = 1/fs; % Sampling period

N = fs*((tf-ti)); % Number of samples

n = (0:1:N-1);
t = n*Ts;
A = 1.1;

fc = 0.5*mega;

c = A*cos(t.*(fc*2*pi));

[time, carrier] = delimit_window(t, c, 0, 5*micro);

% Converter para subplot
figure(1)
plot(time./micro, carrier)
xlabel('Time (us)');
ylabel('c(t)');
title('Plot of the Carrier c(t)')

c_ft = fftshift(fft(c))*Ts;

v = (-1/2:1/N:1/2-1/N);

f = v.*fs;

[neg_freq_values, neg_freq_spectrum] = delimit_window(f, c_ft, -2.01*mega, -1.8*mega);

[pos_freq_values, pos_freq_spectrum] = delimit_window(f, c_ft, 1.8*mega, 2.01*mega);

% Converter para subplot
figure(2)
% Negative part of the spectrum
subplot(2, 2, 1);
plot(neg_freq_values./mega, abs(neg_freq_spectrum));
xlabel('Frequency (MHz)');
ylabel('Magnitude');
title('Negative Part of the Spectrum');

% Positive part of the spectrum
subplot(2, 2, 2);
plot(pos_freq_values./mega, abs(pos_freq_spectrum));
xlabel('Frequency (MHz)');
ylabel('Magnitude');
title('Positive Part of the Spectrum')


m = sinc((t - 100*micro)*mega);
m_ft = fftshift(fft(m))*Ts;

[time, message] = delimit_window(t, m, 90*micro, 110*micro);

% Converter para subplot
figure(3)
plot(time./micro, message)
xlabel('Time (us)');
ylabel('m(t)');
title('Plot of the message m(t)')

[frequency, message_freq_spectrum] = delimit_window(f, m_ft, -2.01*mega, 2.01*mega);

% Converter para subplot
figure(4)
plot(frequency./mega, abs(message_freq_spectrum));
xlabel('Frequency (MHz)');
ylabel('Magnitude');
title('Message Spectrum');

[ bw ] = half_power_bandwidth(f, m_ft);

disp(['Half-Power Bandwidth: ', num2str(bw), ' Hz']);

mc = c .* m;

[time, modulated_message] = delimit_window(t, mc, 90*micro, 110*micro);

% Converter para subplot
figure(5)
plot(time./micro, modulated_message)
xlabel('Time (us)');
ylabel('mc(t)');
title('Plot of the product message m(t) * c(t)')

mc_ft = fftshift(fft(mc))*Ts;

[frequency, mod_m_freq_spectrum] = delimit_window(f, mc_ft, -5.01*mega, 5.01*mega);

% Converter para subplot
figure(6)
plot(frequency./mega, abs(mod_m_freq_spectrum));
xlabel('Frequency (MHz)');
ylabel('Magnitude');
title('Modulated Message Spectrum');

dm = c .* mc * (2/A^2); % correction factor

dmc_ft = fftshift(fft(dm))*Ts;

[frequency, dmod_m_freq_spectrum] = delimit_window(f, dmc_ft, -6.01*mega, 6.01*mega);

% Converter para subplot
figure(7)
plot(frequency./mega, abs(dmod_m_freq_spectrum));
xlabel('Frequency (MHz)');
ylabel('Magnitude');
title('Demodulated Message Spectrum');

% Initialize the filter
rect_filter = f*0;

rect_filter(find(f >= -2.0*mega & f <= 2.0*mega)) = 1;
%rect_filter = ones(1, N);
[frequency, filter_freq_spectrum] = delimit_window(f, rect_filter, -6.01*mega, 6.01*mega);

figure(8);

% Plot the signal on the left y-axis
plot(frequency./mega, abs(dmod_m_freq_spectrum), 'b', 'LineWidth', 1.5); % Plot signal in red
hold on;

% Create a second y-axis using the first plot as reference
ax1 = gca; % Current axes
% Create a second y-axis using the first plot as reference
ax1_pos = get(ax1, 'Position'); % Position of first axes
ax2 = axes('Position', ax1_pos, 'XAxisLocation', 'bottom', 'YAxisLocation', 'right', 'Color', 'none');

% Plot the filter on the right y-axis
line(frequency./mega, abs(filter_freq_spectrum), 'Parent', ax2, 'Color', 'r', 'LineWidth', 1.5); % Plot filter in blue

% Add labels and title
xlabel(ax1, 'Frequency (Hz)');
ylabel(ax1, 'Signal Magnitude');
ylabel(ax2, 'Filter Magnitude');
title(ax1, 'Signal and Filter in Frequency Domain');

% Set the colors of the axes to match the plot colors
%set(ax1, 'YColor', 'r');
%set(ax2, 'YColor', 'b');

% Add grid to the first axis (optional)
grid on;

% Add legend
legend(ax1, 'Signal');
legend(ax2, 'Filter');
hold off;

recuperated_message_ft = dmc_ft .* rect_filter;


[frequency, rec_message_freq_spectrum] = delimit_window(f, recuperated_message_ft, -6.01*mega, 6.01*mega);

% Converter para subplot
figure(9)
plot(frequency./mega, abs(rec_message_freq_spectrum));
xlabel('Frequency (MHz)');
ylabel('Magnitude');
title('Post LPF Message Spectrum');

rec_message = ifft(fftshift(recuperated_message_ft./Ts));

[time, recuperated_message] = delimit_window(t, rec_message, 90*micro, 110*micro);
[time, original_message] = delimit_window(t, m, 90*micro, 110*micro);

% Converter para subplot
figure(10)
plot(time./micro, real(recuperated_message), 'b', 'LineWidth', 1.5);
hold on;
plot(time./micro, real(original_message), 'r');
hold off;
xlabel('Time (\mus)');
ylabel('Amplitude');
title('Original Message and Demodulated Message');
legend( 'Demodulated Message dm(t)', 'Original Message m(t)');
grid on;

corr_matrix = corrcoef(m, rec_message);
similarity = corr_matrix(1, 2);

disp(['Correlation Coefficient between m(t) and recovered m(t): ', num2str(real(similarity))]);
