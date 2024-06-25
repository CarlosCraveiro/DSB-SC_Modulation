% Project: SEL0616 - Communication Principles Final Assignment
% Filename: main.m
% Description: It performs an DSB-SC modulation in a sinc function (the message) 
%
% Author(s): Carlos Craveiro - USP ID 12547187
%            Ivan Pancheniak - USP ID 12624224
%            Beatriz Aimee   - USP ID 12547934
%
% Created:  2024-24-06
% Modified: 2024-25-06
%
% This file is part of the SEL0616 - Communication Principles Final Assignment project.
%
% main.m is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Authors Observation: The titles and legends of the plots are in Brazilian Portuguese, 
%                      as they were used as part of the final report.
%

%% Defines general parameters and samples the Carrier Wave

% Defines the units modifiers
micro = 10^(-6);
mega = 10^(6);

% Verifies if fc is alrealdy defined to allow for user defined fc (carrier frequency)
if exist('fc', 'var')
    printf('Using fc already defined. [fc = %.1f MHz]\n', fc/mega);
else
    fc = 2*mega;
    printf('Using fc = %.1f MHz\n', fc/mega);
end

ti = 0; % Begin of the computational time window
tf = 200*micro; % End of the computational time window

fs = 50*mega; % Sampling frequency
Ts = 1/fs; % Sampling period

N = fs*((tf-ti)); % Number of samples

n = (0:1:N-1); % Defines n
t = n*Ts; % Defines time vector
A = 1.1; % Carrier amplitude

c = A*cos(t.*(fc*2*pi)); % Samples the carrier wave

purple = [0.58, 0, 0.82];
magenta = [0.82,0,0.50];

%% The first plot : Carrier Wave Time Plot

[time, carrier] = delimit_window(t, c, 0, 5*micro); % Delimits the Plot Window

figure(1)
title_text = sprintf('Portadora c(t), fc = %.1f MHz', fc/mega);
file_path = sprintf('../plots/01_fc_%d_kHz_carrier_time_plot.png', fc/1000);
plot(time./micro, carrier, 'Color', purple, 'LineWidth', 1.5);
xlabel('Tempo (\mus)');
ylabel('Amplitude');
title(title_text)
axis([0 5]);
grid on;
print(file_path, '-dpng');

%% Performs the Fourrier Transform of the Carrier Wave

c_ft = fftshift(fft(c))*Ts; % Obtain samples of the fourier transform

v = (-1/2:1/N:1/2-1/N); % Defines the admentional frequency vector

f = v.*fs; % Defines the frequency vector


%% The second plot : Carrier Wave Frequency Plot
% Delimits the Plot Windows
[neg_freq_values, neg_freq_spectrum] = delimit_window(f, c_ft, (-fc - 0.01*mega), (-fc + 0.2*mega));
[pos_freq_values, pos_freq_spectrum] = delimit_window(f, c_ft, (fc - 0.2*mega), (fc + 0.01*mega));

figure(2)
title_text = sprintf('Espectro negativo, fc = %.1f MHz', fc/mega);
title_text2 = sprintf('Espectro positivo, fc = %.1f MHz', fc/mega);
file_path = sprintf('../plots/02_fc_%d_kHz_carrier_freq_plot.png', fc/1000);

% Negative part of the spectrum
subplot(2, 2, 1);
plot(neg_freq_values./mega, abs(amp_normalize(neg_freq_spectrum)), 'Color', purple, 'LineWidth', 1.5);
xlabel('Frequência (MHz)');
ylabel('Magnitude Normalizada');
title(title_text);
axis([(-fc/mega - 0.01) (-fc/mega + 0.2)]);
grid on;

% Positive part of the spectrum
subplot(2, 2, 2);
plot(pos_freq_values./mega, abs(amp_normalize(pos_freq_spectrum)), 'Color', purple, 'LineWidth', 1.5);
xlabel('Frequência (MHz)');
ylabel('Magnitude Normalizada');
title(title_text2);
axis([(fc/mega - 0.2), (fc/mega + 0.01)]);
grid on;

print(file_path, '-dpng');

%% Defines the messagem sinc(x) centered at time 100 microsseconds (us)

m = sinc((t - 100*micro)*mega); % Samples the message

%% The third plot : Message Time Plot
[time, message] = delimit_window(t, m, 90*micro, 110*micro); % Delimits the Plot Window

figure(3)
title_text = sprintf('Mensagem m(t), fc = %.1f MHz', fc/mega);
file_path = sprintf('../plots/03_fc_%d_kHz_message_time_plot.png', fc/1000);
plot(time./micro, message, 'Color', purple, 'LineWidth', 1.5);
xlabel('Tempo (\mus)');
ylabel('Amplitude');
title(title_text);
axis([90 110]);
grid on;
print(file_path, '-dpng');

%% Performs the Fourrier Transform of the Message
m_ft = fftshift(fft(m))*Ts; % Obtain samples of the fourier transform 

%% The fourth plot : Message Frequency Plot
% Delimits the Plot Window
[frequency, message_freq_spectrum] = delimit_window(f, m_ft, -2.0*mega, 2.0*mega); 

figure(4)
title_text = sprintf('Espectro da mensagem m(t), fc = %.1f MHz', fc/mega);
file_path = sprintf('../plots/04_fc_%d_kHz_message_freq_plot.png', fc/1000);
plot(frequency./mega, abs(amp_normalize(message_freq_spectrum)), 'Color', purple, 'LineWidth', 1.5);
xlabel('Frequência (MHz)');
ylabel('Magnitude Normalizada');
title(title_text);
axis([-2.0 +2.0]);
grid on;
print(file_path, '-dpng');

%% Calculates the Half-Power Bandwidth
[ bw ] = half_power_bandwidth(f, m_ft);
disp(['Half-Power Bandwidth: ', num2str(bw), ' Hz']);

%% Modulates the message with the carrier
mc = c .* m;

%% The fifth plot : Modulated Message Time Plot
[time, modulated_message] = delimit_window(t, mc, 90*micro, 110*micro); % Delimits the Plot Window

figure(5)
title_text = sprintf('Mensagem m(t) modulada com portadora c(t), fc = %.1f MHz', fc/mega);
file_path = sprintf('../plots/05_fc_%d_kHz_modulated_message_time_plot.png', fc/1000);
plot(time./micro, modulated_message, 'Color', purple, 'LineWidth', 1.5)
xlabel('Tempo (\mus)');
ylabel('Amplitude');
title(title_text);
axis([90 110]);
grid on;
print(file_path, '-dpng');

%% Performs the Fourrier Transform of the Modulated Message
mc_ft = fftshift(fft(mc))*Ts; % Obtain samples of the fourier transform


%% The sixth plot : Modulated Message Frequency Plot
% Delimits the Plot Window
[frequency, mod_m_freq_spectrum] = delimit_window(f, mc_ft, -5.0*mega, 5.0*mega);

figure(6)
title_text = sprintf('Espectro da mensage m(t) modulada com portadora c(t), fc = %.1f MHz', fc/mega);
file_path = sprintf('../plots/06_fc_%d_kHz_modulated_message_freq_plot.png', fc/1000);
plot(frequency./mega, abs(amp_normalize(mod_m_freq_spectrum)), 'Color', purple, 'LineWidth', 1.5);
xlabel('Frequência (MHz)');
ylabel('Magnitude Normalizada');
title(title_text);
axis([-5.0 +5.0]);
grid on;
print(file_path, '-dpng');


%% Demodulates the message with the carrier and corrects amplitude
% Corrects for the amplitude of the carrier and scales the demodulated message by 2
% Due to the fact that the raw demodulated message is 0.5*(A^2)*[*m(t) + cos(2*pi*(2*fc)*t)]
correction_factor = (2/A^2);
dmc = c .* mc * correction_factor; % demodulates and applies the correction factor

dmc_ft = fftshift(fft(dmc))*Ts; % Obtain samples of the fourier transform 

%% The seventh plot : Demodulated Message Frequency Plot
% Delimits the Plot Window
[frequency, dmod_m_freq_spectrum] = delimit_window(f, dmc_ft, -6.0*mega, 6.0*mega);

figure(7)
title_text = sprintf('Espectro da mensage demodulada dm(t), fc = %.1f MHz', fc/mega);
file_path = sprintf('../plots/07_fc_%d_kHz_demodulated_message_freq_plot.png', fc/1000);
plot(frequency./mega, abs(amp_normalize(dmod_m_freq_spectrum)), 'Color', purple, 'LineWidth', 1.5);
xlabel('Frequência (MHz)');
ylabel('Magnitude Normalizada');
title(title_text);
axis([-6.0 +6.0]);
grid on;
print(file_path, '-dpng');

%% Defines the Low Pass Filter (a Rect filter in the frequency domain)
rect_filter = f*0; % Defines the filter and initialize it with 0s

rect_filter(find(f >= -2.0*mega & f <= 2.0*mega)) = 1; % Defines the passband

%% The eighth plot : Demodulated Message + LPF Frequency Plot
% Delimits the Plot Window
[frequency, filter_freq_spectrum] = delimit_window(f, rect_filter, -6.0*mega, 6.0*mega);

figure(8);
title_text = sprintf('Espectro da mensagem demodulada e do filtro, fc = %.1f MHz', fc/mega);
file_path = sprintf('../plots/08_fc_%d_kHz_demodulated_message_and_filter_freq_plot.png', fc/1000);
plot(frequency./mega, abs(amp_normalize(dmod_m_freq_spectrum)), 'Color', purple, 'LineWidth', 1.5); % Plot signal
hold on;
plot(frequency./mega, abs(amp_normalize(filter_freq_spectrum)), 'Color', magenta, 'LineWidth', 1.5); % Plot filter in red

% Add labels, title and scale
xlabel('Frequência (MHz)');
ylabel('Magnitude Normalizada');
legend('Mensagem demodulada', 'Filtro Passa Baixa');
title(title_text);
axis([-6.0 +6.0]);
grid on;
print(file_path, '-dpng');
hold off;

%% Recovers the message applying the filter to the demodulated message
recovered_message_ft = dmc_ft .* rect_filter; % applies the low pass filter

rec_message = ifft(fftshift(recovered_message_ft./Ts)); % uses ifft to recover the samples of the message

%% The nineth plot : Demodulated Message + LPF Frequency Plot
% Delimits the Plots Window
[time, recovered_message] = delimit_window(t, rec_message, 90*micro, 110*micro);
[time, original_message] = delimit_window(t, m, 90*micro, 110*micro);

figure(9)
title_text = sprintf('Mensagem demodulada após FPB, fc = %.1f MHz', fc/mega);
file_path = sprintf('../plots/09_fc_%d_kHz_demodulated_message_post_LPF_time_plot.png', fc/1000);
plot(time./micro, real(recovered_message), 'Color', purple, 'LineWidth', 1.5);
hold on;
plot(time./micro, real(original_message), 'Color', magenta);
xlabel('Tempo (\mus)');
ylabel('Amplitude');
title(title_text);
legend('Mensagem recuperada dm(t)', 'Messagem original m(t)');
axis([90 110]);
grid on;
print(file_path, '-dpng');
hold off;

%% Calculates the similarity between the two messages
corr_matrix = corrcoef(m, rec_message);
similarity = corr_matrix(1, 2);

disp(['Correlation Coefficient between m(t) and recovered dm(t): ', num2str(real(similarity))]);
