function [ bw ] = half_power_bandwidth(frequency, signal_freq_spectrum)
    % Calculates the Half-power bandwidth
    
    % Find the peak magnitude
    magnitude = abs(signal_freq_spectrum);
    peak_magnitude = max(magnitude);
    
    % Find the -3 dB threshold
    threshold = peak_magnitude / sqrt(2);
      
    % Find indices where the magnitude is greater than or equal to the threshold
    indices = find(magnitude >= threshold);
      
    % Determine the frequency range (half-power bandwidth)
    f_start = frequency(indices(1));
    f_end = frequency(indices(end));
  
    % Calculate the bandwidth
    bw = f_end - f_start;
end
