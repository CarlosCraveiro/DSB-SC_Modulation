% Project: SEL0616 - Communication Principles Final Assignment
% Filename: half_power_bandwidth.m
% Description: It computes the half power bandwidth of a given signal
%
% Author(s): Carlos Craveiro - USP ID 12547187
%            Ivan Pancheniak - USP ID 12624224
%            Beatriz Aimee   - USP ID 12547934
%
% Created:  2024-24-06
% Modified: 2024-24-06
%
% This file is part of the SEL0616 - Communication Principles Final Assignment project.
%
% half_power_bandwidth.m is free software: you can redistribute it and/or modify
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
function [ bw ] = half_power_bandwidth(frequency, signal_freq_spectrum)
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
