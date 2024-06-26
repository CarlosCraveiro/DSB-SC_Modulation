% Project: SEL0616 - Communication Principles Final Assignment
% Filename: delimit_window.m
% Description: It delimits a window in the input signal (cuts the signal)
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
% delimit_window.m is free software: you can redistribute it and/or modify
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
function [x_cropped, y_cropped] = delimit_window(x, y, window_start, window_end)
    indexes = find(x >= window_start & x <= window_end);
    x_cropped = x(indexes);
    y_cropped = y(indexes);
end
