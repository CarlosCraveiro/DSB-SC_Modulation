% Project: SEL0616 - Communication Principles Final Assignment
% Filename: amp_normalize.m
% Description: It normalizes the signal with its absolute maximum value
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
% amp_normalize.m is free software: you can redistribute it and/or modify
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
function [ vector_normalized ] = amp_normalize( vector )
    vector_normalized = vector./max(abs(vector));
end
