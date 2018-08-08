% Copyright 2017 Lukas Lang
%
% This file is part of OFMT.
%
%    OFMT is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    OFMT is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with OFMT.  If not, see <http://www.gnu.org/licenses/>.
function m = meanangle(v)
%MEANANGLE Computes the circular mean angle.
%
%   m = MEANANGLE(v) takes a vector of angles in the interval
%   [0, 360] and returns the mean circular angle m according to
%
%   https://en.wikipedia.org/wiki/Mean_of_circular_quantities
%
%   See https://rosettacode.org/wiki/Averages/Mean_angle#MATLAB_.2F_Octave
%
%   v is a vector of size [n, 1].
%
%   m is a scalar in [0, 360).

m = angle(mean(exp(1i * v)));

end