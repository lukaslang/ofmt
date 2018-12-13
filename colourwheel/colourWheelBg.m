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
function cw = colourWheelBg
%COLOURWHEELBG Creates a colour wheel with white background.

% Create colour wheel.
cw = colourWheel;
[m, n, ~] = size(cw);

% Create and normalise coordiantes.
[X, Y] = meshgrid(1:m, 1:n);
X = 2*X / m - 1;
Y = 2*Y / m - 1;

% Get indices.
idx = sqrt(X.^2 + Y.^2) >= 0.96;
cwr = cw(:, :, 1);
cwg = cw(:, :, 2);
cwb = cw(:, :, 3);
cwr(idx) = 255;
cwg(idx) = 255;
cwb(idx) = 255;
cw = cat(3, cwr, cwg, cwb);

end