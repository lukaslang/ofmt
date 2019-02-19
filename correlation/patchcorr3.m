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
function c = patchcorr3(v1, v2, m, n, s)
%PATCHCORR3 Returns the correlation within a given patch.
%
%   PATCHCORR3(v1, v2, m, n, s) takes cell arrays v1 and v2 that represent
%   all patches of a vector field (v1, v2), the image size [m, n], and an
%   odd integer s that represents the patch size, and returns a
%   correlation measure c.
%
%   v1, v2 are cell arrays. Each entry is a matrix of size [s, s].
%   [m, n] is the image size.
%   s is an odd integer >= 1.
%   c is a matrix of size [m, n]. 

% Normalise each vector.
len = cellfun(@(x, y) sqrt(x.^2 + y.^2 + 1), v1, v2, 'UniformOutput', false);
v1 = cellfun(@(x, y) x ./ y, v1, len, 'UniformOutput', false);
v2 = cellfun(@(x, y) x ./ y, v2, len, 'UniformOutput', false);

% Extract centre value for each patch.
c = (s - 1) / 2 + 1;
c1 = cellfun(@(x) repmat(x(c, c), [s, s]), v1, 'UniformOutput', false);
c2 = cellfun(@(x) repmat(x(c, c), [s, s]), v2, 'UniformOutput', false);

avgangerr = @(a, b, c, d) mean(a(:).*c(:) + b(:).*d(:) + 1);
c = cellfun(@(u, v, p, q) avgangerr(u, v, p, q), c1, c2, v1, v2, 'UniformOutput', true);
c = reshape(c, [m, n]);

end