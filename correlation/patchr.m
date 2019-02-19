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
function c = patchr(v1, v2, m, n)
%PATCHR Returns the correlation within a given patch.
%
%   PATCHR(v1, v2, m, n) takes cell arrays v1 and v2 that represent
%   all patches of a vector field (v1, v2), the image size [m, n], and an
%   odd integer s that represents the patch size, and returns a
%   an image c that is the mean resultant vector length for the patch at a
%   each pixel.
%
%   v1, v2 are cell arrays. Each entry is a matrix of size [s, s].
%   [m, n] is the image size.
%   c is a matrix of size [m, n]. 

% Compute length of mean resultant vector for each patch.
c = cellfun(@(p, q) circ_r(cart2pol(p(:), -q(:))), v1, v2, 'UniformOutput', true);
c = reshape(c, [m, n]);

end