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
function [Dx, Dy] = deriv2dc(m, n, hx, hy)
%DERIV2DC Creates central first difference matrices for 2D.
%
%   [Dx, Dy] = DERIV2DC(m, n, hx, hy) takes the number of columns m, the 
%   number of rows n, and spatial scaling parameters hx and hy, and 
%   creates first order central difference matrices Dx and Dy with Neumann
%   boundary conditions.
%
%   The gradient of a matrix f of size [n, m] is then given by
%   
%       grad(f) = cat(3, f*Dx', Dy*f).
%
%   For v = (v1, v2), the adjoint -div(v) is then given by
%   
%       -div(v) = v1*Dx + Dy'*v2.
%
%   m and n are integers.
%   hx and hy are scalars.
%   Dx is a matrix of size [m, m].
%   Dy is a matrix of size [n, n].
%
%   Note that m, n > 1.

v1 = -ones(m, 1) / (2*hx);
v2 = ones(m, 1) / (2*hx);
Dx = spdiags([v1, v2], [-1, 1], m, m);
Dx(1, :) = 0;
Dx(m, :) = 0;

v1 = -ones(n, 1) / (2*hy);
v2 = ones(n, 1) / (2*hy);
Dy = spdiags([v1, v2], [-1, 1], n, n);
Dy(1, :) = 0;
Dy(n, :) = 0;

end