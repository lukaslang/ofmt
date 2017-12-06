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
function [Dx, Dy, Dz] = vecderiv3dfw(m, n, t, hx, hy, hz)
%VECDERIV3DFW Creates vectorised forward first difference matrices for 3D.
%
%   [Dx, Dy, Dz] = VECDERIV3DFW(m, n, t, hx, hy, hz) takes the number of 
%   columns m, the number of rows n, and the number of time steps t, and 
%   spatial scaling parameters hx, hy, and hz, and creates first order 
%   forward difference matrices Dx Dy, and Dz with Neumann boundary 
%   conditions.
%
%   The gradient of a matrix f in vector form, i.e. of size [n*m*t, 1], is 
%   then given by
%   
%       grad(f) = [Dx*f, Dy*f, Dz*f].
%
%   For v = (v1, v2, v3), the adjoint -div(v) is then given by
%   
%       -div(v) = Dx'*v1 + Dy'*v2 + Dz'*v3.
%
%   m, n, and t are integers.
%   hx, hy, and hz are scalars.
%   Dx, Dy, and Dz are matrices of size [n*m*t, n*m*t].
%
%   Note that m, n, t > 1.

v1 = -ones(t, 1) / hz;
v1(t) = 0;
v2 = ones(t, 1) / hz;
Dz = spdiags([v1, v2], [0, 1], t, t);
Dz = kron(Dz, speye(m*n));

[Dx, Dy] = vecderiv2dfw(m, n, hx, hy);
Dx = kron(speye(t), Dx);
Dy = kron(speye(t), Dy);

end