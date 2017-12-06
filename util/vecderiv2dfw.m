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
function [Dx, Dy] = vecderiv2dfw(m, n, hx, hy)
%VECDERIV2DFW Creates vectorised forward first difference matrices for 2D.
%
%   [Dx, Dy] = VECDERIV2DFW(m, n, hx, hy) takes the number of columns m, the 
%   number of rows n, and spatial scaling parameters hx and hy, and 
%   creates first order forward difference matrices Dx and Dy with Neumann
%   boundary conditions.
%
%   The gradient of a matrix f in vector form, i.e. size [n*m, 1] is then 
%   given by
%   
%       grad(f) = [Dx*f, Dy*f].
%
%   For v = (v1, v2), the adjoint -div(v) is then given by
%   
%       -div(v) = Dx'*v1 + Dy'*v2.
%
%   m and n are integers.
%   hx and hy are scalars.
%   Dx is a matrix of size [n*m, n*m].
%   Dy is a matrix of size [n*m, n*m].
%
%   Note that m, n > 1.

[Dx, Dy] = deriv2dfw(m, n, hx, hy);
Dx = kron(Dx, speye(n));
Dy = kron(speye(m), Dy);

end