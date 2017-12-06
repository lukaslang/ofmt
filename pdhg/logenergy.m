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
function logenergy(iter, p)
%LOGENERGY Plots the energy of the current iterate.
%
%   LOGENERGY(iter, p) takes the current iteration iter and an instance of
%   type pdproblem and logs the current energy.
%
%   iter is a positive integer.
%   p of type pdproblem.

% Log current iteration and energy.
fprintf('Iter=%i, E(x)=%g.\n', iter, p.eval());

end