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
function pdres = pdresidual(p, pprev, tau, sigma)
%PDRESIDUAL Returns the primal-dual residual.
%
%   PDRESIDUAL(p, pprev, tau, sigma) takes two instances of type pdproblem,
%   parameters tau and sigma, and returns the primal-dual residual.
%
%   p, ppref are of type pdproblem.
%   tau, sigma > 0 are scalars.
%
%   See Goldstein, Esser, and Baraniuk. Adaptive primal-dual hybrid
%   gradient methods for saddle-point problems, arXiv preprint 
%   arXiv:1305.0546, 2013.

pres = (pprev.primal - p.primal)/tau - p.applyAdjoint(pprev.dual - p.dual);
dres = (pprev.dual - p.dual)/sigma - p.applyOperator(pprev.primal - p.primal);
pdres = (sum(abs((pres(:)))) + sum(abs((dres(:))))) / numel(p.primal);

end