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
classdef pdproblem < handle
    %PDPROBLEM An abstract problem class.
    %   Provides a type for an abstract primal-dual saddle point problem
    %   of the form
    %
    %       min_x f(Kx) + g(x),
    %
    %   where K is a continuous linear operator, f and g are bounded,
    %   convex, and lower-semicontinuous functions.
    methods(Abstract)
        v = eval(obj);
        % EVAL Returns f(Kx) + g(x).
        
        updatePrimal(obj, tau);
        % UPDATEPRIMAL Performs update step to primal variables.
        
        updateDual(obj, sigma);
        % UPDATEDUAL Performs update step to dual variables.
        
        v = solution(obj);
        % SOLUTION Returns the current solution in the expected dimensions.
    end
end