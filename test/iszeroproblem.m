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
classdef(Sealed) iszeroproblem < pdproblem
    %ISZEROPROBLEM Implements a dummy problem.
    
    properties(GetAccess = public, SetAccess = private)
        x;
    end
    
    methods
        function obj = iszeroproblem(x)
            obj.x = x;
        end
        
        function v = eval(obj)
            if(obj.x == 0)
                v = 0;
            else
                v = Inf;
            end
        end
        
        function updatePrimal(~, ~)
            % Does nothing.
        end
        
        function updateDual(~, ~)
            % Does nothing.
        end
        
        function v = solution(obj)
            v = obj.x;
        end
    end
end