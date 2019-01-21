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
classdef(Sealed) vecof3dl2tv < pdproblem
    %VECOF3DL2TV Computes a vector field v for an image sequence f.
    %
    %   v is a matrix of size [n*m*t, 2].
    %   y is matrices of size [n*m*t, 4].
    %   alpha > 0 is a scalar.
    %   Dx, Dy, and Dt are matrices of size [n*m*t, n*m*t].
    %   fx, fy, and ft are vectors of length n*m*t.
    %   n, m are integers.
    properties(GetAccess = public, SetAccess = private)
        v;
        y;
        alpha, beta;
        Dx, Dy, Dt;
        fx, fy, ft;
        n, m, t;
    end
    
    methods
        function o = vecof3dl2tv(fx, fy, ft, alpha, beta, Dx, Dy, Dt, v, y)
            o.alpha = alpha;
            o.beta = beta;
            o.ft = ft(:);
            o.fx = fx(:);
            o.fy = fy(:);
            
            % Initialise operator K.
            [o.n, o.m, o.t] = size(fx);
            o.Dx = Dx;
            o.Dy = Dy;
            o.Dt = Dt;
            
            % Initialise primal variables.
            o.v = v;
            
            % Initialise dual variables.
            o.y = y;
        end

        function v = eval(o)
            % Compute data term + regularisation term.
            d = sum((o.ft + dot([o.fx, o.fy], o.v, 2)).^2);
            r1 = (o.Dx * o.v(:, 1)).^2 + (o.Dy * o.v(:, 1)).^2;
            r2 = (o.Dx * o.v(:, 2)).^2 + (o.Dy * o.v(:, 2)).^2;
            rt = sum((o.Dt * o.v(:, 1)).^2 + (o.Dt * o.v(:, 2)).^2);
            r = sum(sqrt(r1 + r2));
            v = (d + o.alpha * r + o.beta*rt) / numel(o.fx);
        end
        
        function updatePrimal(o, tau)
            % Apply adjoint of K to dual variables and update primal.
            x = o.v - tau * o.applyAdjoint(o.y);
            
            % Compute proximal mapping.
            b1 = x(:, 1) / tau - o.fx .* o.ft;
            b2 = x(:, 2) / tau - o.fy .* o.ft;
            
            % Compute matrix entries.
            c1 = 1 / tau + o.fx .^2;
            c2 = o.fx .* o.fy;
            c3 = 1 / tau + o.fy .^2;
            denom = (c1 .* c3 - c2 .^2);
            
            x(:, 1) = (b1 .* c3 - c2 .* b2) ./ denom;
            x(:, 2) = (b2 .* c1 - c2 .* b1) ./ denom;
            
            o.v = 2 * x - o.v;
        end
        
        function updateDual(o, sigma)
            % Apply K to primal variables.
            o.y = o.y + sigma * o.applyOperator(o.v);

            % Compute proximal maps.
            len = sqrt(sum((o.y(:, [1, 2, 4, 5]) / o.alpha).^2, 2));
            o.y(:, [1, 2, 4, 5]) = o.y(:, [1, 2, 4, 5]) ./ max(1, len);
            
            o.y(:, [3, 6]) = o.beta * o.y(:, [3, 6]) ./ (o.beta + sigma);
        end
        
        function x = primal(o)
            x = o.v;
        end
        
        function y = dual(o)
            y = o.y;
        end
        
        function y = applyOperator(o, x)
            y(:, 1) = o.Dx * x(:, 1);
            y(:, 2) = o.Dy * x(:, 1);
            y(:, 3) = o.Dt * x(:, 1);
            y(:, 4) = o.Dx * x(:, 2);
            y(:, 5) = o.Dy * x(:, 2);
            y(:, 6) = o.Dt * x(:, 2);
        end
        
        function x = applyAdjoint(o, y)
            x(:, 1) = o.Dx' * y(:, 1) + o.Dy' * y(:, 2) + o.Dt' * y(:, 3);
            x(:, 2) = o.Dx' * y(:, 4) + o.Dy' * y(:, 5) + o.Dt' * y(:, 6);
        end
        
        function [v1, v2] = solution(o)
            v1 = reshape(o.v(:, 1), o.n, o.m, o.t);
            v2 = reshape(o.v(:, 2), o.n, o.m, o.t);
        end
    end 
end