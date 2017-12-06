% Copyright 2017 Lukas Lang
%
% This file is part of JSME.
%
%    JSME is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    JSME is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with JSME.  If not, see <http://www.gnu.org/licenses/>.
classdef(Sealed) denoise3dl2tv < pdproblem
    %DENOISE3DL2TV Implements image denoising with L2 data term and TV
    %regularisation term.
    
    properties(GetAccess = public, SetAccess = private)
        f, fb, fdelta;
        alpha, beta;
        Dx, Dy, Dt;
        y;
        n, m, t;
    end
    
    methods
        function o = denoise3dl2tv(fdelta, alpha, beta)
        %DENOISE3DL2TV Takes matrix f of size [n, m and a scalar alpha > 0.
            o.fdelta = fdelta(:);
            o.alpha = alpha;
            o.beta = beta;
            o.f = fdelta(:);
            o.fb = fdelta(:);
            
            % Initialise operator K.
            [o.n, o.m, o.t] = size(fdelta);
            [o.Dx, o.Dy, o.Dt] = vecderiv3dfw(o.m, o.n, o.t, 1, 1, 1);
            
            % Initialise dual variables.
            o.y = zeros(o.n*o.m*o.t, 3);
        end
        
        function v = eval(o)            
            % Compute partial derivatives.
            fx = o.Dx * o.f;
            fy = o.Dy * o.f;
            ft = o.Dt * o.f;
            tv = sum(sqrt(fx.^2 + fy.^2));
            % Compute data term + regularisation term.
            v = (0.5 * sum((o.f - o.fdelta).^2) + o.alpha * tv + o.beta * sum(ft.^2)) / numel(o.fdelta);
        end
        
        function updatePrimal(o, tau)
            fold = o.f;
            % Apply adjoint of K to dual variables.
            k = o.Dx' * o.y(:, 1) + o.Dy' * o.y(:, 2) + o.Dt' * o.y(:, 3);
            o.f = (o.f - tau * (k - o.fdelta)) / (1 + tau);
            o.fb = 2 * o.f - fold;
        end
        
        function updateDual(o, sigma)
            % Apply K to primal variables.
            o.y(:, 1) = o.y(:, 1) + sigma * o.Dx * o.fb;
            o.y(:, 2) = o.y(:, 2) + sigma * o.Dy * o.fb;
            len = sqrt((o.y(:, 1)/o.alpha).^2 + (o.y(:, 2)/o.alpha).^2);
            o.y(:, 1) = o.y(:, 1) ./ max(1, len);
            o.y(:, 2) = o.y(:, 2) ./ max(1, len);
            
            o.y(:, 3) = o.y(:, 3) + sigma * o.Dt * o.fb;
            o.y(:, 3) = o.beta * o.y(:, 3) ./ (o.beta + sigma);
        end
        
        function v = solution(o)
            v = reshape(o.f, o.n, o.m, o.t);
        end
    end
end