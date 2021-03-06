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
classdef(Sealed) pdhg
    %PDHG Primal Dual Hybrid Gradient Method.
    %   Encapsulates parameters for PDHGM and provides run method.
    
    properties(GetAccess = public, SetAccess = private)
        tau, sigma;
        termeval;
        termhandle;
        logeval;
        loghandle;
    end
    
    methods
        function o = pdhg(tau, sigma, termeval, termhandle, logeval, loghandle)
        % PDHG Takes parameters tau, sigma, and options.
        %   tau and sigma are scalars > 0.
        %   termeval an integer specifying the check for termination.
        %   termhandle a function handle implementing a termination criterion.
        %   logeval an integer specifying the logging interval in iterations.
        %   loghandle a function handle called every logeval iterations.
        %
        %   Both termhandle and loghandle are functions of the form
        %   x(iter, p, pprev, tau, sigma), where iter is the current 
        %   iteration, p and pprev are of type pdproblem, and tau, sigma
        %   are parameters.
            o.tau = tau;
            o.sigma = sigma;
            o.termeval = termeval;
            o.termhandle = termhandle;
            o.logeval = logeval;
            o.loghandle = loghandle;
        end
        
        function stats = run(o, p)
        % RUN Runs the algorithm for a specified problem instance.
        %   p is of type pdproblem.
            tic;
            iter = 1;
            term = false;
            while(~term)
                % Save previous iterate.
                pprev = p.copy;
                
                % Update primal variables.
                p.updatePrimal(o.tau);

                % Update dual variables.
                p.updateDual(o.sigma);

                % Call log function.
                if(mod(iter, o.logeval) == 0)
                    o.loghandle(iter, p, pprev, o.tau, o.sigma);
                end
                
                % Call term function.
                if(mod(iter, o.termeval) == 0)
                    term = o.termhandle(iter, p, pprev, o.tau, o.sigma);
                end
                
                % Increase iteration count.
                iter = iter + 1;
            end
            % Save stats.
            stats.t = toc;
            stats.iter = iter - 1;
            
            % Plot summary.
            fprintf('Terminated after %i iterations.\n', iter - 1);
            fprintf('Elapsed time is %g seconds.\n', stats.t);
            fprintf('Approx. %g seconds/iteration.\n', stats.t / (iter - 1));
            fprintf('Final energy is E(x)=%g.\n', p.eval());
        end
    end
end