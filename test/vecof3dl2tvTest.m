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
function tests = vecof3dl2tvTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function evalTest(testCase)

f = zeros(20, 30, 10);
[n, m, t] = size(f);

% Compute partial derivatives.
fvec = f(:);
[Dx, Dy, ~] = vecderiv3dc(m, n, t-1, 1, 1, 1);
fx = reshape(Dx * fvec(1:n*m*(t-1)), n, m, t-1);
fy = reshape(Dy * fvec(1:n*m*(t-1)), n, m, t-1);
[~, ~, Dt] = vecderiv3dfw(m, n, t, 1, 1, 1);
ft = reshape(Dt * fvec, n, m, t);
ft = ft(:, :, 1:end-1);

% Initialise primal and dual variables.
v = zeros(n * m * (t-1), 2);
y = zeros(n * m * (t-1), 6);

% Create derivative operators.
[Dx, Dy, Dt] = vecderiv3dfw(m, n, t - 1, 1, 1, 1);

p = vecof3dl2tv(fx, fy, ft, 1, 1, Dx, Dy, Dt, v, y);
verifyEqual(testCase, p.eval(), 0);

end