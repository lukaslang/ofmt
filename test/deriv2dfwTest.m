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
function tests = deriv2dfwTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

m = 4;
n = 3;
hx = 2;
hy = 2;

Lx = [-1, 1, 0, 0;
       0,-1, 1, 0;
       0, 0,-1, 1;
       0, 0, 0, 0] / hx;

Ly = [-1, 1, 0;
       0,-1, 1;
       0, 0, 0] / hy;
   
[Dx, Dy] = deriv2dfw(m, n, hx, hy);
verifyEqual(testCase, full(Dx), Lx);
verifyEqual(testCase, full(Dy), Ly);

m = 1;
n = 1;
hx = 2;
hy = 2;

Lx = 0 / hx;
Ly = 0 / hy;
   
[Dx, Dy] = deriv2dfw(m, n, hx, hy);
verifyEqual(testCase, full(Dx), Lx);
verifyEqual(testCase, full(Dy), Ly);

end

function calcDerivativeOfZeroMatrixTest(testCase)

% Set parameters.
hx = 2;
hy = 2;

% Create matrix.
m = 4;
n = 5;
f = ones(n, m);

% Create derivative matrices.
[Dx, Dy] = deriv2dfw(m, n, hx, hy);

% Compute derivatives.
fx = f*Dx';
fy = Dy*f;

verifyEqual(testCase, fx, zeros(n, m));
verifyEqual(testCase, fy, zeros(n, m));

end

function adjointTest(testCase)

% Set parameter.
hx = 2;
hy = 2;

% Create matrix.
m = 4;
n = 5;
f = rand(n, m);

% Create derivative matrices.
[Dx, Dy] = deriv2dfw(m, n, hx, hy);

% Create second matrix.
v = rand(n, m);

% Compute derivatives.
fx = f*Dx';
fy = Dy*f;

% Compute divergence.
vx = v*Dx;
vy = Dy'*v;

% Verify adjoint property.
verifyEqual(testCase, fx(:)' * v(:), f(:)' * vx(:), 'AbsTol', 1e-15);
verifyEqual(testCase, fy(:)' * v(:), f(:)' * vy(:), 'AbsTol', 1e-15);

end