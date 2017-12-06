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
function tests = vecderiv3dfwTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

m = 15;
n = 10;
t = 3;

hx = 1;
hy = 1;
hz = 1;

[Dx, Dy, Dz] = vecderiv3dfw(m, n, t, hx, hy, hz);
verifyEqual(testCase, size(Dx), [n*m*t, n*m*t]);
verifyEqual(testCase, size(Dy), [n*m*t, n*m*t]);
verifyEqual(testCase, size(Dz), [n*m*t, n*m*t]);

end

function calcDerivativeOfZeroMatrixTest(testCase)

% Set parameters.
hx = 1;
hy = 2;
hz = 3;

% Create matrix.
m = 4;
n = 5;
t = 3;
f = ones(n, m, t);

% Create derivative matrices.
[Dx, Dy, Dz] = vecderiv3dfw(m, n, t, hx, hy, hz);

% Compute derivatives.
fx = reshape(Dx*f(:), n, m, t);
fy = reshape(Dy*f(:), n, m, t);
fz = reshape(Dz*f(:), n, m, t);

verifyEqual(testCase, fx, zeros(n, m, t));
verifyEqual(testCase, fy, zeros(n, m, t));
verifyEqual(testCase, fz, zeros(n, m, t));

end

function adjointTest(testCase)

% Set parameter.
hx = 2;
hy = 2;
hz = 2;

% Create matrix.
m = 4;
n = 5;
t = 3;
f = rand(n, m, t);

% Create derivative matrices.
[Dx, Dy, Dz] = vecderiv3dfw(m, n, t, hx, hy, hz);

% Create second matrix.
v = rand(n, m, t);

% Compute derivatives.
fx = Dx*f(:);
fy = Dy*f(:);
fz = Dz*f(:);

% Compute divergence.
vx = Dx'*v(:);
vy = Dy'*v(:);
vz = Dz'*v(:);

% Verify adjoint property.
verifyEqual(testCase, fx' * v(:), f(:)' * vx, 'AbsTol', 1e-15);
verifyEqual(testCase, fy' * v(:), f(:)' * vy, 'AbsTol', 1e-15);
verifyEqual(testCase, fz' * v(:), f(:)' * vz, 'AbsTol', 1e-15);

end