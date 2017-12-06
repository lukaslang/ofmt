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
function tests = derivperformanceTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function derivative2dTest(testCase)

ntests = 100;

f = ones(500, 250);
[n, m] = size(f);

[Dx, Dy] = deriv2dfw(m, n, 1, 1);
tic;
for k=1:ntests
    x = f * Dx';
    y = Dy * f;
end
toc;

f = ones(n*m, 1);
[Dx, Dy] = vecderiv2dfw(m, n, 1, 1);
tic;
for k=1:ntests
    x = Dx * f;
    y = Dy * f;
end
toc;

end

function derivative3dTest(testCase)

ntests = 100;

f = ones(512, 256, 10);
[n, m, t] = size(f);

[Dx, Dy] = deriv2dfw(m, n, 1, 1);
tic;
for k=1:ntests
    for l=1:t
        x = f(:, :, l) * Dx';
        y = Dy * f(:, :, l);
    end
    z = cat(3, f(:, :, 2:end) - f(:, :, 1:end-1), zeros(n, m));
end
toc;

f = ones(n*m*t, 1);
[Dx, Dy, Dt] = vecderiv3dfw(m, n, t, 1, 1, 1);
tic;
for k=1:ntests
    x = Dx * f;
    y = Dy * f;
    z = Dt * f;
end
toc;

end