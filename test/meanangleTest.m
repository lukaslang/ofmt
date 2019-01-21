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
function tests = meanangleTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

angles = [pi/2, -pi/2]';

[mangle, r] = meanangle(angles);
verifyEqual(testCase, mangle, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, r, 0, 'AbsTol', 1e-12);

angles = [pi, pi]';

[mangle, r] = meanangle(angles);
verifyEqual(testCase, mangle, pi, 'AbsTol', 1e-12);
verifyEqual(testCase, r, 1, 'AbsTol', 1e-12);

end