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
function p = pixelsize(groupname, dataset)
%PIXELSIZE Returns the pixel size for a specified dataset.

switch groupname
    case '01_control'
        keySet = {'02_006', '02_009', '02_011', '02_013', '04_002', '04_004', '04_010', '04_015'};
        valueSet = [0.303, 0.217, 0.23, 0.303, 0.188, 0.257, 0.269, 0.298];
    case '02_capu[EY12344]'
        keySet = {'16_005', '16_009', '17_010', '17_014', '17_025', '17_028', '17_031', '17_034', '17_040'};
        valueSet = [0.223, 0.168, 0.24, 0.276, 0.303, 0.284, 0.303, 0.286, 0.225];
    case '03_khc[23]'
        keySet = {'01_011', '04_005', '04_007', '04_011', '04_015', '06_005', '06_009', '06_012', '06_015', '06_022'};
        valueSet = [0.217, 0.257, 0.257, 0.269, 0.12, 0.281, 0.283, 0.265, 0.265, 0.228];
    case '04_grk[26B]-grk[2E12]'
        keySet = {'09_014', '09_022', '09_024'};
        valueSet = [0.214, 0.255, 0.255];
    case '05_capu[EY12344]-khc[17]_DOUBLE'
        keySet = {'10_020', '10_023', '11_031', '11_036', '05_004', '05_003', '05_005', '05_008', '05_020'};
        valueSet = [0.335, 0.412, 0.378, 0.227, 0.22, 0.276, 0.304, 0.28, 0.304];
    case '06_capu[EY12344]-khc[17]_CONTROL'
        keySet = {'10_013', '10_017', '10_027', '11_004', '11_018', '14_007', '01_015', '01_017', '01_021', '01_023'};
        valueSet = [0.265, 0.335, 0.381, 0.374, 0.289, 0.37, 0.303, 0.304, 0.248, 0.244];
    case '07_khc[27]'
        keySet = {'12_014', '12_018', '12_022', '12_029', '12_033', '12_037', '13_008', '13_012', '13_018', '13_022'};
        valueSet = [0.39, 0.297, 0.297, 0.434, 0.303, 0.307, 0.324, 0.299, 0.299, 0.299];
    otherwise
        warning('No pixel size set for sequence: %s/%s.\n', groupname, dataset);
        p = 1;
        return;
end

assert(length(keySet) == length(unique(keySet)));
map = containers.Map(keySet, valueSet);

if(~map.isKey(dataset))
    warning('No pixel size set for sequence: %s/%s.\n', groupname, dataset);
    p = 1;
    return;
else

% Find pixel size for given dataset.
p = map(dataset);

end