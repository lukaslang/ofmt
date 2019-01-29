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
function outlier = isoutlier(groupname, dataset)
%ISOUTLIER Returns true if dataset is marked as outlier.

switch groupname
    case '05_capu[EY12344]-khc[17]_DOUBLE'
        keySet = {'11_036'};
        valueSet = true;
    case '07_khc[27]'
        keySet = {'12_037'};
        valueSet = true;
    otherwise
        outlier = false;
        return;
end

assert(length(keySet) == length(unique(keySet)));
map = containers.Map(keySet, valueSet);

outlier = map.isKey(dataset);

end