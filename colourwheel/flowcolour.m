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
%
% This file incorporates work covered by the following copyright and  
% permission notice:
%
% Copyright 2007, Deqing Sun.
%
%                         All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of the author and Brown University not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.
%
% THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
% INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
% PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR BROWN UNIVERSITY BE LIABLE FOR
% ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
%
% Moreover, this file is based on code shared by Hendrik Dirks.
function img = flowcolour(flow, varargin)

if nargin > 1
    bw = varargin{1};
else
    bw = 5;
end

% Get magnitude.
u = flow(:, :, 1);
v = flow(:, :, 2);
len = sqrt(u.^2+v.^2);

% Histogram adjustment by magnitude.
srad = sort(len(:));

% Threshold.
thresh = srad(round(numel(srad) * 0.99));
u(len>thresh) = u(len > thresh) ./ len(len > thresh) * thresh;
v(len>thresh) = v(len > thresh) ./ len(len > thresh) * thresh;

len = hypot(u, v);
maxlen = max(len(:));
u = u / maxlen;
v = v / maxlen;

% Add colour-coding at boundary.
[m, n] = size(u);

topbottom = repmat(-1:2/(n-1):1,bw,1);
left = -ones(m+2*bw,bw);
right = ones(m+2*bw,bw);
u = [left, [topbottom; u; topbottom], right];

top = -ones(bw,n);
bottom = ones(bw,n);
leftright = repmat((-1:2/(m - 1 + 2 * bw):1)', 1, bw);
v = [leftright, [top; v; bottom], leftright];

% Compute colour.
img = computeColour(u, v);  
