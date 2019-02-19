function img = flowToColorV2noBoundary(flow, varargin)

%  flowToColor(flow,boundaryThickness,maxFlow) flowToColor color codes flow field,
%  displayes a boundary layer of boundaryThickness pixels, normalize
%  based on specified value, 
%
%  flowToColor(flow,boundaryThickness) flowToColor color codes flow field,
%  displayes a boundary layer of boundaryThickness pixels, normalize
%  based on maximum flow present otherwise 
% 
%  flowToColor(flow) flowToColor color codes flow field, normalize
%  based on maximum flow present otherwise 
%
%
%
%   FlowToColorV2:
%   2015-06-29: Dynamic scaling of vector length
%   2015-06-16: Added color-coded boundary layer
%   Contact: Hendrik Dirks, hendrik.dirks@wwu.de
%
%
%   
%   According to the c++ source code of Daniel Scharstein 
%   Contact: schar@middlebury.edu
%
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-10-31 18:33:30 (Wed, 31 Oct 2006) $

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

UNKNOWN_FLOW_THRESH = 1e9;
UNKNOWN_FLOW = 1e10;% 

if size(flow,3) ~= 2
    error('flowToColor: image must have two bands');    
end;    

u = flow(:,:,1);
v = flow(:,:,2);

rad = sqrt(u.^2+v.^2);

%automatic adjustment
%analyse histogram of vector length
srad = sort(rad(:));

%choose threshold such that only 1% of flow vectors above
thresh = srad(round(numel(srad)*0.99));

if nargin > 1
    bw = varargin{1};
else
    bw = 5;
end


u(rad>thresh) = u(rad>thresh)./rad(rad>thresh)*thresh;
v(rad>thresh) = v(rad>thresh)./rad(rad>thresh)*thresh;

normV = sqrt(u.^2+v.^2);
maxNormV = max(normV(:));

u = u / maxNormV;
v = v / maxNormV;

% compute color
img = computeColour(u, v);  
