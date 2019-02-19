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
function outputoverlay(resultfolder, name, v1, v2, u, seg, secx, secy)

    % Get size of image.
    m = secx(end)-secx(1)+1;
    n = secy(end)-secy(1)+1;

    % Compute mean flow in region.
    meanv1 = mean(v1, 3) .* ~(seg > 0);
    meanv2 = mean(v2, 3) .* ~(seg > 0);
    meanv1 = cat(2, zeros(m, n), meanv1(secx, secy));
    meanv2 = cat(2, zeros(m, n), meanv2(secx, secy));
    
    meanv1 = meanv1(11:end-10 , 11:end-10);
    meanv2 = meanv2(11:end-10 , 11:end-10);
    col = flowcolour(cat(3, meanv1, meanv2), 10);

    % Corrected transparency (alternatively use imfuse).
    folder = fullfile(resultfolder, name);
    mkdir(folder);
    
    % Create frame.
    h = figure(1);
    alph = ~(seg(secx, secy) > 0);
    alph = cat(2, zeros(m, n), 0.5 * alph);
    alph(1:10, :) = 1;
    alph(end-9:end, :) = 1;
    alph(:, 1:10) = 1;
    alph(:, end-9:end) = 1;
    for l=1:size(u, 3)
        cla;
        img = u(secx, secy, l);
        img = cat(2, img, img);
        imagesc(img);
        hold on;
        a = imagesc(col);
        set(a, 'AlphaData', alph);
        colormap gray;
        axis image;
        axis off;
        truesize(h, [m, 2 * n]);
        caxis([0, 1]);
        export_fig(h, fullfile(folder, sprintf('%.3i.png', l)), '-png', '-a1', '-transparent', '-native');
    end
    close(h);
end