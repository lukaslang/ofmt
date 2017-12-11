function exportImage( image,filepath,varargin )
    if (nargin > 2)
        for i=1:2:nargin-2
            if (strcmp(varargin{i},'colormap'))
                cm = varargin{i+1};
            end
            
            if (strcmp(varargin{i},'limits'))
                lm = varargin{i+1};
            end
            
            if (strcmp(varargin{i},'v1'))
                v1 = varargin{i+1};
            end
            if (strcmp(varargin{i},'v2'))
                v2 = varargin{i+1};
            end
            if (strcmp(varargin{i},'vScale'))
                vScale = varargin{i+1};
            end
            if (strcmp(varargin{i},'velocityParam'))
                velocityParam = varargin{i+1};
            end
            if (strcmp(varargin{i},'limitImagesc'))
                limitImagesc = varargin{i+1};
            end
            if (strcmp(varargin{i},'type'))
                type = varargin{i+1};
            end
            if (strcmp(varargin{i},'plotArgument'))
                plotArgument = varargin{i+1};
            end
            if (strcmp(varargin{i},'plotBackground'))
                plotBackground = varargin{i+1};
            end
            if (strcmp(varargin{i},'showColorbar'))
                showColorbar = varargin{i+1};
            end
            if (strcmp(varargin{i},'imageOverlay'))
                imageOverlay = varargin{i+1};
            end
            if (strcmp(varargin{i},'format'))
                format = varargin{i+1};
            end
            %exporter is able to plot different contours, just paste a cell
            %array with (distanceFunction,color,lineWidth)
            if (strcmp(varargin{i},'plotContour'))
                plotContour = varargin{i+1};
            end
        end
    end
    
    if (~exist('cm','var'))
        cm = jet;
    end

    if (~exist('lm','var'))
        lm = [0 0 size(image,2) size(image,1)];
        %lm = [0 0 500 500];
    end
    
    if (~exist('velocityParam','var'))
        velocityParam = 'r';
    end
%     if (~exist('limitImagesc','var'))
%         limitImagesc = [0 1];
%     end
    if (~exist('vScale','var'))
        vScale = 1;
    end
    if (~exist('type','var'))
        type = 'image';
    end
    if (~exist('plotArgument','var'))
        plotArgument = '';
    end
    if (~exist('plotBackground','var'))
        plotBackground = 0;
    end
    
    if (~exist('showColorbar','var'))
        showColorbar = 0;
    end
    
    if (~exist('format','var'))
        %format = '-depsc';
        format = '-dpng';
        
    end
    
    thisFigure = figure('visible','off');clf;
    
    if (strcmp(type,'image'))
        if (exist('limitImagesc','var'))
            imagesc(image,limitImagesc);
        else
            imagesc(image);
        end
        axis image;
        set(gca,'TickDir','out');
        colormap(cm);
        if (showColorbar) 
            colorbar('location','South');
        end
    elseif (strcmp(type,'plot'))
        if (sum(abs(plotBackground(:)))~=0)
            imagesc(plotBackground);axis image;hold on;
            colormap(cm);
        end
        plot(image(1,:),image(2,:),plotArgument);
    end
    
    if (exist('plotContour','var'))
        for i=1:numel(plotContour)
            cont = plotContour{i};
            hold on;
            contour(cont{1},[0 0],cont{2},'LineWidth',cont{3});
        end
    end
    
    if (exist('v1','var'))
        %create meshgrid
        x2 = 1:vScale:size(v1,1);
        y2 = 1:vScale:size(v1,2);

        [X2,Y2] = meshgrid(y2,x2);

        hold on;quiver(X2,Y2,v1(x2,y2),v2(x2,y2),0,velocityParam)
        
        set(gca,'XLim',[0,size(image,2)],'YLim',[0,size(image,1)]);
    end
    
    if (exist('imageOverlay','var'))
        hold on;
        h=imshow(cat(3, ones(size(image)), zeros(size(image)), zeros(size(image))));
        set(h,'AlphaData', imageOverlay,'AlphaDataMapping','none');
    end

    %set(gcf, 'Position', lm,'OuterPosition', lm, 'Color', 'w');
    set(gcf, 'Position', lm, 'Color', 'w');
    
    set(gcf, 'PaperPositionMode','auto')   %# WYSIWYG
    set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);

    print(thisFigure,filepath,format,'-r200');
    axis off
    clf;
    close(thisFigure);
end

