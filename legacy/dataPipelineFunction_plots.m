function result = dataPipelineFunction_plots(folder)
    result = 1;

    %% load data

    if (~exist([folder,'mat']))
        mkdir([folder,'mat']);
    end
    if (~exist([folder,'images']))
        mkdir([folder,'images']);
    end
    
    if (exist([folder,'mat',filesep,'rawData.mat'],'file'))
        load([folder,'mat',filesep,'rawData.mat'],'imageStack');
    else
        folderContent = dir([folder,filesep,'*.tif']);

        numImages = 1;
        for i=1:numel(folderContent)
            if (~folderContent(i).isdir)
                rawImage = imread([folder,folderContent(i).name]);
                rawImage = im2double(rawImage);

                imageStack(:,:,numImages) = rawImage;

                numImages = numImages + 1;
            end
        end

        save([folder,'mat',filesep,'rawData.mat'],'imageStack');

        clear numImages rawImage i folderContent;
    end

    %%
    if (exist([folder,'mat',filesep,'results.mat'],'file'))
        load([folder,'mat',filesep,'results.mat'],'u')
        load([folder,'mat',filesep,'results.mat'],'v')
    else
        error('No results found in folder %s. Run flow computation first.', folder);
    end

    
    %% Plot and output flow.
    for i=1:size(imageStack,3)
        singleField = squeeze(v(:,:,i,:));

        %figure(1);imagesc(imageStack(:,:,i),[0,1]);axis image;colormap(gray);
        %figure(2);imagesc(u(:,:,i),[0,1]);axis image;colormap(gray);
        %figure(3);imagesc(flowToColorV2(singleField));axis image;colormap(gray);

        exportImage(imageStack(:,:,i),[folder,'images',filesep,'input',num2str(i),'.png'],'colormap',gray);
        exportImage(u(:,:,i),[folder,'images',filesep,'reconstructed',num2str(i),'.png'],'colormap',gray);
        exportImage(flowToColorV2(singleField),[folder,'images',filesep,'flowField',num2str(i),'.png'],'colormap',gray);
    end

    %%
    clear alpha beta gamma i mainIterations mainJoint numFrames singleField temporalSmoothness verbose; 
    %% go through slices and normalize
    mean1 = u(:,:,1);mean1 = sum(mean1(:)) / numel(mean1);

    for i=2:size(u,3)
        meanCurr = u(:,:,i);meanCurr = sum(meanCurr(:)) / numel(meanCurr);
        u(:,:,i) = u(:,:,i) * mean1/meanCurr;
    end
    clear mean1 meanCurr;

    %% create averaged flow field
    averagedMotionField = squeeze(sum(v,3) / size(v,3));

    exportImage(flowToColorV2(averagedMotionField),[folder,'images',filesep,'averagedFlowField.png']);


    %% segmentationMapOuter

    segmentationMapOuter = im2double((imread([folder,'images',filesep,'segmentationMap.png'])));
    %segmentationMapOuter(segmentationMapOuter<1) = 0;
    segmentationMapOuter = imresize(segmentationMapOuter,[size(u,1),size(u,2)]);
    segmentationMapOuter(segmentationMapOuter>0) = 1;
    segmentationMapOuter = segmentationMapOuter == 0;

    figure(10);imagesc(segmentationMapOuter);colormap(gray);axis image;title('SegmentationMap')

    %% create images with velocities which correlate to intensities above a certain threshold
    threshold = 0.2;

    vCopy = v;

    for i=2:size(v,3)
        imageMask = ((u(:,:,i) >= threshold)) .* (segmentationMapOuter == 1);
        imageMask = (imageMask ==0);

         vPart1 = v(:,:,i,1);vPart1(imageMask) = 0;
         vPart2 = v(:,:,i,2);vPart2(imageMask) = 0;

        vCopy(:,:,i,1) = vPart1;
        vCopy(:,:,i,2) = vPart2;
    end

    %%


    averagedMotionField = squeeze(sum(vCopy,3));

    averagedMotionField1 = averagedMotionField(:,:,1);
    averagedMotionField2 = averagedMotionField(:,:,2);

    averagedMotionField1(segmentationMapOuter == 0) = 0;
    averagedMotionField2(segmentationMapOuter == 0) = 0;

    averagedMotionFieldNew = cat(3,averagedMotionField1,averagedMotionField2);

    figure(801);imagesc(flowToColorV2(averagedMotionFieldNew,5,1));axis image;
    exportImage(flowToColorV2(averagedMotionFieldNew,5,1),[folder,'images',filesep,'flowAveragedOutside.png']);


    %% create images with velocities which correlate to intensities above a certain threshold

    vCopy = v;

    segmentationMapInner = segmentationMapOuter == 0;
    figure(10);imagesc(segmentationMapInner);colormap(gray);axis image;title('SegmentationMap')

    for i=1:size(v,3)
        vPart1 = v(:,:,i,1);vPart1(segmentationMapInner==0) = 0;
        vPart2 = v(:,:,i,2);vPart2(segmentationMapInner==0) = 0;

        vCopy(:,:,i,1) = vPart1;
        vCopy(:,:,i,2) = vPart2;
    end

    averagedMotionField = squeeze(sum(vCopy,3) / size(vCopy,3));

    %%

    exportImage(flowToColorV2(averagedMotionField),[folder,'images',filesep,'averagedFlowFieldInner.png']);

    %%
    %transform averaged field into polar coordinates

    component1 = averagedMotionField(:,:,1);
    component1 = component1(segmentationMapInner);

    component2 = averagedMotionField(:,:,2);
    component2 = -component2(segmentationMapInner);%we have to invert the second component!


    [angle,~] = cart2pol(component1,component2);

    angle(angle==0) = [];


    angle = mod(angle, 2*pi);
    angle = rad2deg(angle);

    [angleHist,angleCenters] = hist(angle,360);

    figure(5);polar(1,1);
    print([folder,'images',filesep,'angleOverview'],'-dpng','-r0')

    meanAngle = mean(angle);
    % Modified since we don't have Curve Fitting Toolbox.
    %figure(6);plot(angleCenters,smooth(angleHist),angleCenters,ones(1,numel(angleCenters))*meanAngle);
    figure(6);plot(angleCenters,angleHist,angleCenters,ones(1,numel(angleCenters))*meanAngle);



    print([folder,'images',filesep,'anglePlot'],'-dpng','-r0')

    figure(7);rose(degtorad(angle),100);
    print([folder,'images',filesep,'angleRoseDiagramFine'],'-dpng','-r0')

    figure(8);rose(degtorad(angle),8);
    print([folder,'images',filesep,'angleRoseDiagramCoarse'],'-dpng','-r0')
    
    
    [tout,rout] = rose(degtorad(angle),100);
    figure(9);polar(tout,rout/numel(angle))
    print([folder,'images',filesep,'angleRoseDiagramFineRelative'],'-dpng','-r0')
    
    [tout,rout] = rose(degtorad(angle),8);
    figure(10);polar(tout,rout/numel(angle))
    print([folder,'images',filesep,'angleRoseDiagramCoarseRelative'],'-dpng','-r0')
    

end