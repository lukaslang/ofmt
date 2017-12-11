function result = dataPipelineFunction_flow(folder)
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
    verbose = 1; %can be 0,1 or 2
    alpha = 0.005; %regularizer u
    beta = 0.01; % regularizer v
    gamma = 0.1; % optical flow weight
    temporalSmoothness = 0.5; 
    mainIterations = 2;
    numFrames = size(imageStack,3);

    imageStack = imageStack(:,:,1:numFrames);
    
    if (exist([folder,'mat',filesep,'results.mat'],'file'))
        load([folder,'mat',filesep,'results.mat'],'u')
        load([folder,'mat',filesep,'results.mat'],'v')
    else
        mainJoint = jointModelLargeScale(imageStack ,alpha,beta,gamma,'temporalSmoothness',temporalSmoothness,'verbose',verbose);
        mainJoint.verbose = 1;

        mainJoint.init;

        mainJoint.numMainIt = mainIterations;
        %% run


        mainJoint.run;

        %%
        u = mainJoint.u;
        v = mainJoint.v;



        %% save data
        save([folder,'mat',filesep,'results.mat'],'u','v')
        
        for i=1:size(imageStack,3)
            singleField = squeeze(v(:,:,i,:));

            %figure(1);imagesc(imageStack(:,:,i),[0,1]);axis image;colormap(gray);
            %figure(2);imagesc(u(:,:,i),[0,1]);axis image;colormap(gray);
            %figure(3);imagesc(flowToColorV2(singleField));axis image;colormap(gray);

            exportImage(imageStack(:,:,i),[folder,'images',filesep,'input',num2str(i),'.png'],'colormap',gray);
            exportImage(u(:,:,i),[folder,'images',filesep,'reconstructed',num2str(i),'.png'],'colormap',gray);
            exportImage(flowToColorV2(singleField),[folder,'images',filesep,'flowField',num2str(i),'.png'],'colormap',gray);

        end
    
    end
end