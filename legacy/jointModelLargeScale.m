classdef jointModelLargeScale < handle
    %JOINTSUPERRESOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
         numFrames
         dims
         factor
         alpha %weight for tv u
         listAlpha
         beta %weight for tv v
         betaStart %initial weight for tv v
         gamma %parameter for coupling
         numMainIt;
         numMinMainIt;
         imageSequence
         mainU
         mainV
         numWarpTerm
         u %image sequence
         uOld
         v %flow sequence
         vOld
         temporalSmoothness
         doWarping %use warping in the optical flow problem to massively enhance the result
         medianFiltering %use median filtering in the optical flow problem to massively enhance the result
         doGradientConstancy %use gradient constancy in the optical flow problem
         opticalFlowTerm %specifies the type of optical flow term. Either the classical or the warping for large scale
         regularizerTermV %specifies the regularizer for the flow field v, can be 'TV', 'Huber' or 'L2'
         
         tolV
         tolMain
         
         verbose
         
         gtU
         gtV
         
         listK %list of operators, one for each image in imageSequence
    end
    
    methods
        function obj = jointModelLargeScale(imageSequence,alpha,beta,gamma,varargin)
            vararginParser;
            
            obj.imageSequence = imageSequence;
            
            obj.alpha = alpha;
            obj.beta = beta;
            obj.gamma = gamma;
            obj.dims = size(imageSequence);
            obj.dims = obj.dims(1:2);
            
            obj.numMainIt = 3;
            obj.numFrames = size(imageSequence,3);
            
            obj.u = zeros([obj.dims,obj.numFrames]);
            obj.v = zeros([obj.dims,obj.numFrames,2]);
            
            if (exist('gtU','var'))
                obj.gtU = gtU;
            else
                obj.gtU = 0;
            end
            
            if (exist('gtV','var'))
                obj.gtV = gtV;
            else
                obj.gtV = 0;
            end
            
            if (exist('verbose','var'))
                obj.verbose = verbose;
            else
                obj.verbose = 1;
            end
            
            if (exist('betaStart','var'))
                obj.betaStart = betaStart;
            else
                obj.betaStart = 0.01;
            end
            
            if (exist('temporalSmoothness','var'))
                obj.temporalSmoothness = temporalSmoothness;
            else
                obj.temporalSmoothness = 0;
            end
            
            if (exist('doWarping','var'))
                obj.doWarping = doWarping;
            else
                obj.doWarping = 1;
            end
            
            if (exist('medianFiltering','var'))
                obj.medianFiltering = medianFiltering;
            else
                obj.medianFiltering = 1;
            end
            
            if (exist('doGradientConstancy','var'))
                obj.doGradientConstancy = doGradientConstancy;
            else
                obj.doGradientConstancy = 1;
            end
            
            if (exist('tolV','var'))
                obj.tolV = tolV;
            else
                obj.tolV = 1e-6;
            end
            
            if (exist('tolMain','var'))
                obj.tolMain = tolMain;
            else
                obj.tolMain = 1e-5;
            end
			
            if (exist('opticalFlowTerm','var'))
                obj.opticalFlowTerm = opticalFlowTerm;
            else
                obj.opticalFlowTerm = 'warping';
            end
            
            if (exist('regularizerTermV','var'))
                obj.regularizerTermV = regularizerTermV;
            else
                obj.regularizerTermV = 'TV';
            end
            
            
            
            if (exist('listK','var'))
                obj.listK = listK;
            else %create list of identity operators
                
                dims = size(imageSequence);
                nPx = dims(1) * dims(2);
                
                for i=1:size(imageSequence,3)
                    obj.listK{i} = identityOperator(nPx);
                end
            end
            
            if (exist('listAlpha','var'))
                obj.listAlpha = listAlpha;
            else
                for i=1:size(imageSequence,3)
                    obj.listAlpha{i} = alpha;
                end
            end
            
            if (exist('numMinMainIt','var'))
                obj.numMinMainIt = numMinMainIt;
            else
                obj.numMinMainIt = 4;
            end
            
            
        end
        
        function init(obj)
            if (obj.verbose > 0) %show on light verbose level
                disp('Calculating initial u');
            end
            
            %create initial u and v
            initialU = flexBox;
            initialU.params.tryCPP = 1;

            for i=1:obj.numFrames
                %add primal variable for each image
                initialU.addPrimalVar(obj.dims);

                %add data term
                initialU.addTerm(L2dataTermOperator(1,obj.listK{i},obj.imageSequence(:,:,i)),i);

                %add tv term for each primal var
                initialU.addTerm(L1gradientIso(obj.listAlpha{i},obj.dims),i);

                if (i>1 && obj.temporalSmoothness > 0)
                    idOp = speye(prod(obj.dims));
                    
                    initialU.addTerm(L2operator(obj.temporalSmoothness,2,{idOp,-idOp}),[i,i-1]);
                end
            end

            initialU.runAlgorithm;
            
            for j=1:obj.numFrames
                obj.u(:,:,j) = initialU.getPrimal(j);

                if (obj.verbose > 1) %show only on massive verbose level
                    figure(100+j);clf;imagesc(obj.u(:,:,j));axis image;colormap(gray);title(['Initial Image #',num2str(j)]);drawnow;
                end
            end

            clear initialU;

            if (obj.verbose > 0)
                disp('Calculating initial velocity fields')
            end
            
            %% init real motion estimator
            if (strcmp(obj.opticalFlowTerm,'classic'))
                imageDiscretization = 'regularCentral';
            else
                imageDiscretization = 'interpolated';
            end
            
            obj.mainV = motionEstimatorClass(obj.u,obj.tolV,obj.beta,'verbose',obj.verbose,'doGradientConstancy',obj.doGradientConstancy,'medianFiltering',obj.medianFiltering,'doWarping',obj.doWarping,'imageDiscretization',imageDiscretization,'regularizerTerm',obj.regularizerTermV);
            obj.mainV.init;
            
            obj.solveV;
            
            %% create flexBox object for u
            obj.mainU = flexBox;

            obj.mainU.params.showPrimals = 100;
            obj.mainU.params.tol = 1e-6;
            obj.mainU.params.verbose = obj.verbose;
            obj.mainU.params.tryCPP = 1;
    
            %add for each frame data term and tv term
            for i=1:obj.numFrames
                %add primal variable for each image
                obj.mainU.addPrimalVar(obj.dims);

                %add data term
                %obj.mainU.addTerm(L2dataTermOperator(1,identityOperator(numel(obj.imageSequence(:,:,i))),obj.imageSequence(:,:,i)),i);
                obj.mainU.addTerm(L2dataTermOperator(1,obj.listK{i},obj.imageSequence(:,:,i)),i);
                
                %add tv term for each primal var
                obj.mainU.addTerm(L1gradientIso(obj.listAlpha{i},obj.dims),i);
                
                if (i>1 && obj.temporalSmoothness > 0)
                    %idOp = speye(prod(obj.dims));
                    %obj.mainU.addTerm(L2operator(obj.temporalSmoothness,2,{idOp,-idOp}),[i,i-1]);
                end
            end

            %create optical flow term connecting subsequent images; this is only a placeholder. The operators will be overwritten by updateFlexBoxU
            for i=1:obj.numFrames - 1
                obj.mainU.addTerm(L1operatorAniso(obj.gamma,2,{zeroOperator(prod(obj.dims)),zeroOperator(prod(obj.dims))}),[i,i+1]);
                obj.numWarpTerm(i) = numel(obj.mainU.duals);
            end
			
			%assign correct operators to optical flow term
			obj.updateFlexBoxU;
            
            if (obj.verbose > 0)
                disp('Initialization finished');
            end
        end
        
        function solveU(obj)
            obj.mainU.runAlgorithm;

            % extract solution and store into matrix u
            for j=1:obj.numFrames
                obj.u(:,:,j) = obj.mainU.getPrimal(j);
                if (obj.verbose > 1)
                    figure(100+j);imagesc(obj.u(:,:,j),[0,1]);axis image;colormap(gray);drawnow
                end
            end
        end
        
        function updateFlexBoxU(obj)
		
			if (strcmp(obj.opticalFlowTerm,'classic'))
				gradientXf = gradientOperator(obj.dims,1,'discretization','forward');
                gradientXb = gradientOperator(obj.dims,1,'discretization','backward');
                gradientY = (gradientXf.matrix + gradientXb.matrix)/2;
                
                gradientYf = gradientOperator(obj.dims,2,'discretization','forward');
                gradientYb = gradientOperator(obj.dims,2,'discretization','backward');
				gradientX = (gradientYf.matrix + gradientYb.matrix)/2;
				
				nPx = prod(obj.dims);
				
				for j=1:obj.numFrames-1
					 %extract flow field between i and i+1
					singleField1 = reshape(obj.v(:,:,j,1),[nPx,1]);
					singleField2 = reshape(obj.v(:,:,j,2),[nPx,1]);
					
					idOp = speye(prod(obj.dims));
					
					ofOp = idOp + spdiags(singleField1,0,nPx,nPx)*gradientX + spdiags(singleField2,0,nPx,nPx) * gradientY;
					
					obj.mainU.duals{obj.numWarpTerm(j)}.operator{1} = -idOp;
					obj.mainU.duals{obj.numWarpTerm(j)}.operatorT{1} = -idOp';

					obj.mainU.duals{obj.numWarpTerm(j)}.operator{2} = ofOp;
					obj.mainU.duals{obj.numWarpTerm(j)}.operatorT{2} = ofOp';
				end
			else %warping
				for j=1:obj.numFrames-1
					 %extract flow field between i and i+1
					singleField = squeeze(obj.v(:,:,j,:));

					%create warping operator forward and backward
					warp = warpingOperator(obj.dims,singleField);

					%find out of range warps in each of the operators and set the
					%corresponding line in the other operator also to zero
					marker = sum(abs(warp),2) == 0;
					marker = marker > 0;
					
					idOp = -speye(prod(obj.dims));

					warp(marker,:) = 0;
					idOp(marker,:) = 0;
					
					obj.mainU.duals{obj.numWarpTerm(j)}.operator{1} = idOp;
					obj.mainU.duals{obj.numWarpTerm(j)}.operatorT{1} = idOp';

					obj.mainU.duals{obj.numWarpTerm(j)}.operator{2} = warp;
					obj.mainU.duals{obj.numWarpTerm(j)}.operatorT{2} = warp';
				end
			end
        end
        
        function errorConstraint = getErrorInConstraint(obj)
            gradientXf = gradientOperator(obj.dims,1,'discretization','forward');
            gradientXb = gradientOperator(obj.dims,1,'discretization','backward');
            gradientY = gradientXf.matrix + gradientXb.matrix;

            gradientYf = gradientOperator(obj.dims,2,'discretization','forward');
            gradientYb = gradientOperator(obj.dims,2,'discretization','backward');
            gradientX = gradientYf.matrix + gradientYb.matrix;
            
            nPx = prod(obj.dims);
                
            for j=1:obj.numFrames-1
                %extract flow field between i and i+1
                singleField1 = reshape(obj.v(:,:,j,1),[nPx,1]);
                singleField2 = reshape(obj.v(:,:,j,2),[nPx,1]);

                idOp = speye(prod(obj.dims));

                ofOp = idOp + spdiags(singleField1,0,nPx,nPx)*gradientX + spdiags(singleField2,0,nPx,nPx) * gradientY;
                    
                im1 = obj.mainU.getPrimal(j);
                im2 = obj.mainU.getPrimal(j+1);

                errorConstraint(:,:,j) = abs(reshape(-idOp * im1(:) + ofOp * im2(:),size(im1)));
            end
            
            %errorConstraint = sum(errorConstraint(:)) / (nPx * (obj.numFrames-1))
        end
        
        function showErrorInConstraint(obj)
            
            errorConstraint = obj.getErrorInConstraint;
            
            for j=1:obj.numFrames-1
                figure(123);imagesc(errorConstraint(:,:,j));axis image;colormap(gray);colorbar;
                pause
            end
        end
           
        function solveV(obj)
            obj.mainV.resetImages(obj.u);
            
            if (obj.doWarping)
                obj.mainV.runPyramid;
            else
                obj.mainV.runLevel(1);
            end
            
            obj.v = obj.mainV.getResult;
            
            tmp = obj.v(:,:,:,1);
            obj.v(:,:,:,1) = obj.v(:,:,:,2);
            obj.v(:,:,:,2) = tmp;
            
            if (obj.verbose > 1)
                for j=1:obj.numFrames-1
                    figure(200+j);imagesc(flowToColorV2(cat(3,obj.v(:,:,j,1),obj.v(:,:,j,2))));axis image;
                end %todo: add comparison with ground-truth flow
            end
        end
        
        function [errorU,errorV] = calculateErrors(obj)
            errorU = -1;
            errorV = -1;
            
            %image error
            if (~isscalar(obj.gtU))
                % extract solution and store into matrix u
                for j=1:obj.numFrames
                    obj.u(:,:,j) = obj.mainU.getPrimal(j);
                end
                
                errorU = 0;
                for j=1:obj.numFrames
                    %squared error
                    err = (obj.u(:,:,j) - obj.gtU(:,:,j)).^2;

                    %cut out inner window
                    err = err(20:end-20,20:end-20);
                    
                    figure(124);imagesc(err);colorbar;
                    figure(123);imagesc(err);colorbar;
                    %sum up
                    err = sum(sqrt(err(:))) / numel(err);


                    errorU = errorU + err;
                end
                errorU = errorU / obj.numFrames; %average;
            end

            %flow error
            if (~isscalar(obj.gtV))
                errorV = 0;
                for j=1:obj.numFrames-1
                    field = squeeze(obj.v(:,:,j,:));
                    fieldGT = squeeze(obj.gtV(:,:,j,:));
                    
                    figure(121);imagesc(flowToColorV2(field));
                    figure(122);imagesc(flowToColorV2(fieldGT));

                    err = absoluteError(field,fieldGT);
                    errorV = errorV + err;
                end
                errorV = errorV / (obj.numFrames-1);
            end
        end
        
        function residual = calculateMainResidual(obj)
            residualU = (obj.uOld - obj.u).^2;
            residualU = sqrt(sum(residualU(:))) / numel(residualU);
            
            residualV = (obj.vOld - obj.v).^2;
            residualV = sqrt(sum(residualV(:))) / numel(residualV);
            
            residual = residualU + residualV;
        end
        
        function [energy,ofTerm] = calculateEnergy(obj)
            
            gradientXf = gradientOperator(obj.dims,1,'discretization','forward');
            gradientXb = gradientOperator(obj.dims,1,'discretization','backward');
            gradientY = (gradientXf.matrix + gradientXb.matrix)/2;

            gradientYf = gradientOperator(obj.dims,2,'discretization','forward');
            gradientYb = gradientOperator(obj.dims,2,'discretization','backward');
            gradientX = (gradientYf.matrix + gradientYb.matrix)/2;
            
            energy = 0;
            ofTerm = 0;
            for i=1:obj.numFrames
                u1 = obj.u(:,:,i);

                kuterm = 0.5*(u1 - obj.imageSequence(:,:,i)).^2;
                tvterm = obj.alpha * sqrt((gradientX * u1(:)).^2 + (gradientY * u1(:)).^2);

                energy = energy + sum(kuterm(:)) + sum(tvterm(:));

                if (i < obj.numFrames)
                    u2 = obj.u(:,:,i+1);

                    ux = reshape(gradientX*u2(:),size(u2));
                    uy = reshape(gradientY*u2(:),size(u2));

                    v1 = obj.v(:,:,i,1);
                    v2 = obj.v(:,:,i,2);

                    tvv1term = obj.beta * sqrt((gradientX * v1(:)).^2 + (gradientY * v1(:)).^2);
                    tvv2term = obj.beta * sqrt((gradientX * v2(:)).^2 + (gradientY * v2(:)).^2);

                    energy = energy + sum(tvv1term(:)) + sum(tvv2term(:));

                    ofconstraintTmp = u2 - u1 + ux .* v1 + uy .* v2;
                    ofTerm = ofTerm + sum(abs(ofconstraintTmp(:))) / numel(ofconstraintTmp);
                end
            end
        end
        
        function run(obj)
            
            residual = Inf;
            iteration = 1;
            while (iteration <= obj.numMainIt && residual > obj.tolMain) || iteration <= obj.numMinMainIt
                obj.uOld = obj.u;
                obj.vOld = obj.v;
                
                %% solve u problem
                if (obj.verbose > 0)
                    disp('Solving problem for u');
                end
                obj.solveU;
                %% solve problem v
                if (obj.verbose > 0)
                    disp('Solving problem for v');
                end
                obj.solveV;

                %% update warping operators in flexBox
                if (obj.verbose > 0)
                    disp('Updating operators in flexBox');
                end
                obj.updateFlexBoxU;
                
                %% update residual
                residual = obj.calculateMainResidual;
                iteration = iteration + 1;
                
                
            end
            
            disp('Finished');
        end
        
    end
    
end


