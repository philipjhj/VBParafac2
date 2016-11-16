classdef analysisVBParafac2 < handle
    
    properties
        
        root
        root_figures
        root_csv
        
        testResults
        countArray
        ELBOArray
        testOpts
        
    end
    
    
    methods
%         
        function obj = analysisVBParafac2(root)
            
            if isempty(root)
                obj.root = 'output/';
            else
                obj.root = root;
            end
            
            obj.root_figures = strcat(obj.root,'figures/');
            obj.root_csv = strcat(obj.root,'csv/');
            
        end
        %
        %
        
        % # Generic functions
        %
%                 function loadTestResults(obj
%                     % Generic function to load results
%         
%         
%         
%                 end
        %
        %         function saveTestResults
        %             % Save
        %
        %         end
        %
        
        function obj = computeTestResults(obj,testDir)
            
            
            % filename = 'ARD_tests_dim_20_20_10_4_pMethod_parafac2svd_mEsti_2_ARD_avg_datasetRNG_1_initRNG_1*.mat';
            %testDir='/media/data/DataAndResults/Thesis/output/results/results_ARD_tests/';
            
            files=dir(testDir);
            files=cat(1,{files.name})';
            
            %pMethods={'parafac2svd','vonmises'};
            %ARDMethods={'max','avg','off'};
            
            opt = struct;
            
            % Find options from
            % <data_info>__<opt1>_<val1>_..._<optN>_<valN>.mat
            % TODO: save <data_info>
            for file_i = 1:numel(files)
                if regexp(files{file_i},'.*\.mat')
                    
                    optionsString = regexp(files{file_i},'(?!.*__)(?!_).*','match');
                    optionsStruct = regexp(optionsString,'(.+?)(?:_|\.)','tokens');
                    
                    all_options=cat(1,optionsStruct{:}{:});
                    
                    
                    for i = 1:2:numel(all_options)
                        if ~isfield(opt,all_options{i})
                            opt.(all_options{i}) = {all_options{i+1}};
                        else
                            opt.(all_options{i}){end+1} = all_options{i+1};
                        end
                    end
                end
            end
            
            % Determine number of unique values for each option
            opt_field = fields(opt);
            nOpts = numel(opt_field);
            dims = zeros(1,nOpts);
            
            for opt_i = 1:nOpts
                optUniq.(opt_field{opt_i}) = unique(opt.(opt_field{opt_i}));
                dims(opt_i) = numel(optUniq.(opt_field{opt_i}));
            end
            
            
            fitarray = zeros([2 dims]);
            countarray = zeros(dims);
            ELBOarray = zeros(dims);
            
            
            % Sort numbers correctly
            if isfield(optUniq,'initRNG')
                [~,idx]=sort(cellfun(@str2num,optUniq.initRNG));
                optUniq.initRNG=optUniq.initRNG(idx);
            end
            
            if isfield(optUniq,'datasetRNG')
                [~,idx]=sort(cellfun(@str2num,optUniq.datasetRNG));
                optUniq.datasetRNG=optUniq.datasetRNG(idx);
            end
            
            
            % Find Best Solution over repeated tests
            
            % Compute fit
            file_j = 1; 
            for file_i = 1:numel(files)
                if regexp(files{file_i},'.*\.mat')
                    load(strcat(testDir,files{file_i}));
                    
                    [fit, fit_true]=myModel.Parafac2Fit;
                    %                         disp(fit)
                    %                         disp(fit_true)
                    %
                    fit_idx = zeros(1,nOpts);
                    for opt_i = 1:nOpts
                        fit_idx(opt_i) = find(ismember(optUniq.(opt_field{opt_i}),opt.(opt_field{opt_i}){file_j}));
                    end
                    
                    
                    fitarray(1,fit_idx(1),fit_idx(2),fit_idx(3),fit_idx(4),fit_idx(5)) = fit;
                    fitarray(2,fit_idx(1),fit_idx(2),fit_idx(3),fit_idx(4),fit_idx(5)) = fit_true;
                    countarray(fit_idx(1),fit_idx(2),fit_idx(3),fit_idx(4),fit_idx(5)) = countarray(fit_idx(1),fit_idx(2),fit_idx(3),fit_idx(4),fit_idx(5))+1;
                    ELBOarray(fit_idx(1),fit_idx(2),fit_idx(3),fit_idx(4),fit_idx(5)) = myModel.ELBO_chain(myModel.data.iter-1);
                    
                    file_j=file_j+1; % Only increment this if file is approved
                end
            end
            
            obj.testResults = fitarray;
            obj.countArray = countarray;
            obj.ELBOArray = ELBOarray;
            obj.testOpts = optUniq;
%             obj.testValues = 
            
        end
        
        
        function [X,A,C,F,P,dims] = loadParafac2(obj,Parafac2obj)
            % Function to load model fit based on class type
            if isa(Parafac2obj,'varBayesModelParafac2')
                X = Parafac2obj.data.X;
                A = Parafac2obj.qDist.qA.mean;
                C = Parafac2obj.qDist.qC.mean;
                F = Parafac2obj.qDist.qF.mean;
                P = Parafac2obj.qDist.qP.mean;
            elseif isa(Parafac2obj,'normalParafac2')
                X = Parafac2obj.X;
                A = Parafac2obj.A;
                C = Parafac2obj.C;
                F = Parafac2obj.F;
                P = Parafac2obj.P;
            end
            
            dims = [size(X) size(A,2)];
            
        end
        
        % # Plot functions
        
        function saveFig(obj,plotC,outdir,plotname,format)
            
            if nargin < 5
                format = '.png';
            end
            
            exportPath = strcat(obj.root_figures,outdir,'/');
            exportName = strcat(exportPath,plotname);
            
            figPath = strcat(obj.root_figures,outdir,'/plotdata/');
            if ~exist(exportPath,'dir')
                mkdir(exportPath)
                mkdir(figPath)
            end
            
            plotC.export(strcat(exportName,format))
            savefig(strcat(figPath,plotname))
        end
        
        function formatPlot(obj,plotC,plotType)
            
             switch plotType
                 case 'ElutionProfile'
                     plotC.BoxDim = [8, 8];
             end
            
            
        end
        
        
        
        function plotARDtest(obj)
            
            val_ELBO = max(obj.ELBOArray,[],5);
            
            % [i1,i2,i3,i4,i5]=ind2sub(size(obj.ELBOArray),find(ismember(obj.ELBOArray(:),nonzeros(val_ELBO(:)))));
            idx_max_ELBO = find(ismember(obj.ELBOArray(:),nonzeros(val_ELBO(:))));
            
            for s = size(obj.testResults,1)
            
            fit1 = obj.testResults(s,:,:,:,:,:);
            fit1 = fit1(idx_max_ELBO);
            
            fit1_final = zeros(size(obj.ELBOArray));
            fit1_final(idx_max_ELBO) = fit1;
            fit1_final = sum(fit1_final,5);
            fit1_final = sum(fit1_final,4)./sum(fit1_final~=0,4);
            

            ARD_choice = [1 2 3];
            for pMethod = 1:2
                figure
                for ARD = ARD_choice
                    
                    dataMethod = fit1_final(pMethod,:,ARD);
                    dataMethod = reshape(dataMethod,size(dataMethod,1),numel(dataMethod)/size(dataMethod,1));
                    
                    plot(2:7,dataMethod')
                    hold on
                    axis tight
                    
                        
                    %         end
                end
                hold off
                
                plt = Plot();
                plt.BoxDim = [8, 8];
%                 plt.LineStyle = {'-.','-.'};
%                 plt.Colors = {'r','b','g'};
                plt.YLim = [70 100];
                plt.Legend = obj.testOpts.ARD(ARD_choice);
                plt.LegendLoc = 'northeastoutside';
                plt.Title = obj.testOpts.pMethod(pMethod);
                
                set(gca,'position',[0.075 0.075 0.875 0.875],'units','normalized');
                set(gca,'position',[0.075 0.075 0.875 0.875],'units','normalized');
                
               
%                 obj.saveFig(plt,'ARDTestPlot',['ARDTestPlot_',obj.testOpts.pMethod{pMethod}])
                
            end
            
            
            
            end
        end
        
        
        function plotReconElutionProfiles(obj,Parafac2obj,plotname,saveFlag)
            % Plot of Elution Profiles per reconstructed component
            
            %             set(0,'DefaultFigureWindowStyle','docked')
            
            if nargin < 4
                saveFlag = 0;
            end
            
            
            [X,A,C,F,P,dims] = loadParafac2(obj,Parafac2obj);
            
            I=dims(1);
            J=dims(2);
            K=dims(3);
            M=dims(4);
            
            Xmin = min(min(X,[],3),[],1);
            Xmax = max(max(X,[],3),[],1);
            Xidx = 1:numel(Xmin);
            
            % Plotting
            figure
            
            fill([Xidx fliplr(Xidx)],[Xmin, fliplr(Xmax)],'r')
            
            % Formatting
            yl = ylim;
            yMax = yl(2);
            axis tight
            
            plt = Plot();
            obj.formatPlot(plt,'ElutionProfile')
            plt.YLim = yl;
            
            ylim(yl);
            
            % Saving
            if saveFlag
                obj.saveFig(plt,plotname,'dataProfile');
            end
            
            for m = 1:M
               
                figure
                reconM = zeros(I,J,K);
                for k = 1:K
                    reconM(:,:,k) = A(:,m)*C(k,m)*F(:,m)'*P(:,:,k)';
                end
                reconMmin = min(min(reconM,[],3),[],1);
                reconMmax = max(max(reconM,[],3),[],1);
                idx = 1:numel(reconMmin);
                fill([idx fliplr(idx)],[reconMmin, fliplr(reconMmax)],'r')
%                 plot(reconMmax,'.b')
                yl = ylim;
                axis tight
                
                plt = Plot();
                obj.formatPlot(plt,'ElutionProfile')
                plt.YLim = [yl(1) yMax];
                
                title(['Component ',num2str(m)])
                if saveFlag
                    obj.saveFig(plt,plotname,['Component ',num2str(m)]);
                end
            end
        
            %             set(0,'DefaultFigureWindowStyle','normal')
            
            if saveFlag
                % Make montage of components
                plotpath = [' ',obj.root_figures,plotname];
                system(strcat('montage',plotpath,'/Com*.png -geometry 250x250 ',plotpath,'/plotdata/all_comps.png'));
            end
            
        end
            
    end
    
    
end