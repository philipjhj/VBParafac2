classdef analysisVBParafac2 < handle
    
    properties
        
        % Savepaths for output
        root
        root_figures
        root_csv
        
        
        
        % Meta of loaded data
        results_path
        test_title
        
        test_dir % strcat of results_path and test_title
        nFiles
        
        % Values from tests
        fitArray
        countArray
        ELBOArray
        
        % Best ELBO info
        max_ELBO_lin_idx
        max_ELBO_sub_idx
        
        max_ELBO_lin_idx_failed
        max_ELBO_sub_idx_failed
        
        max_ELBO_filenames
        max_ELBO_filenames_failed
        
        
        % Information about loaded tests
        testOpts % Struct of options and possible values
        testOpts_names % Option names
        testOpts_count % Number of options
        testOpts_dims % dimension of tests (not including data sets or initializations)
        
        models_count % Models with highest ELBO
        
        data_names
        data_count
        
        % Results
        sortOrder = [2 1 3];
        
        AllResults
        nFoundComponents_low_threshold
        nFoundComponents_high_threshold
        
        % Plot settings
        fontsize = 15;
        
        % Data
        X
        
        % Estimated
        A
        C
        F
        P
        dims
        
        Xrecon % Reconstruction of X
        
        % True if simulated, ground truth (best estimated) if real data
        Atrue
        Ctrue
        Ftrue
        Ptrue
        
        
    end
    
    
    methods
        %
        function obj = analysisVBParafac2(root)
            
            if nargin < 1 || isempty(root)
                obj.root = 'output/';
            else
                obj.root = root;
            end
            
            obj.root_figures = strcat(obj.root,'figures/');
            obj.root_csv = strcat(obj.root,'csv/');
            
        end
        
        
        function obj = computeTableELBOAll(obj,resultsPath,testTitle)
            
            obj.results_path = resultsPath;
            obj.test_title = testTitle;
            
            obj.test_dir = strcat(obj.results_path,obj.test_title,'/');
            
            
            ResultsFile_name = strcat(obj.results_path,obj.test_title,'.mat');
            
            files=dir(obj.test_dir);
            files=cat(1,{files.name})';
            nFiles_temp = numel(files);
            
            if (exist(ResultsFile_name,'file')==2)
                resultsData = load(ResultsFile_name,'obj');
                myObj = resultsData.obj;
                
                
                % ---------------------------------------------------------
                % REDO collection if new number of files
                % ---------------------------------------------------------
                if myObj.nFiles ~= nFiles_temp
                    collect_results = 1;
                else
                    collect_results = 0;
                end
            else
                collect_results = 1;
            end
            
            if collect_results
                
                obj.nFiles = nFiles_temp;
                
                opt = struct;
                
                % ---------------------------------------------------------
                % Find options from
                % <data_name>__<opt1>_<val1>_..._<optN>_<valN>.mat
                % ---------------------------------------------------------
                
                data_names = cell(1,obj.nFiles);
                for file_i = 1:obj.nFiles
                    if regexp(files{file_i},'.*\.mat')
                        
                        fullname = regexp(files{file_i},'(.*)__(.*)','tokens');
                        
                        optionsString = fullname{1}{2};
                        optionsStruct = regexp(optionsString,'(.+?)(?:_|\.)','tokens');
                        
                        data_names{file_i} = fullname{1}{1};
                        
                        all_options=cat(1,optionsStruct{:});
                        
                        
                        for i = 1:2:numel(all_options)
                            if ~isfield(opt,all_options{i})
                                opt.(all_options{i}) = {all_options{i+1}};
                            else
                                opt.(all_options{i}){end+1} = all_options{i+1};
                            end
                        end
                    end
                end
                
                data_names_all = data_names(~cellfun('isempty',data_names));
                
                data_names = unique(data_names_all);
                
                % Determine number of unique values for each option
                opt_field = fields(opt);
                nOpts = numel(opt_field);
                dims_opts = zeros(1,nOpts);
                
                for opt_i = 1:nOpts
                    optUniq.(opt_field{opt_i}) = unique(opt.(opt_field{opt_i}));
                    dims_opts(opt_i) = numel(optUniq.(opt_field{opt_i}));
                end
                
                % Sort numbers correctly
                if isfield(optUniq,'initRNG')
                    [~,idx]=sort(cellfun(@str2num,optUniq.initRNG));
                    optUniq.initRNG=optUniq.initRNG(idx);
                end
                
                if isfield(optUniq,'datasetRNG')
                    [~,idx]=sort(cellfun(@str2num,optUniq.datasetRNG));
                    optUniq.datasetRNG=optUniq.datasetRNG(idx);
                end
                
                
                
                obj.data_count = numel(data_names);
                obj.data_names = data_names; %#ok<*PROPLC>
                
                
                
                fitarray = zeros([obj.data_count 2 dims_opts]);
                countarray = zeros([obj.data_count dims_opts]);
                ELBOarray = -realmax*ones([obj.data_count dims_opts]);
                
                % Create table over all tests and repetitions
                
                file_j = 1;
                for file_i = 1:obj.nFiles
                    
                    if regexp(files{file_i},'.*\.mat')
                        load(strcat(obj.test_dir,files{file_i}));
                        
                        if obj.data_count ~= 1
                            data_idx = find(ismember(obj.data_names,data_names_all{file_j}));
                        else
                            data_idx = 1;
                        end
                        
                        % Find option setting for current file
                        Opt_idx = zeros(1,nOpts);
                        for opt_i = 1:nOpts
                            Opt_idx(opt_i) = find(ismember(optUniq.(opt_field{opt_i}),opt.(opt_field{opt_i}){file_j}));
                        end
                        
                        [fit, fit_true]=myModel.Parafac2Fit;
                        
                        % LinIdx(size,sub_idx) computes linear index from vector input
                        
                        fitarray(LinIdx([obj.data_count 2 dims_opts],[data_idx 1 Opt_idx])) = fit;
                        fitarray(LinIdx([obj.data_count 2 dims_opts],[data_idx 2 Opt_idx])) = fit_true;
                        
                        idx_linear = LinIdx([obj.data_count dims_opts],[data_idx Opt_idx]);
                        
                        countarray(idx_linear) = countarray(idx_linear)+1;
                        ELBOarray(idx_linear) = myModel.ELBO_chain(myModel.data.iter-1);
                        
                        file_j=file_j+1; % Only increment this if file is approved
                    end
                    
                end
                
                obj.fitArray = squeeze(fitarray);
                obj.countArray = squeeze(countarray);
                obj.ELBOArray = squeeze(ELBOarray);
                obj.testOpts = optUniq;
                obj.testOpts_names = fields(obj.testOpts);
                obj.testOpts_count = numel(obj.testOpts_names);
                
                save(ResultsFile_name,'obj');
                
                
            else
                
                data = load(ResultsFile_name,'obj');
                obj = data.obj;
            end
        end
        
        
        function obj = find_max_ELBO_lin_idx(obj)
            % Does not work correctly if only one permutation of parameters
            % is given
            
            if ismember('initRNG',obj.testOpts_names)
                val_max_ELBO = max(obj.ELBOArray,[],ndims(obj.ELBOArray));
                obj.max_ELBO_lin_idx = find(ismember(obj.ELBOArray(:),val_max_ELBO(val_max_ELBO(:)~=-realmax)));
                obj.models_count = numel(obj.max_ELBO_lin_idx);
                
                obj.max_ELBO_lin_idx_failed = find(ismember(val_max_ELBO(:),val_max_ELBO(val_max_ELBO(:)==-realmax)));
                
            else
                
            end
        end
        
        function obj = find_max_ELBO_sub_idx(obj)
            
            if isempty(obj.max_ELBO_lin_idx)
                obj.find_max_ELBO_lin_idx;
            end
            
            nDimELBOarray = ndims(obj.ELBOArray);
            
            myIdx = cell(1,nDimELBOarray);
            myIdx_failed = cell(1,nDimELBOarray-1);
            
            sizeELBOarray = size(obj.ELBOArray);
            
            [myIdx{:}]=ind2sub(sizeELBOarray,obj.max_ELBO_lin_idx);
            [myIdx_failed{:}]=ind2sub(sizeELBOarray(1:end-1),obj.max_ELBO_lin_idx_failed);
            
            obj.max_ELBO_sub_idx = sortrows(cat(2,myIdx{:}),[obj.sortOrder 4:nDimELBOarray]);
            obj.max_ELBO_sub_idx_failed = sortrows(cat(2,myIdx_failed{:}),...
                [obj.sortOrder 4:nDimELBOarray-1]);
            
        end
        
        % Function to get max_ELBO filesname (SEE FUNCTION BELOW)
        function obj = find_max_ELBO_filenames(obj)
            
            filenames_temp = cell(obj.models_count,1);
            
            
            for i = 1:obj.models_count
                
                if obj.data_count > 1
                    data_idx = obj.max_ELBO_sub_idx(i,1);
                    myString = obj.data_names{data_idx};
                    k=2;
                else
                   myString = obj.data_names{1};
                   k= 1;
                end
                
                myString = strcat(myString,'_');
                
                for j = 1:obj.testOpts_count
                    if numel(obj.testOpts.(obj.testOpts_names{j})) > 1
                        myString = strcat(myString,'_',obj.testOpts_names{j},...
                            '_',obj.testOpts.(obj.testOpts_names{j})(...
                            obj.max_ELBO_sub_idx(i,k)));
                        k = k+1;
                    else
                        myString = strcat(myString,'_',obj.testOpts_names{j},...
                            '_',obj.testOpts.(obj.testOpts_names{j})(1));
                    end
                end
                
%                 disp(myString)
                
                 filenames_temp{i} = strcat(myString,'.mat');
                
            end
            obj.max_ELBO_filenames = cat(1,filenames_temp{:});
        end
        
        function dispBestRuns(obj)
            for i = 1:obj.models_count
                disp(obj.max_ELBO_filenames{i})
%                 disp(obj.AllResults{i}.nFoundComponents)
            end
        end
        
        
        function obj = computeTableResultsFull(obj)
            
            
            
            ResultsFile_name = strcat(obj.results_path,obj.test_title,'_TableResultsFull.mat');
            
            if (exist(ResultsFile_name,'file')==2)
            
                
                resultsData = load(ResultsFile_name,'obj');
                obj = resultsData.obj;
                
            else
                
                obj.find_max_ELBO_lin_idx;
                obj.find_max_ELBO_sub_idx;
                
                obj.find_max_ELBO_filenames;
                
                
                obj.AllResults = cell(obj.models_count,1);
                
                
                for i = 1:obj.models_count
                    load(strcat(obj.test_dir,obj.max_ELBO_filenames{i}))
                    
                    %
                    myModel.compute_reconstruction;
                    
                    obj.AllResults{i}.ELBO = myModel.ELBO_chain(end);
                    obj.AllResults{i}.fit = myModel.Parafac2Fit;
                    
                    
                    
                    %
                    %                 disp(repmat(['-'],1,100))
                    %                 disp(obj.max_ELBO_filenames{i})
                    %                 disp(repmat(['-'],1,100))
                    %                 fprintf('\n\n')
                    %
                    %
                    
                    
                    Xrecon_m_full = myModel.data.Xrecon_m;
                    Xtrue_m_full = myModel.data.Xtrue_m;
                    
                    mEsti = myModel.data.M;
                    mTrue = myModel.data.Mtrue;
                    
                    
                    obj.AllResults{i}.congruenceCompare = zeros(mEsti,mTrue);
                    obj.AllResults{i}.corrCompare = zeros(mEsti,mTrue);
                    %
                    for m = 1:mEsti
                        %                     disp(repmat(['-'],1,100))
                        %                     disp(['Component ',num2str(m),' compared to the true'])
                        %                     disp(repmat(['-'],1,100))
                        
                        Xrecon_m = Xrecon_m_full(:,:,:,m);
                        
                        for mtrue = 1:mTrue
                            %     disp(mtrue)
                            Xtrue_m = Xtrue_m_full(:,:,:,mtrue);
                            %     disp(congruenceScore(Xrecon_m(:),Xtrue_m(:)))
                            congruence_m = congruenceScore(Xrecon_m(:),Xtrue_m(:));
                            corr_m = corr(Xrecon_m(:),Xtrue_m(:));
                            
                            obj.AllResults{i}.congruenceCompare(m,mtrue) = congruence_m;
                            obj.AllResults{i}.corrCompare(m,mtrue) = corr_m;
                            
                            
                            %                         fprintf('%20.10f',congruence_m);
                            %                         fprintf('/');
                            %                         fprintf('%.10f',corr_m);
                            %
                            %     fprintf('\n')
                        end
                        %                     fprintf('\n')
                    end
                    %                 fprintf('\n\n')
                    %
                    
                    
                end
                save(ResultsFile_name,'obj');
            
            end
            
        end
        
        function obj = computeTableFoundComponents(obj)
        
            obj.testOpts_dims = max(obj.max_ELBO_sub_idx(:,1:3));
            unique_sub_idx = sortrows(unique(obj.max_ELBO_sub_idx(:,1:3),'rows'),obj.sortOrder);
            
            nFoundComponents_temp1 = -1*ones(prod(obj.testOpts_dims),10);
            nFoundComponents_temp2 = -1*ones(prod(obj.testOpts_dims),10);
            
            for i = 1:obj.models_count
                
%                 rowIdx = LinIdx(obj.testOpts_dims,obj.max_ELBO_sub_idx(i,1:3));
                rowIdx = ismember(unique_sub_idx,obj.max_ELBO_sub_idx(i,1:3),'rows');
                colIdx = obj.max_ELBO_sub_idx(i,4);
                
                mTrue=size(obj.AllResults{i}.congruenceCompare,2);
                
                nFoundComponents_temp1(rowIdx,colIdx) = sum(sum(obj.AllResults{i}.congruenceCompare>=0.85));
                nFoundComponents_temp2(rowIdx,colIdx) = sum(sum(obj.AllResults{i}.congruenceCompare>=0.95));
                
            end
%             


            obj.nFoundComponents_low_threshold = nFoundComponents_temp1;
            obj.nFoundComponents_high_threshold = nFoundComponents_temp2;
            
        end
        
        
        
        function obj = loadParafac2(obj,Parafac2obj)
            % Function to load model fit based on class type
            if isa(Parafac2obj,'varBayesModelParafac2')
                obj.X = Parafac2obj.data.X;
                obj.A = Parafac2obj.qDist.qA.mean;
                obj.C = Parafac2obj.qDist.qC.mean;
                obj.F = Parafac2obj.qDist.qF.mean;
                obj.P = Parafac2obj.qDist.qP.mean;
            elseif isa(Parafac2obj,'normalParafac2')
                obj.X = Parafac2obj.X;
                obj.A = Parafac2obj.A;
                obj.C = Parafac2obj.C;
                obj.F = Parafac2obj.F;
                obj.P = Parafac2obj.P;
            end
            
            % I x J x K x M
            obj.dims = [size(obj.X) size(A,2)];
            
        end
        
        
        
        
        
        
        
        
        
        % -----------------------------------------------------------------
        % ######### PLOT FUNCTIONS #################
        % -----------------------------------------------------------------
        
        
        function saveFig(obj,outdir,plotname,format)
            
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
            
%             plotC.export(strcat(exportName,format)) % for pubPlot package
            
            saveas(gcf,strcat(exportName,format))
            savefig(strcat(figPath,plotname))
        end
        
        function formatPlot(obj,plotType,labels)
            
            switch plotType
                case 'ElutionProfile'
%                     plotC.BoxDim = [8, 8];
                case 'TableOfResults'
                    
                    % Labels
                    if isfield(labels,'xtitle')
                        xlabel(labels.xtitle,'interpreter','latex')
                    end
                    if isfield(labels,'ytitle')
                        ylabel(labels.ytitle,'interpreter','latex')
                    end
                    
                    if isfield(labels,'ynames')
                        set(gca,'TickLabelInterpreter','latex')
                        set(gca,'YTickLabel',labels.ynames)
                    end
                    if isfield(labels,'xnames')
                        set(gca,'TickLabelInterpreter','latex')
                        set(gca,'XTickLabel',labels.xnames)
                    end
                    
                    if isfield(labels,'plottitle')
                        title(labels.plottitle,'interpreter','latex')
                    end
                    
                    
                    set(gcf,'Position',[0 0 700 1240]);
                    set(gca,'tickdir','out');
                    set(gca,'fontsize',obj.fontsize);
%                     colormap(gca,(copper));
                    tightfig;
            end
            
            
        end
        
        
        function plotTableResults(obj,type)
            
            
            if nargin<2
                type = 'congruence';
            end
            
            if strcmp(type,'congruence')
                obj.computeTableFoundComponents;
                table{1} = obj.nFoundComponents_low_threshold;
                table{2} = obj.nFoundComponents_high_threshold;
                savenames = {'_low_thres','_high_thres'};
            end
            
            for t = 1:numel(table)
                %Plot - Tests x Data
                imagesc(table{t})
                colorbar
                axis image
                yticks(1:size(table{t},1))
                xticks(1:size(table{t},2))
                
                % Make labels for tests
                testConfig = sortrows(...
                    unique(obj.max_ELBO_sub_idx(:,1:3),'rows'),obj.sortOrder);
                ynames = cell(1,size(testConfig,1));
                
                for i = 1:size(testConfig,1)
                    
                    myString = '';
                    for j = 1:3
                        myString = sprintf('%s %s',myString,...
                            obj.testOpts.(obj.testOpts_names{j}){...
                            testConfig(i,j)});
                    end
                    ynames{i} = myString;
                end
                
                % ynames = unique(cat(1,ynames{:}));
                
                labels.plottitle = ' ';
                labels.xtitle = 'Data sets';
                labels.ynames = ynames;
                
                obj.formatPlot('TableOfResults',labels)
                
                plotnamedir = strcat(obj.test_title,'_result_table_',type);
                plotnamefile = strcat(plotnamedir,type,savenames{t});
                obj.saveFig(plotnamedir,plotnamefile)
            end
        end
        
        function plotARDtest(obj)
            
            val_ELBO = max(obj.ELBOArray,[],5);
            
            % [i1,i2,i3,i4,i5]=ind2sub(size(obj.ELBOArray),find(ismember(obj.ELBOArray(:),nonzeros(val_ELBO(:)))));
            idx_max_ELBO = find(ismember(obj.ELBOArray(:),nonzeros(val_ELBO(:))));
            
            for s = size(obj.fitArray,1)
                
                fit1 = obj.fitArray(s,:,:,:,:,:);
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
            
            
            obj.loadParafac2(Parafac2obj);
            
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