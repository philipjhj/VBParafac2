
dataset ='Cheese';
mydir=strcat('data/',dataset,' data/');
files=dir(mydir);

% myModel=varBayesModelParafac2;
% 1600 combinations in tests below
% i=0;
estimationARD = {'avg','max'};
estimationP = {'vonmises','parafac2svd'};

parfor nfile = 1:numel(files)
    if regexp(files(nfile).name, regexptranslate('wildcard','Int*'))
        for P = 1:length(estimationP)
            for ARD = 1:length(estimationARD)
                for init = 1:10
                    disp(strcat(num2str(nfile),'/70 - ',num2str((ARD-1+2*(P-1))*10+init),'/40'))
                    pause(0.5)
                    %for replicate = 1:5
                    %for I = [10 20]
                    %   for K = [3 6]
                    %      for congruence = [0.4 0.8]
                    for M = 10
                        pathname = 'output_RealData/';
                        interval = files(nfile).name(1:end-4);
                        filename = sprintf('%s_%s_estimationP_%s_estimationARD_%s_M_%d_init_%d',...
                            dataset,interval,estimationP{P},estimationARD{ARD},M,init);%
                        % pathname = 'output/';
                        % filename = sprintf('estimationP_%s_estimationARD_%s_I_%d_K_%d_M_%d_congruence_%.1f_dataset_%d_init_%d',...
                        %     estimationP{P},estimationARD{ARD},I,K,M,congruence,replicate,init);
                        savestr = strcat(pathname,filename,'.mat');
                        % disp(exist(savestr,'file') == 2)
                        % i=i+1
                        if ~(exist(savestr,'file') == 2)
                            % I=20;
                            %                             J=I;
                            % K=6;
                            % M=3;
                            Mesti = M;
                            
                            %                             options.dimensions = [I J K M];
                            %                             options.initMethod = 'kiers';
                            %                             options.congruence = congruence;
                            %                             options.precision = [1e12 1e-3];
                            
                            %                             rng(replicate);
                            %                             data = varBayesModelParafac2.generateDataFromModel(options);
                            
                            
                            filepath = strcat(mydir,files(nfile).name);
                            myData=load(filepath);
                            
                            I_no = regexp(interval,'\d*','Match');
                            
                            
                            
%                             eval(strcat('data = permute(I',I_no{1},',[2 1 3]);'))
                            
                            data = permute(myData.(strcat('I',I_no{1})),[2 1 3]);

                            myModel=varBayesModelParafac2(data,Mesti);
                            
                            % myModel=varBayesModelParafac2(Y,100);
                            
                            myModel.opts.verbose = 0;
                            myModel.opts.debugFlag = 0;
                            myModel.opts.estimationP= estimationP{P};
                            % myModel.opts.estimationP = 'parafac2svd';
                            myModel.opts.estimationARD = estimationARD{ARD};
                            myModel.opts.matrixProductPrSlab = 'gpu';
                            myModel.opts.nActiveComponents = 'hard';
                            myModel.opts.showIter = 1;
                            myModel.opts.rngInput = init;%'shuffle';
                            
                            % data set; rng(3)
                            % seed; 1461191309
                            % iter (avg); 441
                            % iter (max); 44
                            
                            % myModel.opts.maxTime = 1;
                            
                            
                            
                            
                            %myModel.qDist.SNR
                            % clc
                            
                            myModel.qDist.opts.activeParams = {'qA','qF','qC','qP','qAlpha','qSigma'};
                            % myModel.qDist.activeParams_opt = {'qC','qAlpha'};
                            
                            
                            % clc
                            
                            % myModel.data.iter = myModel.data.iter-1;
                            % myModel.restartqDist;
                            myModel.opts.maxTime = 60*60*23.5;
%                             toptic=tic;
                            myModel.computeVarDistribution;
                            
                            %myModel.qDist.SNR
                            %
%                             fprintf('estimationP: %s; estimationARD: %s; init: %d; replicate: %d; I: %d; K: %d; congruence: %.1f; M: %d; ## fit: %f; time: %f\n',...
%                                 estimationP{P},estimationARD{ARD},init,replicate,I,K,congruence,M, myModel.Parafac2Fit,toc(toptic))
                            %
                            
                            %                                 disp(savestr)
                            parsave_named(savestr,{'myModel'},myModel)
                        end
                    end
                end
            end
        end
    end
end