I=50;
J=I;
K=10;
M=4;

options.dimensions = [I J K M];
options.initMethod = 'kiers';
% options.initMethod = 'generative';
options.congruence = 0.4;
% 1e4 1e-3 i ARD tests
options.precision = [1e2 1e-6];

% [1e4 1e-8] creates problems for qC



noiseTypes={'homoscedastic','heteroscedastic'};
for SNR = 2:2:10
    for noiseTypeIDX = 1:2
        for rngSeed=2:10
            rng(rngSeed)
            options.SNR = SNR;
            options.noiseType = noiseTypes{noiseTypeIDX};

            data = varBayesModelParafac2.generateDataFromModel(options);

            filepath=sprintf(['/media/data/DataAndResults/VBParafac2paper/data/' ...
                'simulatedTestSNR/simulatedData_SNR_%s_noiseType_%s_%d'],...
                num2str(options.SNR),options.noiseType,rngSeed);
            
            save(filepath,'data','options')

        end
    end
end

