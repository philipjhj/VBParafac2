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

noiseTypes={'homoscedastic',
            'heteroscedastic'
            };
        
rng(1)
options.SNR = 2;
options.noiseType = noiseTypes{1};

data = varBayesModelParafac2.generateDataFromModel(options);

save('example_data','data','options')

