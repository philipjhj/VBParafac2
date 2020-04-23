I=50;
J=I;
K=100;

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
    
for M = 20 %1:8        
    for rng_seed=1
        rng(rng_seed)
        options.SNR = 10;
        options.noiseType = noiseTypes{1};
        options.dimensions(4) = M;

        data = varBayesModelParafac2.generateDataFromModel(options);

        %save(['example_data_M_',num2str(M),'_seed_',num2str(rng_seed)],'data','options')
        save('sampled_data_many_components_test','data','options')
    end
end