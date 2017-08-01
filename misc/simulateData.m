I=150;
J=I;
K=20;
M=4;

options.dimensions = [I J K M];
options.initMethod = 'kiers';
% options.initMethod = 'generative';
options.congruence = 0.4;
% 1e4 1e-3 i ARD tests
options.precision = [1e2 1e-6];
options.SNR = 5;
options.noiseType = 'homoscedastic';
% [1e4 1e-8] creates problems for qC

filepath=sprintf(['/media/data/DataAndResults/VBParafac2paper/data/' ...
    'simulatedTest/simulatedSmallTestData_SNR_%s_noiseType_%s_'],...
    num2str(options.SNR),options.noiseType);

for rngSeed=1:4
rng(rngSeed)
data = varBayesModelParafac2.generateDataFromModel(options);

save(strcat(filepath,num2str(rngSeed)),'data','options')
end

