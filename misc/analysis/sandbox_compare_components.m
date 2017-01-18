


pMethod = 'parafac2svd';
mEsti = '7';
ARD = 'off';
datasetRNG = '1';
initRNG = '1';


myDir = '/media/data/DataAndResults/Thesis/output/results/results_ARD_tests/';
myFile = ['ARD_tests_dim_20_20_10_4__pMethod_',pMethod,...
    '_mEsti_',mEsti,'_ARD_',ARD,'_datasetRNG_',datasetRNG,'_initRNG_',initRNG,'.mat'];


load(strcat(myDir,myFile))

myModel.compute_reconstruction;

disp(repmat(['-'],1,100))
disp(['pMethod: ',pMethod,'mEsti: ',mEsti,'ARD: ',ARD,'datasetRNG: ',datasetRNG,'initRNG: ',initRNG])
disp(repmat(['-'],1,100))    
fprintf('\n\n')



for m = 1:myModel.data.M
disp(repmat(['-'],1,100))
disp(['Component ',num2str(m),' compared to the true'])
disp(repmat(['-'],1,100))    

Xrecon_m = myModel.data.Xrecon_m(:,:,:,m);

for mtrue = 1:myModel.data.Mtrue
%     disp(mtrue)
    Xtrue_m = myModel.data.Xtrue_m(:,:,:,mtrue);
%     disp(congruenceScore(Xrecon_m(:),Xtrue_m(:)))
    fprintf('%20.10f',congruenceScore(Xrecon_m(:),Xtrue_m(:)));
    fprintf('/');
    fprintf('%.10f',corr(Xrecon_m(:),Xtrue_m(:)));
    
%     fprintf('\n')
end
fprintf('\n')
end
fprintf('\n\n')


%%
% fill3(X,Y,Z,C)

X = myModel.data.X;
I = myModel.data.I;
J = myModel.data.J;
K = myModel.data.K;


% pValues = [0.025,0.25,0.35,0.65,0.75,0.975];
pValues = [0.025,0.975];
confintervals = quantile(myModel.data.X,pValues,3);
nConfIntervals = size(confintervals,3)-1;
alphaValues = linspace(0.35,0.3,(nConfIntervals+1)/2);
alphaValues = [alphaValues fliplr(alphaValues(1:end-1))];

nColors=8;
myColors = lines(nColors);

% figure
clf
hold on

for k = 1:K
    for i = 1:I
        %         plot3(1:J,repmat(i,1,J),X(i,:,k))
        
        for cinf = 1:nConfIntervals
            
%             if sum(sum(confintervals(i,:,cinf:cinf+1))) > 100*J*K
                xFill = [1:J fliplr(1:J)];
                zFill = [confintervals(i,:,cinf) fliplr(confintervals(i,:,cinf+1))];
                %         keyboard
                fill3(xFill,repmat(i,1,2*J),zFill,myColors(mod(i-1,nColors)+1,:),...
                    'facealpha',alphaValues(cinf),'LineWidth',2,'EdgeColor',myColors(mod(i-1,nColors)+1,:))
%                 keyboard
%             end
        end
    end
end
xlabel('retention time')
ylabel('Mass Spectrum')
zlabel('Elution')

% view(3)
view(215.3,37.2)
grid on
set(gca,'Color',[0.9,0.9,0.9])


%%
clf
for m = 1:myModel.data.M

    subplot(1,2,m)
    
    for k = 1:myModel.data.K
        for i = 1:myModel.data.I
            plot3(1:myModel.data.J,repmat(k,1,myModel.data.J),myModel.data.Xrecon_m(i,:,k,m))
            hold on
        end
    end
    xlabel('retention time')
    ylabel('k')
    zlabel('Elution')
end