% plot elution profiles
% close all

[A,F,C,P]=parafac2(permute(I2,[2 1 3]),10);
%%
% X = myModel.data.X;
% A = myModel.qDist.qA.mean;
% C = myModel.qDist.qC.mean;
% F = myModel.qDist.qF.mean;
% P = myModel.qDist.qP.mean;

% X = permute(I2,[2 1 3]);

[I,J,K]=size(X);
M = size(A,2);

f = figure;%('visible','off');

cols = hsv(M);

legNames = cell(1,M);

plotRows = 4;

% for k = 1:K
subplot(plotRows,3,1:2)
% hold on
% plot(X(:,:,k)')
% end
Xmin = min(min(X,[],3),[],1);
Xmax = max(max(X,[],3),[],1);
Xidx = 1:numel(Xmin);
fill([Xidx fliplr(Xidx)],[Xmin, fliplr(Xmax)],'r','FaceAlpha',0.3)
hold on
plot(Xmax,'.b')


% pause(1)
% 
yl = ylim;
yMax = yl(2);
axis tight
ylim(yl);

for m = 1:M
    subplot(plotRows,3,2+m)
%     figure(1+m)
    hold on
    legNames{m} = ['Comp ',num2str(m)];
    reconM = zeros(I,J,K);
    for k = 1:K
%         reconM(:,:,k) = A(:,m)*C(k,m)*F(:,m)'*P(:,:,k)';
        reconM(:,:,k) = A(:,m)*C(k,m)*F(:,m)'*P{k}';
        %plot((reconM)','color',cols(m,:));
    end
    reconMmin = min(min(reconM,[],3),[],1);
        reconMmax = max(max(reconM,[],3),[],1);
        idx = 1:numel(reconMmin);
        fill([idx fliplr(idx)],[reconMmin, fliplr(reconMmax)],'r')
        plot(reconMmax,'.b')
    yl = ylim;
    axis tight
    ylim([yl(1) yMax]);
    title(['Component ',num2str(m)])
%     pause(1)
end

% print -djpeg test.jpg
% close(f)

% legend(plts(1,:),legNames)
%%


f = figure;%('visible','off');

alpha = 0.4;
cols = hsv(M);

legNames = cell(1,M);

for k = 1:K
hold on
plot(X(:,:,k)')
end
Xmin = min(min(X,[],3),[],1);
Xmax = max(max(X,[],3),[],1);
Xidx = 1:numel(Xmin);
% fill([Xidx fliplr(Xidx)],[Xmin, fliplr(Xmax)],'r','FaceAlpha',alpha)
hold on
plot(Xmax,'.b')

yl = ylim;
yMax = yl(2);
axis tight
ylim(yl);

for m = 1:M
    m = 11-m;
    legNames{m} = ['Comp ',num2str(m)];
    reconM = zeros(I,J,K);
    for k = 1:K
        reconM(:,:,k) = A(:,m)*C(k,m)*F(:,m)'*P(:,:,k)';
%         reconM(:,:,k) = A(:,m)*C(k,m)*F(:,m)'*P{k}';
        %plot((reconM)','color',cols(m,:));
    end
    reconMmin = min(min(reconM,[],3),[],1);
        reconMmax = max(max(reconM,[],3),[],1);
        idx = 1:numel(reconMmin);
        fill([idx fliplr(idx)],[reconMmin, fliplr(reconMmax)],cols(m,:),'FaceAlpha',alpha)

end

% print -djpeg test.jpg
% close(f)

% legend(plts(1,:),legNames)
