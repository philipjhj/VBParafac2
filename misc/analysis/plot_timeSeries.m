



% Plot component-wise contribution

FtPt=mtimesx(permute(data.Ftrue,[2 1 3]),permute(data.Ptrue,[2 1 3]));

%%
M=size(data.Atrue,2);
K = size(data.Xtrue,3);

Xrecon_m = zeros([size(data.Xtrue) M]);
for m = 1:M
    for k = 1:K
    Xrecon_m(:,:,k,m)=data.Atrue(:,m)*diag(data.Ctrue(k,m))*FtPt(m,:,k);
    end   
end

%%

Z = permute(squeeze(Xrecon_m(:,:,1,1)),[1 2]);


bar3(Z,'stacked')

%%
I=size(Xrecon_m,1);
J=size(Xrecon_m,2);


cols=autumn(9);
for l = 1:5
for m = 1:3
    clf
    for i=1:I
        for k=1:K
            plot3(repmat(k,1,I)*2+0.1*i,1:I,Xrecon_m(i,:,k,m),'LineWidth',2,'Color',cols(m+5,:))
            hold on
        end
    end
    xlabel('K')
    xlim([1.5 10.5])
    xticks(2.5:2:(1+2*K))
    view(-20,45)
     xticklabels(1:K)
ylabel('J')
title(['Component #',num2str(m)])
grid on
zlim([-1.2 1.2])
    pause(3)
end

end

%% Full Data Plot
K=4
clf
    for i=1:I
        for k=1:K
            if k ==1
                m=1;
            else
                m=3;
            end
            plot3(repmat(k,1,I)*2+0.1*i,1:I,data.Xtrue(i,:,k),'LineWidth',2,'Color',cols(m+5,:))
            hold on
        end
    end
    xlabel('Sample #')
    xlim([1.5 10.5])
    xticks(2.5:2:(1+2*K))
    view(-15,45)
     xticklabels(1:K)
ylabel('Retention Time')

title('Full Data')
grid on
zlim([-1.2 1.2])

%% One Sample Plot
K=1
clf
    for i=1:I
        for k=1:K
            plot3(repmat(k,1,I)*2+0.05*i,1:I,data.Xtrue(i,:,k),'LineWidth',2,'Color',cols(1+5,:))
            hold on
        end
    end
    zlabel('Intensity')
    xlim([2 2.6])
    xlabel('Mass Spectrum')
    view(-15,45)
     xticklabels([])
ylabel('Retention Time')
title('Full Data')
grid on
zlim([-1.2 1.2])