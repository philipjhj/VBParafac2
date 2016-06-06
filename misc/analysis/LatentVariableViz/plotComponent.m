function plotComponent(W,Z,mask,mask_dim,varargin)
% W has dimension V x q, where each column is a spatial map
% Z has dimension q x T x subjects,
% mask is the binary mask as a vector
% mask_dim is the dimensionality of the measured 3D volume
% V is the number of voxels
% T is the number of timepoints
% q is the number of components

% Parse arguments and check if parameter/value pairs are valid 
paramNames = {'inConvention','outConvention','TR','FontSize','LineWidth','save',...
              'Position'};
defaults   = {'Radiological', 'Neurological',2.49,    18    , 2 ,'',...
              [0 0 800 425]};

[inputConvention, outputConvention, TR, f_size, l_size,save_loc,fig_position]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});

%TODO customizable
%Non trivial
viewSetting = [50 60; 130 60; 270 0]; % ViewSetting affect zoomLevel
zoomLevel = [5/3,5/3,1.3]; % Zoom level affects display of s
s = {'','','Back Side'}; %Note, first and second argument plots weirdly

%%Validate custom input
%Radiological convention (patients left is on the right side)
%Neurological convention (patients left is on the left side)
assert(any(strcmpi(inputConvention,{'Radiological','Neurological'})),'Input convention unsupported')
assert(any(strcmpi(outputConvention,{'Radiological','Neurological'})),'Output convention unsupported')

% Show average if multiple subjects are given
Zmean = mean(Z,3);
Zsd = sqrt(var(Z,0,3));
mask = reshape(mask,mask_dim);
T=size(Z,2);
%% Scaling
threshold = [-.9:.1:-.1 , 0.1:0.1:.9];
W = W/max(abs(W(:)));

%% Plotting spatial maps and temporal activation
timepoints = (1:T)*TR;
plotArea = [1 4; 2 5; 3 6];

for j = 1:size(W,2) %Illustrate every component
    
    figure('Position',fig_position)
    %% Plotting S views
    for i = 1:3
        subplot(3,3,plotArea(i,:))
        
        pa = patch(isosurface(mask)); %Plot underlying 3D brain
        set(pa,'edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.2);
        axis off;
        axis equal;
        axis tight;
        view(viewSetting(i,:));
        for t = 1:length(threshold)
            m = mask;
            dispColor = 'Green';
            if threshold(t) > 0
                m(m==1)=W(:,j)>threshold(t);
            else
                dispColor = 'Red';
                m(m==1)=W(:,j)<threshold(t);
            end
            
            %If input convention is different from output convension
            if ~strcmpi(inputConvention,outputConvention)
                m = m(end:-1:1,:,:); %Flip x-axis
            end
            pb = patch(isosurface(m)); %Add colored brain
            set(pb,'edgecolor','none','facecolor',dispColor...
                ,'facealpha',abs(threshold(t)));
        end
        camzoom(zoomLevel(i));
        set(gca,'Fontsize',f_size)
        title(s{i},'FontSize',f_size+2)
        
        % Offset to the left
        pos = get(gca, 'Position');
        xoffset = -0.06;
        if i ==1
            pos(1) = pos(1) + xoffset;
            pos(3) = pos(3) - 0.5*xoffset;
        elseif i==2
            pos(1) = pos(1) +0.5*xoffset;
            pos(3) = pos(3) - 0.5*xoffset;
        elseif i==3
            pos(3) = pos(3) -.5* xoffset;
        end
        set(gca, 'Position', pos)
    end
        
    %% Plot temporal activation Z
    subplot(3,3,7:9)
    hold on
    plot(timepoints,Zmean(j,:),'Color','b','LineWidth',l_size);
    if size(Z,3) > 1
        plot(timepoints,Zmean(j,:)+2*Zsd(j,:),':','Color','r','LineWidth',l_size-.5);
        plot(timepoints,Zmean(j,:)-2*Zsd(j,:),':','Color','r','LineWidth',l_size-.5);
    end
    xlabel('time (s)')
    ylabel('Activation')
    
    set(gca,'fontsize',f_size)
    if nargin <= 7
        %str = sprintf('Latent Component %i: alpha. %3.2f:',j,alpha(j));%,corr_val(j));
        str = sprintf('Latent Component %i',j);
        title(str,'FontSize',f_size+2);
    end
    
    %Move slightly to the left
    pos = get(gca, 'Position');
    xoffset = -0.04;
    pos(1) = pos(1) + xoffset;
    pos(3) = pos(3) -2* xoffset;
    set(gca, 'Position', pos)
    
    %Save file
    if ~isempty(save_loc)
        results_name = sprintf('_d%i',j);
        print(sprintf('%s%s',save_loc,results_name),'-dpng');
        print(sprintf('%s%s',save_loc,results_name),'-depsc')
    end
    
    
end

end