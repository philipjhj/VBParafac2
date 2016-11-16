
datasets = {'Apple','Aroma','Cheese','Wine'};


dimArray = zeros(1,3);
    k=0;

for i = 1:4
    
    mydir=strcat('/media/data/DataAndResults/Thesis/data/dataBro/Models and data/',datasets{i},' data/');
    files=dir(mydir);
    %
    % all_logp = cell(3,10,2);
    % count = zeros(3,10);
    % count_all = zeros(3,10,20000);
    
    
    for file = files'
        if regexp(file.name, regexptranslate('wildcard','Int*'))
%             disp(file.name)
            load(strcat(file.folder,'/',file.name))
            
            I_no = regexp(file.name,'\d*','Match');

            eval(strcat('I_dim=size(I',I_no{1},');'))
            
            
            dimArray = dimArray + I_dim;
            k=k+1;
            %         filepath = strcat(mydir,file.name);
            %         m = matfile(filepath);
            %         syntheticdatasetscript = m.syntheticdatasetscript;
            %         fliptype=m.fliptype;
            %         load(filepath)
            %         if isempty(all_logp{fliptype+1,syntheticdatasetscript,1})
            %             all_logp{fliptype+1,syntheticdatasetscript,1} = zeros(1,20000);
            %             all_logp{fliptype+1,syntheticdatasetscript,2} = zeros(1,20000);
            %          end
        end
    end
    
end

% 
% ans =
% 
%    57.1419  380.3613   53.2323
