function expmnts = characterizeDocument(f,func)
global fig_save plt
addpath(genpath('chronux_2_11'));
%if function is not specified use default
if ~exist('func')
    func=@parameterCharacterization;
end

if ~isempty(fig_save)
    plt=true;
end

[status,sheets]=xlsfinfo(f);


sheetMask=true(size(sheets));
% for i=1:length(strMask)
%     for j=1:length(sheets)
%         disp(regexpi(sheets{j},strMask{i}));
%         if ~isempty(regexpi(sheets{j},strMask{i}))
%             sheetMask(j)=true;
%         end
%     end
% end
expmnts=containers.Map();
[pathstr,name,ext] = fileparts(f);
if isempty(ext)
    ext='.xlsx';
end
timestamp=regexprep(regexprep(datestr(now),'( |-)','_'),':','.');
nameout=[pathstr filesep() name '_' func2str(func) '_'  timestamp];
if ~isempty(fig_save)
    mkdir(nameout)
end
for i=1:length(sheets)
if sheetMask(i)
    s=sheets{i};
    expmntSet=splitSheet(f,s);
    for j=1:length(expmntSet)
    expmntSet{j}.resps={};
    for k=1:size(expmntSet{j}.fluo,2)
        try
        resp=func(expmntSet{j}.time,expmntSet{j}.fluo(:,k));
        
        if ~isempty(resp) && ~isempty(fieldnames(resp))
            fn=fieldnames(resp);
            str=[fn{1} ' : ' num2str(resp.(fn{1}))];
            for l=2:length(fn);
                fnl=fn{l};
                try
                str=[str ': ' fnl ' : ' num2str(resp.(fnl))];
                catch ex
                    disp(ex)
                end
            end
                disp(str);
        else
            expmntSet{j}.resps{k}=[];
        end
        
            if ~isempty(fig_save)
                saveas(1,[nameout filesep() name '.' regexprep(regexprep(s,'( |-)','_'),':','.') '.block_' int2str(j) '.roi_' int2str(k) '.' fig_save],fig_save)
            end
        
        expmntSet{j}.resps{k}=resp;
        catch ex
            disp(ex)
           expmntSet{j}.resps{k}=[];
        end
    end
    end
     writeSheet([nameout  ext],s,expmntSet);
end
end


%characterize all loaded experiments


end

