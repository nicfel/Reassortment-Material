function [pairwiseDistances, translate] = getPairwiseDistances(filename, startval)
% get the pairwise distances between taxa from a beast run
f = fopen(filename);
translate = cell(0,0);
pairwiseDistances = cell(0,0);
c = 1;
while ~feof(f)
    line = fgets(f);
    if contains(line, 'Translate')
        line = fgets(f);
        while ~contains(line, ';')
            tmp = strsplit(strrep(strtrim(line), ',', ''));
            translate{str2double(tmp{1})} = tmp{2};
            line = fgets(f);
        end        
    end
    if contains(line, 'tree STATE')
        if startval<=c
            tmp = strsplit(strtrim(line));
            ptree = phytreeread(tmp{end});
            pairdist = pdist(ptree, 'Squareform', true);
            % sort the leafnames correctly
            leafs = get(ptree, 'leafnames');
            sort_indices = zeros(size(leafs));
            for i = 1 : length(sort_indices)
                sort_indices(str2double(leafs{i})) = i;
            end
            pairwiseDistances{c-startval+1} = pairdist(sort_indices,sort_indices);
        end
        c = c+1;    
    end
end
    
end