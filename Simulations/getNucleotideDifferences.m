function [nucdiff] = getNucleotideDifferences(filename)
% get the nucleotide differences between samples
f = fopen(filename);
sequence = cell(0,0);
while ~feof(f)
    line = fgets(f);
    if contains(line, 'matrix')
        line = fgets(f);
        while ~contains(line, 'end;')                
            tmp = strsplit(strtrim(line));
            sequence{end+1,1} = strrep(tmp{end}, ';','');
            line = fgets(f);
        end
    end
end
% calculate the number of differences between 2 taxa
nucdiff = zeros(length(sequence),length(sequence));
for a = 1 : length(sequence)-1
    for b = a+1 : length(sequence)
        nucdiff(a,b) = sum(sequence{a}~=sequence{b});        
    end
end

fclose(f);

end