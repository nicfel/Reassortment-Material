clear

% define the segments
segments = {'HA', 'M1', 'NA', 'NP', 'NS1', 'PA', 'PB1', 'PB2'};
segments_file = {'HA', 'MP', 'NA', 'NP', 'NS1', 'PA', 'PB1', 'PB2'};

f = fopen('dates.txt');
used_names =cell(0,0);
while ~feof(f)
    line = strsplit(fgets(f), '_');
    if length(line) ==3
        used_names{end+1,1} = line{1};
    elseif length(line)==5
        used_names{end+1,1} = [line{2} ' ' line{3}];
    else
        used_names{end+1,1} = line{2};
    end
        
end

% set the random number generator
rng(1);

% gets a subsampled set from the sequence data in oridat
filename = dir('oridata/*.fasta');

system('rm -r data');
system('mkdir data');

sequences = cell(0,0);
seqs = cell(0,0);
for f = 1 : length(filename)
    % read in the fasta file
    fasta = fastaread(['oridata/' filename(f).name]);c=1;
    for i = 1 : length(fasta)
        tmp = strsplit(fasta(i).Header, '|');  
        % look if the sequences are in the used name cell
        ind = find(ismember(used_names,tmp{1}));
        if ~isempty(find(ismember(used_names,tmp{1})))
            seqs{end+1,1} = fasta(i).Header;
        elseif ~isempty(strfind(tmp{1}, 'Dummy'))
            seqs{end+1,1} = fasta(i).Header;
        end
    end
end

unique_sequences = unique(seqs);
for i = length(unique_sequences) : -1 : 1
    if ~isempty(strfind(unique_sequences{i}, '2005'))
        unique_sequences(i) = [];
    end
end

%%

for f = 1 : length(filename)
    % read in the fasta file
    fasta = fastaread(['oridata/' filename(f).name]);c=1;
    
    clear Data
    
    seq_name = cell(0,0);
    for i = 1 : length(fasta)
        tmp = strsplit(fasta(i).Header, '|');     
        seq_name{i,1} = fasta(i).Header;
    end
    
    for i = 1 : length(unique_sequences)
        indices = find(ismember(seq_name, unique_sequences{i}));
        if length(indices)==0
            Data(i) = fasta(1);
            Data(i).Header = unique_sequences{i};
            Data(i).Sequence = repmat('N', 1, length(Data(i).Sequence));
            
            tmp = strsplit(Data(i).Header, '|');     
            tmp2 = strsplit(tmp{2}, ' ');
            if ~isempty(strfind(tmp{2}, '-'))
                Data(i).Header = [tmp{1} '|' tmp{2}];
            else
                Data(i).Header = [tmp{1} '|' tmp2{1} '-XX-XX'];
            end
        else       
            takeval = indices(1);
            
            tmp = strsplit(fasta(takeval).Header, '|');     
            tmp2 = strsplit(tmp{2}, ' ');
            if ~isempty(strfind(tmp{2}, '-'))
                Data(i) = fasta(takeval);
                Data(i).Header = [tmp{1} '|' tmp{2}];
            else
                Data(i) = fasta(takeval);
                Data(i).Header = [tmp{1} '|' tmp2{1} '-XX-XX'];
            end
        end
        Data(i).Header = strrep(Data(i).Header, ' ', '');
    end
    fastawrite(['data/' strrep(filename(f).name, 'Inf', 'inf')], Data);
end
