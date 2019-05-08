clear

% define the segments
segments = {'HA', 'M1', 'NA', 'NP', 'NS1', 'PA', 'PB1', 'PB2'};
segments_file = {'HA', 'MP', 'NA', 'NP', 'NS1', 'PA', 'PB1', 'PB2'};


% set the random number generator
rng(1);

% gets a subsampled set from the sequence data in oridat
filename = dir('oridata/*.fasta');

system('rm -r data');
system('mkdir data');

for f = 1 : length(filename)
    % read in the fasta file
    fasta = fastaread(['oridata/' filename(f).name]);c=1;
    for i = 1 : length(fasta)
        tmp = strsplit(fasta(i).Header, '_');     
        tmp2 = strsplit(tmp{2}, '-');
        if length(tmp2)==3
            Data(c) = fasta(i);
            Data(c).Header = [tmp{1} '|' tmp{2}];
            c = c+1;
        end
    end
    fastawrite(['data/' strrep(filename(f).name, 'Inf', 'inf')], Data);
end
