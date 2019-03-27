function [] = getSubsampledDataset(virus, workingdir, dateorder, nrSampled)

cd(workingdir)


% define the segments
segments = {'HA', 'M1', 'NA', 'NP', 'NS1', 'PA', 'PB1', 'PB2'};
segments_file = {'HA', 'MP', 'NA', 'NP', 'NS1', 'PA', 'PB1', 'PB2'};


% set the random number generator
rng(1);
% gets a subsampled set from the sequence data in oridat
filename = dir('oridata/*.fasta');
% read in the fasta file
fasta = fastaread(['oridata/' filename(1).name]);
% get the names and sampling times of all sequences
name = cell(length(fasta),1);
segmentname = cell(length(fasta),1);
year = zeros(length(fasta),1);
for i = 1 : length(fasta)
    tmp = strsplit(fasta(i).Header, '|');     
    if length(tmp)>=5
        name{i} = [tmp{1} '|' tmp{4} '|' tmp{5}];
        segmentname{i} = tmp{3};
        time = strsplit(tmp{4},'/');
        if length(time)>2
            year(i) = str2double(time{end});
        end
    else
        name{i} = 'dont_use';
    end
end


%% get all the unique names
unique_names = unique(name);
indices = zeros(length(unique_names), length(segments));
sampling_year = zeros(length(unique_names), 1);
% get the indices of all 8 segments of each of these sequences
for i = 1 : length(unique_names)
    same_indices = find(ismember(name, unique_names{i}));
    if length(same_indices) == 8        
        for j = 1 : length(segments)
            seg = segmentname{same_indices(j)};
            sampling_year(i) = year(same_indices(j));
            if length(seg)>0
                segment_indice = find(ismember(segments, segmentname{same_indices(j)}));
                indices(i, segment_indice) = same_indices(j);
            end
        end
    end
    if sum(indices(i,:)>0)~=8
        sampling_year(i) = 0;
    end
end    

% get the weights for the subsampling
unique_years = unique(sampling_year(sampling_year>0));

nr_samples = 0;    
% sample nrSampled sequences at random
use_sample = zeros(nrSampled,1);
while nr_samples < nrSampled
    use_year = randsample(length(unique_years),1);
    in_year = find(sampling_year == unique_years(use_year));
    if length(in_year)>0
        use_ind = randsample(length(in_year),1);
        use_sample(nr_samples+1,1) = in_year(use_ind);
        nr_samples = nr_samples+1;
        sampling_year(in_year(use_ind)) = 0;
    end
end

system('rm -r data');
system('mkdir data');
%%
% sort the samples to use
use_sample = sort(use_sample);
% build fasta files for each of the segments
for i = 1 : length(segments)
    clear Data
    for j = 1 : length(use_sample)
        Data(j) = fasta(indices(use_sample(j),i));
        tmp = strsplit(fasta(indices(use_sample(j),i)).Header, '|');
        tmp2 = strsplit(tmp{4}, '/');
        Data(j).Header = strrep([tmp{1} '|' tmp2{dateorder(1)} '-' tmp2{dateorder(2)} '-' tmp2{dateorder(3)}], '''','');
        if str2double(tmp2{1})>12
            error('order of dates is wrong')
        end
    end
    fastawrite(['data/' virus '_' segments_file{i} '.fas'], Data )
end

disp('align sequences...')
for i = 1 : length(segments)
    in = ['data/' virus '_' segments_file{i} '.fas'];
    out = ['data/' virus '_' segments_file{i} '.afas'];
    system(sprintf('./../../Software/muscle3.8.31_i86darwin64 -in %s -out %s -maxiters 1 -diags', in, out));
end
%%
disp('make FastTree trees')
system('rm -r treetime')
keep_leafs = cell(length(segments),1);
for i = 1 : length(segments)
    system('rm tree_file.trees')
    in = ['data/' virus '_' segments_file{i} '.afas'];
    disp('build fastTree');
    system(sprintf('./../../Software/FastTree -gtr -nt < %s > tree_file.nwk', in));
    
    % build a file with the sampling times
    disp('build sampling time map')
    fasta = fastaread(in);
    f = fopen('dates.csv','w');
    fprintf(f, 'name,date\n');
    for j = 1 : length(fasta)
        fprintf(f, '%s,',fasta(j).Header);
        tmp = strsplit(fasta(j).Header, '|');
        time = tmp{end};
        tmp = strsplit(time,'-');
        
        if length(tmp) == 2
            tmp{3} = '15';
            time = [time '-15'];
        end

        
        deztime = (datenum(time,'yyyy-mm-dd')- datenum(tmp{1},'yyyy'))...
                    /(datenum(num2str(str2double(tmp{1})+1),'yyyy')-datenum(tmp{1},'yyyy'))...
                    +str2double(tmp{1});        
        fprintf(f, '%.5f\n', deztime);
    end
    fclose(f);
    disp('build timetree');
    system(sprintf('rm -r timetree/%s',segments_file{i}))
    system(sprintf('/usr/local/bin/treetime --clock-filter 0 --reroot ML --aln %s --tree tree_file.nwk --dates dates.csv --outdir timetree/%s',in,segments_file{i}));
    
end


%%
for i = 1 : length(segments)    
    disp('read in time tree')
    f = fopen(['timetree/' segments_file{i} '/divergence_tree.nexus']);
    while ~feof(f)
        line = fgets(f);
        if length(line)>1000
            tree = line;
        end
    end
    tree = strrep(strtrim(tree), 'Tree tree1=', '');
    tree=regexprep(tree, '\[(.*?)\]', '');
    tree=regexprep(tree, 'NODE_(\d*)', '');


    ptree = phytreeread(tree);
    
    
    
    all_leaves = get(ptree, 'leafnames');
    % get the branch lengths
    branches = get(ptree, 'Distances');
    prunebranch = find(branches>0.1);    
    if length(prunebranch)>0
        largest_set = cell(length(prunebranch),1);
        for prunenr = 1 : length(prunebranch)
            righttree = prune(ptree, prunebranch(prunenr));
            rightleafs = get(righttree, 'leafnames');
            
            if length(all_leaves)-length(rightleafs)>1
                lefttree = subtree(ptree, prunebranch(prunenr));
                % get the leafnames
                leftleafs = get(lefttree, 'leafnames');
                if length(rightleafs) > length(leftleafs)
                    largest_set{prunenr} = rightleafs;
                else
                    largest_set{prunenr} = leftleafs;
                end
            else
                largest_set{prunenr} = rightleafs;
            end
        end
        
        % get the intersection
        intersectset = largest_set{1};
        if length(prunebranch)>1
            for prunenr = 2 : length(prunebranch)
                intersectset = intersect(intersectset, largest_set{prunenr});
            end
        end
        
        
        keep_leafs{i} = intersectset;
    else
        keep_leafs{i} = all_leaves;
    end    
    disp('done')
end
system('rm tree_file.nwk')
system('rm dates.csv')


%% reduce the alignment to only keep leaves that are present in all alignments
% get the intersection
useleaves = keep_leafs{1};
for i = 2 : length(keep_leafs)
    useleaves = intersect(useleaves, keep_leafs{i});
end


disp('reduce taxa...')
system('rm data/*.fasta');
for i = 1 : length(segments)
    infasta = fastaread(['data/' virus '_' segments_file{i} '.afas']);
    for j = length(infasta):-1:1
        ind = find(ismember(infasta(j).Header, useleaves));
        if isempty(ind)
            infasta(j) = [];
        end
    end
    % delete everything before the first ATG
    start = strfind(infasta(1).Sequence, 'ATG');
    for j = 1:length(infasta)        
        if start>1
            infasta(j).Sequence(1:(start-1)) = [];
        end
    end
    % write to file
    fastawrite(['data/' virus '_' segments_file{i} '.fasta'], infasta);
    
end
disp('done')
system('rm data/*.afas');
system('rm data/*.fas');



