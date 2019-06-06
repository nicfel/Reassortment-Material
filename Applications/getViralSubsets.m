% sub sample from datasets
clear

dateorder = [3,1,2];

currdir = cd;

virus = {'h1n1pdm', 'h1n1sea','h3n2', 'infB'};
workingdir = {'H1N1pandemic','H1N1seasonal', 'H3N2', 'InfB'};

% there are only 1037 samples from seasonal h1n1
nrseq = [2000 1037 5000 2000];

for i = 3% : length(virus)
    getSubsampledDataset(virus{i}, workingdir{i}, dateorder, nrseq(i));
    cd(currdir);
end
    
