% sub sample from datasets
clear

dateorder = [3,1,2];

currdir = cd;

virus = {'h1n1pdm', 'h1n1sea','h3n2', 'infB'};
workingdir = {'H1N1pandemic','H1N1seasonal', 'H3N2', 'InfB'};

nrseq = [500 1000 5000 2000];

for i = 4% : length(virus)
    getSubsampledDataset(virus{i}, workingdir{i}, dateorder, nrseq(i));
    cd(currdir);
end
    
