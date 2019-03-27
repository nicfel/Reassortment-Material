% sub sample from datasets
clear

dateorder = [3,1,2];

currdir = cd;

virus = {'h1n1pdm', 'h1n1sea','h3n2', 'h5n1', 'h7n9', 'h9n2'};
workingdir = {'H1N1pandemic','H1N1seasonal', 'H3N2', 'H5N1', 'H7N9', 'H9N2'};

nrseq = [500 1000 1000 500 500 500]

for i = 5% : length(virus)
    getSubsampledDataset(virus{i}, workingdir{i}, dateorder, nrseq(i));
    cd(currdir);
end
    
