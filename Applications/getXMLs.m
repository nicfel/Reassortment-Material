clear

currdir = cd;

virus = {'h1n1pdm', 'h1n1sea','h3n2', 'h5n1', 'h7n9', 'h9n2', 'h3n2old', 'h3n2recent', 'infB', 'infC', 'infD'};
workingdir = {'H1N1pandemic','H1N1seasonal', 'H3N2', 'H5N1', 'H7N9', 'H9N2', 'H3N2', 'H3N2', 'InfB', 'InfC', 'InfD'};

from = [2000, 1990, 1980, 2000, 2010, 1960, 1980, 1995, 1980, 1980, 1980];
to =   [2019, 2010, 2010, 2019, 2019, 2019, 1995, 2010, 2019, 2019, 2019];


temperature = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01, 0.025];
temperature = temperature*2;

% delete all xmls in all working dirs
for i = 1 : length(workingdir)
    system(sprintf('rm -r %s/xmls', workingdir{i}));
    system(sprintf('mkdir %s/xmls', workingdir{i}));
end


for i = 1 : length(virus)
    getXMLallsegments(virus{i}, workingdir{i}, 100, from(i), to(i), temperature(i));
    cd(currdir);
end

