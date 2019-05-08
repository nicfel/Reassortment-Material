clear

currdir = cd;

virus = {'h3n2dum'};
workingdir = {'H3N2'};

from = [1980];
to =   [ 2010];


temperature = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01, 0.025, 0.01];
temperature = temperature*2;

% delete all xmls in all working dirs
for i = 1 : length(workingdir)
    system(sprintf('rm -r %s/dummyxmls', workingdir{i}));
    system(sprintf('mkdir %s/dummyxmls', workingdir{i}));
end


for i = 1 : length(virus)
    getXMLallsegmentsDummy(virus{i}, workingdir{i}, 100, from(i), to(i), temperature(i));
    cd(currdir);
end

