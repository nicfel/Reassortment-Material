clear

currdir = cd;

workingdir = {'H1N1pandemic',...
        'H1N1seasonal',...
        'H3N2',...
        'InfB',...
        'H2N2'};

% delete all xmls in all working dirs
for i = 1 : length(workingdir)
    system(sprintf('rm -r %s/simulation', workingdir{i}));
    system(sprintf('mkdir %s/simulation', workingdir{i}));
end


for i = 1 : length(workingdir)
    convertXmlToSimulation(workingdir{i})
    cd(currdir)
end
