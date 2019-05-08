clear

currdir = cd;

virus = {'h1n1pdm', 'h1n1sea','h3n2', 'h3n2ancient', 'h3n2old', 'h3n2recent', 'h3n2new', 'infB','h2n2'};
workingdir = {'H1N1pandemic','H1N1seasonal', 'H3N2', 'H3N2', 'H3N2', 'H3N2', 'H3N2', 'InfB','H2N2'};

from = [2000, 1990, 1980, 1980, 1990, 2000, 2010, 1980, 1920];
to =   [2019, 2010, 2020, 1990, 2000, 2010, 2020, 2019, 1970];


temperature = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01];
temperature = temperature*2;

% delete all xmls in all working dirs
% for i = 1 : length(workingdir)
%     system(sprintf('rm -r %s/xmls', workingdir{i}));
%     system(sprintf('mkdir %s/xmls', workingdir{i}));
% end




for i = 8 : length(virus)-1
    if i==3
        getXMLallsegmentsDummy(virus{i}, workingdir{i}, 200, from(i), to(i), temperature(i), 1234567);
        cd(currdir)
        getJointCoalallsegments(virus{i}, workingdir{i}, 200, from(i), to(i), temperature(i), 1234567);
    elseif i==7
        getXMLallsegmentsDummy(virus{i}, workingdir{i}, 200, from(i), to(i), temperature(i), 12345678);
        cd(currdir)
        getJointCoalallsegments(virus{i}, workingdir{i}, 200, from(i), to(i), temperature(i), 12345678);
    else
        getXMLallsegmentsDummy(virus{i}, workingdir{i}, 200, from(i), to(i), temperature(i),randi(1000000));
    end
    cd(currdir)
end

% i=9;
% getXMLallsegments(virus{i}, workingdir{i}, 200, from(i), to(i), temperature(i), 1234567);
% cd(currdir)

