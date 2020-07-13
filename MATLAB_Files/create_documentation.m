close all;
clear all;
clc;

% The next line creates all the documentation in HTML format in the folder
% doc. Also, it creates a graph with all the flow.
% The script that creates this documentation in the folder m2html and its
% webpage is https://www.artefact.tk/software/matlab/m2html.

try
    rmdir('doc', 's');
    disp('Removing old version...');
    for i = 1:5
        pause(1);
        disp(['Please wait ',num2str(i),'/5']);
    end
end

m2html('mfiles','./', 'htmldir','doc', 'recursive','on', 'global','on',...
    'download','on','graph','on');