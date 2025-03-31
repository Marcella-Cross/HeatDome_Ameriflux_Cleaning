function structProject = get_TAB_project_configuration(projectPath)
%This file is generated automatically by create_TAB_ProjectFolders.m
projectName = '';
structProject.projectName   = projectName;
structProject.path          = fullfile(projectPath);
structProject.databasePath  = fullfile(structProject.path,'Database');
structProject.sitesPath     = fullfile(structProject.path,'Sites');
structProject.matlabPath    = fullfile(structProject.path,'Matlab');
