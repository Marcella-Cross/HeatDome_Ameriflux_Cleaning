%% US-Me6
% Created by <Marcella> on <Oct. 28, 2024>
% 
% ============================
projectPath = '/Users/Marcella Cross/Desktop/HeatDome_Ameriflux_Cleaning/MicroMet/Projects/My_Micromet/';
structProject=set_TAB_project(projectPath);
siteID = 'ME6';

% Create database from raw data
%% Flux data from EddyPro output files
%
% Input file name
fileName = fullfile(structProject.sitesPath, siteID, 'Flux','AMF_US-Me6_BASE_HH_16-5.csv');

assign_in = 'caller';
varName = [];
dateColumnNum = [2];
timeInputFormat = {'uuuuMMddHHmm'};
colToKeep = [3 Inf];
structType = 1;
inputFileType = 'delimitedtext';
modifyVarNames = 0;
VariableNamesLine = 3;
rowsToRead = [];
isTimeDuration =[];
[~,~,~,outStruct] = fr_read_generic_data_file(fileName,assign_in,...
                 varName, dateColumnNum,timeInputFormat ,colToKeep,structType,inputFileType,modifyVarNames,VariableNamesLine,rowsToRead,isTimeDuration);
%%
% set database path 
databasePath = fullfile(db_pth_root,'yyyy',siteID,'Flux'); 

% Convert outStruct into database 
missingPointValue = NaN; 
timeUnit= '30MIN'; 
structType = 1; 
db_struct2database(outStruct,databasePath,0,[],timeUnit,missingPointValue,structType,1); 

%% 
structProject = set_TAB_project('/Users/mcross/Desktop/FluxProjectPipeline/MicroMet/Projects/My_Micromet/')
stationIDs = [28051];
yearsIn = 2021:2023;
monthsIn  = 1:12;
 
for cntStations = 1:length(stationIDs)
    sID = stationIDs(cntStations);
    pathECCC = fullfile('yyyy','ECCC',num2str(sID));
 
    try
        db_ECCC_climate_station(yearsIn,monthsIn,sID,...
                                fullfile(structProject.databasePath,pathECCC),60);
    catch
        fprintf('Error processing station: %d (year: %d, month: %d)\n',sID,yearsIn,monthsIn(end));
    end
end

%% Plot Data
%Multiple years
yearsIn = 2021:2023;                                    % loading multiple years in one go
pth = biomet_path('yyyy','Me6','MET');                   % find data base path for multiple years, BB2 site
tv = read_bor(fullfile(pth,'clean_tv'),8,[],yearsIn);   % load the time vector (Matlab's datenum format)
tv_dt = datetime(tv,'convertfrom','datenum');           % convert to Matlab's datetime object (use for all new stuff)
x = read_bor(fullfile(pth,'MET_HMP_T_2m_Avg'),[],[],yearsIn); % load MET_CNR4_Net_Avg trace from BB2/MET folder
plot(tv_dt,x)                                           % plot data
grid on; zoom on;
%% ini file creation - this will create an ini file in the working directory
% Move file to Calculation_Procedures\TraceAnalysis_ini\CA3
structSetup.startYear = 2020;
structSetup.startMonth = 1;
structSetup.startDay = 1;
structSetup.endYear = 2022;
structSetup.endMonth = 12;
structSetup.endDay = 31;
structSetup.Site_name = 'ME6';
structSetup.SiteID = 'ME6';
structSetup.allMeasurementTypes = {'MET','Flux'};
structSetup.Difference_GMT_to_local_time = 8;  % local+Difference_GMT_to_local_time -> GMT time
structSetup.outputPath = []; % keep it in the local directory
createFirstStageIni(structSetup)

%% Once first stage ini file has been created:
fr_automated_cleaning([2020, 2021, 2022], 'ME6', 1);

%% Once second stage ini file has been created:
fr_automated_cleaning([2020, 2021, 2022], 'ME6', 2);
%%
fr_automated_cleaning([2020, 2021, 2022], 'ME6', [1 2 7 8]);