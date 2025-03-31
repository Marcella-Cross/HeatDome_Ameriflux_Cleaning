%% Main function for My_Micromet data processing
% Created by <Marcella> on <Oct. 28, 2024>
% Run each file first and then run this section to do phase 1, 2, 3 and
% output as final ameriflux csv file
% ============================
%Code to paste into R: 
%# Define the list of packages
%packages <- c("tidyverse", "caret", "REddyProc", "dplyr", "lubridate", 
            %  "data.table", "fs", "yaml", "rlist", "zoo", "reshape2", 
             % "stringr", "ranger", "ggplot2")

%install.packages(setdiff(packages, installed.packages()[, "Package"]))

%# Load all packages into the library
%lapply(packages, library, character.only = TRUE)

%% Cleaning:
fr_automated_cleaning([2020, 2021, 2022], 'CA3', [1 2 7 8]);
fr_automated_cleaning([2020, 2021, 2022], 'LP1', [1 2 7 8]);
fr_automated_cleaning([2020, 2021, 2022], 'CF1', [1 2 7 8]);
fr_automated_cleaning([2020, 2021, 2022], 'CF2', [1 2 7 8]);
fr_automated_cleaning([2020, 2021, 2022], 'CF3', [1 2 7 8]);
fr_automated_cleaning([2020, 2021, 2022], 'CF4', [1 2 7 8]);
fr_automated_cleaning([2020, 2021, 2022], 'ME2', [1 2 7 8]);
fr_automated_cleaning([2020, 2021, 2022], 'ME6', [1 2 7 8]);
fr_automated_cleaning([2020, 2021, 2022], 'XWR', [1 2 7 8]);
fr_automated_cleaning([2020, 2021, 2022], 'XAB', [1 2 7 8]);
