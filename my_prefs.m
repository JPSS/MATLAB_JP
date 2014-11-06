% default directories
data_dir = [ getenv('HOME') filesep 'Documents' ]; % where data is typically stored
matlab_dir = userpath;
matlab_dir = matlab_dir(1:end-1); % home matlab folder

% variables
scrsz = get(0,'ScreenSize'); % screen size

% default style
set(0, 'defaulttextinterpreter', 'none')
%set(0,'defaultAxesFontName', 'Times-Roman')
%set(0,'defaultTextFontName', 'Times-Roman')
set(0,'defaultlinelinewidth',2)
set(0,'DefaultLineMarkerSize',10)


%testtest