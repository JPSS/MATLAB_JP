display('Compiling dataprocessing functions... please wait');
%fit functions
%mex fit_traces.cpp DataProcessing.cpp;
%mex fit_2step.cpp DataProcessing.cpp;
%mex fit_bleaching.cpp DataProcessing.cpp;

%filter functions
%mex nlf.cpp DataProcessing.cpp;
%mex nlf_alex.cpp DataProcessing.cpp;
%mex nlf_fret.cpp DataProcessing.cpp;

%other
%mex studentTtest.cpp DataProcessing.cpp;
mex find_peaks1d.cpp DataProcessing.cpp;



display('Compiling done');



