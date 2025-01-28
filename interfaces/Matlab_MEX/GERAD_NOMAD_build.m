%% GERAD NOMAD Build for Matlab

% This file will help you compile NOMAD for use with MATLAB. 

% To recompile you will need to get / do the following:

% 1) Get NOMAD
% NOMAD is available from https://www.gerad.ca/nomad
% Complete the download form then download the latest version. Define the
% $NOMAD_HOME environment variable.

% 2) Start Matlab and go into $NOMAD_HOME/examples/interfaces/Matlab_Mex
% The NOMAD MEX Interface is a simple MEX interface written to use NOMAD.

% 3) Compile the MEX File by executing this file in Matlab
%
%
% The code below will build the NOMAD MEX file and set the Matlab path. 

clear nomad

% Current directory
cdir = cd;

% Check and set nomad_home and create variables for path
clear nomad_root nomad_src nomad_src_sgtelib nomad_build_lib;

% Default values
% nameLibNomad = '';
%updateLDFLAGS= '';
%install_name_tool='';

% Get a default directory for Nomad root
cd ..
cd ..
nomad_root = cd;

if ( ~ exist(nomad_root,'dir') )
    error('Cannot access Nomad root directory.');
end
    
% Attempt to access build dir (release, debug)
% This is not robust!
%
nomad_build_lib = [nomad_root filesep 'build' filesep 'release' filesep 'lib' filesep];
if ( ~ exist(nomad_build_lib,'dir') )
    nomad_build_lib = [nomad_root filesep 'build' filesep 'debug' filesep 'lib' filesep];
    if ( ~ exist(nomad_root,'dir') )
        error('Cannot access Nomad build directory (release and debug). Make sure to build nomad and check the build dir path in this script.');
    end
end
    
nomad_src=[nomad_root filesep 'src' filesep];
if ( ~ exist(nomad_src,'dir') )
    error('The default Nomad source directory does not exist. Please make sure that it exists.');
end
nomad_src_sgtelib=[nomad_root filesep 'ext' filesep 'sgtelib' filesep 'src' filesep ];
if ( ~ exist(nomad_src_sgtelib,'dir') )
    error('The default Sgtelib source directory does not exist. Please make sure that it exists.');
end

% Return to base dir
cd(cdir);

if ( strcmp(computer,'PCWIN64') == 1 || strcmp(computer,'PCWIN32') == 1 )
   
    % Nomad library name
    nameLibNomad = 'nomadAlgos.lib nomadUtils.lib nomadEval.lib sgtelib.lib';
    tmpMexLibName = 'nomad_tmp';
    finalMexLibName = 'nomad';
    
    %Compile & MOVE (Windows) ---> recompile Nomad and sgtelib
    post = [' -I.  -I' nomad_src ' -I' nomad_src_sgtelib ' -L' nomad_build_lib ' -lut -output ' nomad_build_lib filesep tmpMexLibName '.' mexext];
    pre = ['mex -v ' countEvalFlag ' -largeArrayDims nomadmex.cpp ' nameLibNomad ];
        
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LINUX AND OSX  ---> use dynamic libraries
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Default update LDFLAGS (linux only)
    updateLDFLAGS= '';
    % Post compilation tool for path to library (osx only)
    install_name_tool='';
    
    % nameLibNomad = 'libnomadAlgos.so libnomadUtils.so libnomadEval.so libsgtelib.so';
    switch(computer)
        case 'GLNX86'
            updateLDFLAGS = 'LDFLAGS=''$LDFLAGS -Wl,-rpath,''''$ORIGIN/../lib/'''' -Wl,-rpath-link,''''../lib/'''' '' ';
        case 'GLNXA64'
            updateLDFLAGS = 'LDFLAGS=''$LDFLAGS -Wl,-rpath,''''$ORIGIN/../lib/'''' -Wl,-rpath-link,''''../lib/'''' '' ';
        case 'MACI64'
            % change library names
            % nameLibNomad = 'libnomadAlgos.dylib libnomadUtils.dylib libnomadEval.dylib libsgtelib.dylib';
            install_name_tool=['install_name_tool -change ' nameLibNomad ' @loader_path/../lib/' nameLibNomad ' ' nomad_build_lib filesep 'nomad.' mexext];
    end
   
    %Compile & Move (Default) --> use shared object library
    post = [' -I.  -I' nomad_src ' -I' nomad_src_sgtelib ' -lut -lnomadAlgos -lnomadUtils -lnomadEval -L' nomad_build_lib ' -output ' nomad_build_lib filesep 'nomad.' mexext ];
    pre =[ 'mex -g -v ' countEvalFlag ' -largeArrayDims nomadmex.cpp ' updateLDFLAGS ];
    
end
    

fprintf('\n------------------------------------------------\n');
fprintf('NOMAD MEX FILE BUILD --- GERAD VERSION \n\n');

%CD to Source Directory
cd 'Source';

try


    eval([pre post])
    
    if ( strcmp(computer,'MACI64') == 1 )
        system(install_name_tool);
    end
    
    if ( strcmp(computer,'PCWIN64') == 1 || strcmp(computer,'PCWIN32') == 1 )
        movefile([nomad_bin filesep tmpMexLibName '.' mexext],[nomad_lib filesep finalMexLibName '.' mexext]);
    end
    
    cd(cdir);
    fprintf('Compilation done!\n');
    fprintf('\n----------------------------------------------------------------------------------------------\n');
    fprintf(' To be able to use the nomad functions, you may need to modify the Matlab path \n');
    qstring = 'To be able to use the nomad functions, you may need to modify the Matlab path. Do you want to update the Matlab path?';
    if ( usejava('desktop'))
    	choice = questdlg(qstring,'Set path','Yes','No','Yes');
    else
        choice = 'Yes';
    end

    if ( strcmp(choice,'Yes') )
        addpath([ cdir filesep 'Functions']);
        addpath(nomad_build_lib);
        fprintf('  ---> The Matlab path has been modified but not saved.\n');
    end
    clear nomad_root nomad_build_lib nomad_src nomad_src_sgtelib cdir post pre updateLDFLAGS qstring choice install_name_tool nameLibNomad;
catch ME
    cd(cdir);
	clear nomad_root nomad_build_lib nomad_src nomad_src_sgtelib cdir post pre updateLDFLAGS qstring choice install_name_tool nameLibNomad;
    error('Error Compiling NOMAD!\n%s',ME.message);
end
