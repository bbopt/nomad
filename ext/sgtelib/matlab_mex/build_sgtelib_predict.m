%% Sgtelib_predict Build for Matlab

% The code below will build the sgtelib_predict MEX file and set the Matlab path.

% Check and set sgtelib_home and create variables for path
clear sgtelib_home sgtelib_src sgtelib_bin ;

sgtelib_predict = 'sgtelib_predict';
sgtelib_predict_source = [ sgtelib_predict '.cpp' ];

sgtelib_home = getenv('SGTELIB_HOME');

% Current directory
cdir = cd;

if ( length(sgtelib_home) > 1)
    warning('The sgtelib_home variable for Matlab is set to %s. The default can be replaced by using the command setenv(''SGTELIB_HOME'',ARG1)! before running the build_sgtelib_predict.m script.',sgtelib_home);
    if ( ~isempty( find(isspace(sgtelib_home),1) ) )
        error('The compilation of sgtelib_predict for Matlab must be performed in the matlab_mex directory. The path should not contain empty space.');
    end
else
    cd ..
    sgtelib_home = cd; 
    if ( ~ exist(sgtelib_home,'dir') )
            error('The default sgtelib directory does not exist. Please set the SGTELIB_HOME variable using the command setenv(''SGTELIB_HOME'',ARG1)! before running the build_sgtelib_predict.m script.');
    end
    cd (cdir);
    
end
  
sgtelib_src=[sgtelib_home filesep 'src' filesep];


% Default values
nameLibSgtelib = '';
updateLDFLAGS= '';
install_name_tool='';

if ( strcmp(computer,'PCWIN64') == 1 || strcmp(computer,'PCWIN32') == 1 )
    sgtelib_bin=[sgtelib_home filesep 'bin' filesep];
    nameLibSgtelib = 'sgtelib.lib';
    post = [' -I.  -I' sgtelib_src ' -L' sgtelib_bin ' -lut -output ' sgtelib_bin sgtelib_predict '.' mexext];
    pre = ['mex -v -largeArrayDims ' sgtelib_predict_source ' ' nameLibSgtelib ];
else
    sgtelib_bin=[sgtelib_home filesep 'lib' filesep];
    nameLibSgtelib = 'libsgtelib.so';

    switch(computer)
        case 'GLNX86'
            updateLDFLAGS = 'LDFLAGS=''$LDFLAGS -Wl,-rpath,''''$ORIGIN/../lib/'''' -Wl,-rpath-link,''''../lib/'''' '' ';
        case 'GLNXA64'
            updateLDFLAGS = 'LDFLAGS=''$LDFLAGS -Wl,-rpath,''''$ORIGIN/../lib/'''' -Wl,-rpath-link,''''../lib/'''' '' ';
        case 'MACI64'
            install_name_tool=['install_name_tool -change ' nameLibSgtelib ' @loader_path/../lib/' nameLibSgtelib ' ' sgtelib_bin filesep 'nomad.' mexext];
    end
    post = [' -I. -I' sgtelib_src ' -lut -lsgtelib -L' sgtelib_bin ' -output ' sgtelib_bin sgtelib_predict '.' mexext ];
    pre =[ 'mex -v -largeArrayDims ' sgtelib_predict_source updateLDFLAGS ];
end

if ( ~ exist([sgtelib_bin nameLibSgtelib],'file') )
        error('The SGTELIB library file %s is not available. Please perform sgtelib project compilation before proceeding.',nameLibSgtelib);      
end

fprintf('\n-------------------------------\n');
fprintf('SGTELIB_PREDIC MEX FILE BUILD \n\n');


try

    eval([pre post])
    
    if ( strcmp(computer,'MACI64') == 1 )
        system(install_name_tool);
    end

    fprintf('Compilation done!\n');
    fprintf('\n----------------------------------------------------------------------------------------------\n');
    fprintf(' To be able to use the sgtelib_predict function, you may need to modify the Matlab path \n');
    qstring = 'To be able to use the sgtelib_predict function, you may need to modify the Matlab path. Do you want to update the Matlab path?';
    choice = questdlg(qstring,'Set path','Yes','No','Yes');
    if ( strcmp(choice,'Yes') )
        addpath(sgtelib_bin);
        fprintf('  ---> The path has been modified but not saved.\n');
    end
    cd (cdir)
    clear cdir sgtelib_home sgtelib_bin sgtelib_src post pre choice ;
    
catch ME
    cd (cdir)
	clear sgtelib_home sgtelib_bin sgtelib_src cdir post pre choice ;
    error('Error Compiling sgtelib_predict !\n%s',ME.message);
end
