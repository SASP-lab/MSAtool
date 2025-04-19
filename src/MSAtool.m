function [results,model] = MSAtool(varargin) %Full version
%function MSAtool(varargin)

    % University of Minnesota Twin Cities
    % Department of Civil, Environmental, and Geo- Engineering
    % Author: Prof. Ketson R. M. dos Santos
    % Email: dossantk@umn.edu
    % Date of last modification: May 5, 2024.
    % Version 2.1.0

    clc
    clear results
    close all

    % Turn warnings off.
    warning('off','all');

    % Instantiate the modules.

    % Read input data module.
    data = MSAread;

    % Pre-processing module.
    preprocess = MSApre;

    % Processing module.
    processor = MSAprocessor;

    % Visualization module.
    viewer = MSAviewer;

    verbose = true; % Verbosity.
    homedir = pwd; % Home directory.
    inputdir = pwd; % Input data directory
    MSAtooldir = pwd; % MSAtool directory.
    if numel(varargin)==1
        argin=varargin{1};
    elseif numel(varargin) == 2
        argin=varargin{1};
        verbose=varargin{2};
    elseif numel(varargin) == 3
        argin=varargin{1};
        verbose=varargin{2};
        homedir = varargin{3};
    elseif numel(varargin) == 4
        argin=varargin{1};
        verbose=varargin{2};
        homedir = varargin{3};
        inputdir = varargin{4};
    else
        error('Error: the number of input variables in read_data is larger than expected')
    end

    if isstr(argin)
        input_file=varargin{1};
        model_input = false;
    elseif isstruct(argin)
        model=varargin{1};
        model_input = true;
    else
        error('MSAtool: error in the input arguments.');
    end
    
    % Read input data.
    if ~model_input
        model = data.msa_read_data(input_file,verbose,homedir,inputdir,MSAtooldir);
    end

    % Construct the model.
    model = preprocess.msa_pre(model);

    % Process the model.
    results = processor.msa_processor(model);
    viewer.msa_viewer(model,results);

    % Show the name of the output file.
    if verbose
        if model.output_config.save_output
            fprintf('\n%s%s\n',['* Output file: ' results.output_file]);
        end
    end
    
end










