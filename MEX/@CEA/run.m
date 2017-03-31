function data = run(obj)
    % CEA.run Runs the compiled mex function
    %   Must have ioinp already generated to run as well as having the
    %   files thermo.lib and trans.lib in the current matlab path. Returns
    %   a struct with all the data in the usual way the output files were
    %   generated on the cearun website.
    %
    % CEA.run Examples
    %   CEAobj = CEA;
    %   data = CEAobj.run();
    %
    % See also:
    % CEA, input_class
    if obj.Debug
        c2 = clock;
    end
    try
        thermoFind = which('thermo.lib');
        [thermoPath,~,~] = fileparts(thermoFind);
    catch
        error('Could not find thermo.lib in current Path')
    end
    try
        transFind = which('trans.lib');
        [transPath,~,~] = fileparts(transFind);
    catch
        error('Could not find trans.lib in current Path')
    end
    if (thermoPath ~= transPath)
        error('thermo.lib and trans.lib must be in the same directory')
    end

    inputFile = 'wrapper';
    if ismac
        inputFile = strcat('/',inputFile);
    elseif isunix
        inputFile = strcat('/',inputFile);
    elseif ispc
        inputFile = strcat('\',inputFile);
    else
        error('Platform not supported')
    end  
    data = cea2(obj.ioinp,inputFile,thermoPath);
    if obj.Debug
        c2 = clock - c2;
        fprintf('\t\ttime to run CEA \t= %16.15e sec \n',c2(end))
    end
    obj.data = data;
    return;
end