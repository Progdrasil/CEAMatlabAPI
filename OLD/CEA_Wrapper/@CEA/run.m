function data = run(obj)
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
    
    if obj.Debug
        c2 = clock;
    end
    % ID OS and run CEA acordingly
    if ismac
        [status,cmdout] = dos(strcat(thermoPath,'/PCEA2.out'));
%         disp(status)
%         disp(cmdout)
    elseif isunix
        [status,cmdout] = dos(strcat(thermoPath,'/PCEA2.out'));
%         disp(status)
%         disp(cmdout)
    elseif ispc
        [status,cmdout] = dos(strcat(thermoPath,'\PCEA2.exe'));
%         disp(status)
%         disp(cmdout)
    else
        disp('Platform not supported')
    end  
    % wait until wrapper.dat exists
    if (status == 0)
        if (cmdout =='')
            fprintf('OK')
            while exist(strcat(obj.ioinp,'.dat'),'file')==0
            end
        else
            error(cmdout);
        end
    end
    if obj.Debug
        c2 = clock - c2;
        fprintf('time to run CEA = %16.15e sec \n',c2(end))
    end
    
    if obj.Debug
        c3 = clock;
    end
    % read output file and create struct
    data = obj.ReadOutput(strcat(obj.ioinp,'.dat'));
    if obj.Debug
        c3 = clock - c3;
        fprintf('time to read output file = %16.15e sec \n',c3(end))
    end
    
    if exist(obj.folder_name,'dir') ~= 7
        mkdir(obj.folder_name);
    end
    
    movefile(strcat(obj.ioinp,'.dat'),strcat(obj.folder_name,'/OF',...
        num2str(obj.OF),'_P',num2str(obj.pressure),'.dat'));
    return;
end