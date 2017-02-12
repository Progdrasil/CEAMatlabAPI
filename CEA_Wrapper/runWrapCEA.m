% function [ output_args ] = runWrapCEA( input_args )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Takes in as input information to write in wrapper.inp file
    % write wrapper.inp file
    % identify operating system
    % run PCEA2.out or PCEA2.exe depending on OS
    % open wrapper.out file created by CEA
    % read data in wrapper.out
    % close wrapper.out
    % return data in a table or struct 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Open wrapper.inp and overwrite
    IOinp = fopen('wrapper.inp','w');
    
    % Write data in wrapper.inp
    fprintf(IOinp,'prob case=wrapper ro equilibrium \n\n');
    fprintf(IOinp,' ! iac problem \n');
    fprintf(IOinp,'p,bar  7000\n');
    fprintf(IOinp,'reac\n');
    fprintf(IOinp,'  fuel  N2H4(L) wt%%=100. t,k=298.15\n');
    fprintf(IOinp,'  oxid  N2H4(L) wt%%=100. t,k=298.15\n');
    fprintf(IOinp,'output    short\n');
    fprintf(IOinp,'output trace=1e-5\n');
    fprintf(IOinp,'end\n');
    
    % Close wrapper.inp
    fclose(IOinp);

    % ID OS and run CEA acordingly
    if ismac
        currentpath = pwd;
        [status,cmdout] = dos(strcat(currentpath,'/PCEA2.out'));
        disp(status)
        disp(cmdout)
    elseif isunix
        currentpath = pwd;
        [status,cmdout] = dos(strcat(currentpath,'/PCEA2.out'));
        disp(status)
        disp(cmdout)
    elseif ispc
        [status,cmdout] = dos('PCEA2.exe');
        disp(status)
        disp(cmdout)
    else
        disp('Platform not supported')
    end  

    IOout = fopen('wrapper.dat','r');
    tline = fgets(IOout) ;
    
    while ischar(tline)
            
        %Finds the line containing the number of the harmonic
        if any(strfind(tline,'NOHARM(I)'),1)
            disp(tline)
            found = tline ;
            ii = 0 ;
            % Separates the line in strings
            tempcell = strread(found, '%s', 'delimiter','=') ;
            out.noharm = tempcell{end}(1) ;
        end

%         if any(strfind(tline,'NO.          OMEGA          RAD./SEC.      HERTZ'),1)
%             tline = fgets(fid);
%             tline = fgets(fid);
%             C_data = textscan(fileID,,'%d %f %f %f');
%             C_data(1);
%         end

        tline = fgets(IOout);
    end
    
    fclose(IOout);
% end

