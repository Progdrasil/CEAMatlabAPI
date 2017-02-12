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
        dos(strcat(currentpath,'/PCEA2.out'))
    elseif isunix
        dos('PCEA2.out')
    elseif ispc
        dos('PCEA2.exe')
    else
        disp('Platform not supported')
    end  



% end

