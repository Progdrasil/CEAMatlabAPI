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
    
    %%% Debug variables
    OF = 3;
    pression = 350;
    presUnt = 'psia';
    supar = 4.84;
    PcPe = 23.8;
    fuel = 'paraffin';
    fuelWt = '100';
    oxid = 'N2O';
    oxidWt = '100';
    
    % Open wrapper.inp and overwrite
    IOinp = fopen('wrapper.inp','w');
    
    % Write data in wrapper.inp
    fprintf(IOinp,'prob case=wrapper ro equilibrium \n\n');
    fprintf(IOinp,' ! iac problem \n');
    fprintf(IOinp,'o/f %g\n',OF);   %oxidiser to fuel ratio
    fprintf(IOinp,'p,psia  %g\n',pression); %pressure
    fprintf(IOinp,'supar %g\n',supar);  %supersonic area ratio
    fprintf(IOinp,'pip %g\n',PcPe);    %Pc/Pe
    fprintf(IOinp,'reac\n');
    fprintf(IOinp,'  fuel  %s wt%%=%g. t,k=298.15\n',fuel,fuelWt); 
    fprintf(IOinp,'  oxid  %s wt%%=%g. t,k=298.15\n',oxid,oxidWt);
    fprintf(IOinp,'output    short\n');
    fprintf(IOinp,'output trace=1e-5\n');
    fprintf(IOinp,'end\n');
    
    % Close wrapper.inp
    fclose(IOinp);

    % ID OS and run CEA acordingly
    if ismac
        currentpath = which('runWrapCEA.m');
        [pathstr,name,ext] = fileparts(currentpath);
        [status,cmdout] = dos(strcat(pathstr,'/PCEA2.out'));
        disp(status)
        disp(cmdout)
    elseif isunix
        currentpath = which('runWrapCEA.m');
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

