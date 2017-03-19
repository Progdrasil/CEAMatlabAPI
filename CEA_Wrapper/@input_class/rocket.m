function inp = rocket(obj)%, OF, pressure, presUnit, supar, PcPe, fuel, fuelWt, fuelTemp, oxid, oxidWt, oxidTemp)
    % examples
    % inp = test.rocket( 3 , 350, 'psia', 4.84 , 23.8 , 'paraffin' , 100, 298.15 , 'N2O', 100, 298.15)
    % inp = test.rocket( 3 , 350 ,'psia', 4.84 , 23.8 , {'paraffin' 'CH4' 'RP-1'} , [50 25 25], [298.15 298.15 298.15], {'N2O' 'O2(L)'}, [75 25],[298.15 90.1])
    % verify wt percentages
    if obj.parent.Debug
        c1 = clock;
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
    
    if sum(obj.parent.fuelWt)<=1
        obj.parent.fuelWt = obj.parent.fuelWt*100;
    end
    if sum(obj.parent.oxidWt)<=1
        obj.parent.oxidWt = obj.parent.oxidWt*100;
    end
    
    
    % Open wrapper.inp and overwrite
    if ismac
        inp = strcat(thermoPath,'/wrapper');
    elseif ispc
        inp = strcat(thermoPath,'\wrapper');
    elseif isunix
        inp = strcat(thermoPath,'/wrapper');
    else
        error('Platform not supported')
    end
    IOinp = fopen(strcat(inp,'.inp'),'w');
    
    % Write data in wrapper.inp
    fprintf(IOinp,'prob case=wrapper ro equilibrium \n\n');
    fprintf(IOinp,' ! iac problem \n');
    fprintf(IOinp,'o/f %g\n',obj.parent.OF);   %oxidiser to fuel ratio
    fprintf(IOinp,'p,%s  %g\n',obj.parent.presUnit,obj.parent.pressure); %pressure
    fprintf(IOinp,'supar %g\n',obj.parent.supar);  %supersonic area ratio
    fprintf(IOinp,'pip %g\n',obj.parent.PcPe);    %Pc/Pe
    fprintf(IOinp,'reac\n');
    if length(obj.parent.fuelWt)>1
        for i = 1:length(obj.parent.fuelWt)
            fprintf(IOinp,'  fuel  %s wt%%=%6.3f t,k=%6.2f\n',obj.parent.fuel{i},...
                obj.parent.fuelWt(i),obj.parent.fuelTemp(i));
        end
    else
        fprintf(IOinp,'  fuel  %s wt%%=%g. t,k=%6.2f\n',obj.parent.fuel,...
            obj.parent.fuelWt,obj.parent.fuelTemp);
    end
    if length(obj.parent.oxidWt)>1
        for i = 1:length(obj.parent.oxidWt)
            fprintf(IOinp,'  oxid  %s wt%%=%6.3f t,k=%6.2f\n',obj.parent.oxid{i},...
                obj.parent.oxidWt(i),obj.parent.oxidTemp(i));
        end
    else
        fprintf(IOinp,'  oxid  %s wt%%=%g. t,k=%6.2f\n',obj.parent.oxid,...
            obj.parent.oxidWt,obj.parent.oxidTemp);
    end
    fprintf(IOinp,'output    short\n');
    fprintf(IOinp,'output trace=1e-5\n');
    fprintf(IOinp,'end\n');
    
    % Close wrapper.inp
    fclose(IOinp);
    if obj.parent.Debug
        c1 = clock - c1;
        fprintf('time to write input file = %16.15e sec \n',c1(end))
    end

    obj.parent.ioinp = inp;
    return;
end
