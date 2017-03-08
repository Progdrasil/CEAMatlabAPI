classdef CEA

    properties (SetAccess = private)
        data
        input

    end
    methods 
        function data = run(obj)
            [pathstr,~,~] = fileparts(currentpath);
            obj.data = dos(strcat(pathstr,'/PCEA2.out'));
            
        end
        function data = get.data(obj)
            data = obj.data;
        end
        function input = rocket(obj,OF,pressure,presUnit, supar, PcPe ...
                , fuel, fuelWt, fuelTemp, oxid, oxidWt, oxidTemp)
            input = sprintf(IOinp,'prob case=wrapper ro equilibrium \n\n');
            input = strcat(input,fprintf(IOinp,' ! iac problem \n'));
            input = strcat(input,fprintf(IOinp,'o/f %g\n',OF));   %oxidiser to fuel ratio
            input = strcat(input,fprintf(IOinp,'p,%s  %g\n',presUnit,pressure)); %pressure
            input = strcat(input,fprintf(IOinp,'supar %g\n',supar));  %supersonic area ratio
            input = strcat(input,fprintf(IOinp,'pip %g\n',PcPe));    %Pc/Pe
            input = strcat(input,fprintf(IOinp,'reac\n'));
            if length(fuelWt)>1
                for i = 1:length(fuelWt)
                    input = strcat(input,fprintf(IOinp,'  fuel  %s wt%%=%6.3f t,k=%6.2f\n',fuel{i},fuelWt(i),fuelTemp(i)));
                end
            else
                input = strcat(input,fprintf(IOinp,'  fuel  %s wt%%=%g. t,k=%6.2f\n',fuel,fuelWt,fuelTemp));
            end
            if length(oxidWt)>1
                for i = 1:length(oxidWt)
                    input = strcat(input,fprintf(IOinp,'  oxid  %s wt%%=%6.3f t,k=%6.2f\n',oxid{i},oxidWt(i),oxidTemp(i)));
                end
            else
                input = strcat(input,fprintf(IOinp,'  oxid  %s wt%%=%g. t,k=%6.2f\n',oxid,oxidWt,oxidTemp));
            end
            input = strcat(input,fprintf(IOinp,'output    short\n'));
            input = strcat(input,fprintf(IOinp,'output trace=1e-5\n'));
            input = strcat(input,fprintf(IOinp,'end\n'));
            
            obj.input = input;
            return;
        end

    end
    
end