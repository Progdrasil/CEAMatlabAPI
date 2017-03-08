classdef input_class < handle
    
    properties (SetAccess = private)
        inp;
        parent;
    end
    
    methods (Access = ?CEA)
        function this = input_class(parent)
            this.inp = 'No inputs created';
            this.parent = parent;
            this.parent.ioinp = this.inp;
        end
    end
    methods
        function input = rocket(obj, OF,pressure,presUnit, supar, PcPe, fuel, fuelWt, fuelTemp, oxid, oxidWt, oxidTemp)
            % examples
            % input = test.rocket( 3 , 350, 'psia', 4.84 , 23.8 , 'paraffin' , 100, 298.15 , 'N2O', 100, 298.15)
            % input = test.rocket( 3 , 350 ,'psia', 4.84 , 23.8 , {'paraffin' 'CH4' 'RP-1'} , [50 25 25], [298.15 298.15 298.15], {'N2O' 'O2(L)'}, [75 25],[298.15 90.1])
            i = 1;
            input{i,1} = sprintf('prob case=wrapper ro equilibrium');
            i = i + 1;
            input{i,1} = sprintf('');
            i = i + 1;
            input{i,1} = sprintf(' ! iac problem');
            i = i + 1;
            input{i,1} = sprintf('o/f %g',OF);   %oxidiser to fuel ratio
            i = i + 1;
            input{i,1} = sprintf('p,%s  %g',presUnit,pressure); %pressure
            i = i + 1;
            input{i,1} = sprintf('supar %g',supar);  %supersonic area ratio
            i = i + 1;
            input{i,1} = sprintf('pip %g',PcPe);    %Pc/Pe
            i = i + 1;
            input{i,1} = sprintf('reac');
            i = i + 1;
            if length(fuelWt)>1
                for j = 1:length(fuelWt)
                    input{i,1} = sprintf('  fuel  %s wt%%=%6.3f t,k=%6.2f',fuel{j},fuelWt(j),fuelTemp(j));
                    i = i + 1;
                end
            else
                input{i,1} = sprintf('  fuel  %s wt%%=%g. t,k=%6.2f',fuel,fuelWt,fuelTemp);
                i = i + 1;
            end
            if length(oxidWt)>1
                for j = 1:length(oxidWt)
                    input{i,1} = sprintf('  oxid  %s wt%%=%6.3f t,k=%6.2f',oxid{j},oxidWt(j),oxidTemp(j));
                    i = i + 1;
                end
            else
                input{i,1} = sprintf('  oxid  %s wt%%=%g. t,k=%6.2f',oxid,oxidWt,oxidTemp);
                i = i + 1;
            end
            input{i,1} = sprintf('output    short');
            i = i + 1;
            input{i,1} = sprintf('output trace=1e-5');
            i = i + 1;
            input{i,1} = sprintf('end');
            
            obj.inp = input;
            obj.parent.ioinp = input;
            return;
        end

    end
end