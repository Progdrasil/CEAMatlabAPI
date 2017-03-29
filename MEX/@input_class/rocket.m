function inp = rocket(obj)%, OF, pressure, presUnit, supar, PcPe, fuel, fuelWt, fuelTemp, oxid, oxidWt, oxidTemp)
    % Generates Input cell array for mex function
    % examples
    % inp = test.rocket( 3 , 350, 'psia', 4.84 , 23.8 , 'paraffin' , 100, 298.15 , 'N2O', 100, 298.15)
    % inp = test.rocket( 3 , 350 ,'psia', 4.84 , 23.8 , {'paraffin' 'CH4' 'RP-1'} , [50 25 25], [298.15 298.15 298.15], {'N2O' 'O2(L)'}, [75 25],[298.15 90.1])
    if obj.parent.Debug
        c1 = clock;
    end
    i = 1;
    inp{i,1} = sprintf('prob case=wrapper ro equilibrium');
    i = i + 1;
    inp{i,1} = sprintf('');
    i = i + 1;
    inp{i,1} = sprintf(' ! iac problem');
    i = i + 1;
    inp{i,1} = sprintf('o/f %g',obj.parent.OF);   %oxidiser to fuel ratio
    i = i + 1;
    inp{i,1} = sprintf('p,%s  %g',obj.parent.presUnit,obj.parent.pressure); %pressure
    i = i + 1;
    inp{i,1} = sprintf('supar %g',obj.parent.supar);  %supersonic area ratio
    i = i + 1;
    inp{i,1} = sprintf('pip %g',obj.parent.PcPe);    %Pc/Pe
    i = i + 1;
    inp{i,1} = sprintf('reac');
    i = i + 1;
    if length(obj.parent.fuelWt)>1
        for j = 1:length(obj.parent.fuelWt)
            inp{i,1} = sprintf('  fuel  %s wt%%=%6.3f t,k=%6.2f',obj.parent.fuel{j},obj.parent.fuelWt(j),obj.parent.fuelTemp(j));
            i = i + 1;
        end
    else
        inp{i,1} = sprintf('  fuel  %s wt%%=%g. t,k=%6.2f',obj.parent.fuel,obj.parent.fuelWt,obj.parent.fuelTemp);
        i = i + 1;
    end
    if length(obj.parent.oxidWt)>1
        for j = 1:length(obj.parent.oxidWt)
            inp{i,1} = sprintf('  oxid  %s wt%%=%6.3f t,k=%6.2f',obj.parent.oxid{j},obj.parent.oxidWt(j),obj.parent.oxidTemp(j));
            i = i + 1;
        end
    else
        inp{i,1} = sprintf('  oxid  %s wt%%=%g. t,k=%6.2f',obj.parent.oxid,obj.parent.oxidWt,obj.parent.oxidTemp);
        i = i + 1;
    end
    inp{i,1} = sprintf('output    short');
    i = i + 1;
    inp{i,1} = sprintf('output trace=1e-5');
    i = i + 1;
    inp{i,1} = sprintf('end');
    if obj.parent.Debug
        c1 = clock - c1;
        fprintf('time to write input string \t= %16.15e sec \n',c1(end))
    end

    obj.parent.ioinp = inp;
    return;
end
