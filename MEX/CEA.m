classdef CEA < handle

    properties (SetAccess = public)
        data;
        ioinp = 'No inputs created';
        input;
        % variables for Input creation
        OF;
        pressure;
        presUnit = 'psia';
        supar;
        PcPe;

    end
    
    properties (SetAccess = ?input_class)
        fuel;
        fuelWt;
        fuelTemp;

        oxid;
        oxidWt;
        oxidTemp;
    end
    methods
        function this = CEA()
            this.input = input_class(this);
        end
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
            data = cea2(obj.ioinp,inputFile,thermoPath);
            obj.data = data;
            return;
        end
        function setFuel(obj, fuels, fuelWeights, fuelTemps)
            if (iscellstr(fuels))
                if (length(fuelWeights) == length(fuelTemps) && length(fuelTemps) == length(fuels))
                    obj.fuel = fuels;
                    obj.fuelWt = fuelWeights;
                    obj.fuelTemp = fuelTemps;
                else
                    error ('The size of the fuels cell, fuelWeights array and fuelTemps array must match')
                end
            elseif (length(fuelWeights) == 1 && length(fuelTemps) == 1)
                obj.fuel = fuels;
                obj.fuelWt = fuelWeights;
                obj.fuelTemp = fuelTemps;
            else
                error ('if there is only one fuel, the fuelWeights and fuelTemps must be scalar')
            end
        end
        function setOxid(obj, oxids, oxidWeights, oxidTemps)
            if (iscellstr(oxids))
                if (length(oxidWeights) == length(oxidTemps) && length(oxidTemps) == length(oxids))
                    obj.oxid = oxids;
                    obj.oxidWt = oxidWeights;
                    obj.oxidTemp = oxidTemps;
                else
                    error ('The size of the oxids cell, oxidWeights array and oxidTemps array must match')
                end
            elseif (length(oxidWeights) == 1 && length(oxidTemps) == 1)
                obj.oxid = oxids;
                obj.oxidWt = oxidWeights;
                obj.oxidTemp = oxidTemps;
            else
                error ('if there is only one oxid, the oxidWeights and oxidTemps must be scalar')
            end
        end

    end
    
end