function setFuel(obj, fuels, fuelWeights, fuelTemps)
    % CEA.setFuel Function to set fuel properties
    %   Assures a cohesif setup of all the fuel properties of the current
    %   CEA objet. All properties of the fuel must have the same index of
    %   the corresponding fuel in it's cell array. Each input vector must 
    %   also have the same size.
    %
    % CEA.setFuel Examples
    %   CEAobj = CEA;
    %   CEAobj.setFuel('paraffin',100,298.15);
    %   CEAobj.setFuel({'paraffin' 'CH4' 'RP-1'} , [50 25 25], [298.15
    %   298.15 298.15]);
    %
    % See also:
    % CEA, CEA.setOxid
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