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