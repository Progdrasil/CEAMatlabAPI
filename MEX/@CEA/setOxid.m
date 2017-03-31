function setOxid(obj, oxids, oxidWeights, oxidTemps)
    % CEA.setOxid Function to set Oxidizer properties
    %   Assures a cohesif setup of all the oxidizer properties of the current
    %   CEA objet. All properties of the oxidizer must have the same index of
    %   the corresponding oxidizer in it's cell array. Each input vector must 
    %   also have the same size.
    %
    % CEA.setOxid Examples
    %   CEAobj = CEA;
    %   CEAobj.setOxid('N2O', 100, 298.15);
    %   CEAobj.setOxid({'N2O' 'O2(L)'}, [75 25],[298.15 90.1]);
    %
    % See also:
    % CEA, CEA.setFuel
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