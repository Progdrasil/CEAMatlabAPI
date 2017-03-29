function setOxid(obj, oxids, oxidWeights, oxidTemps)
    % Function to set oxidizer property's
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