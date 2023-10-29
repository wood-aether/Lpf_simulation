% calculateMAE.m
function mae = calculateMAE(actual, predicted)
    mae = mean(abs(actual - predicted));
end

