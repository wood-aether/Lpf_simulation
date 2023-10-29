function mse = calculateMSE(actual, predicted)
    mse = mean((actual - predicted).^2);
end
