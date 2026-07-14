function y = gaussian(amplitude, mean, sigma, x)
    y = amplitude * exp(-(x - mean).^2 / (2 * sigma^2));
end