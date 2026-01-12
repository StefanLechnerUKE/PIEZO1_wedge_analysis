function [traces] = LoadDataFromMFXFile_V3(myfile)
    % Load the data
    load(myfile);
    
    % Convert and transpose variables
    tid = double(tid)';
    tim = double(tim)';
    itr = double(itr)';
    efo = efo';
    cfr = cfr';
    fbg = fbg';
    
    % Combine data
    data = [loc tid efo cfr fbg tim itr];
    
    % Find rows where column 9 == 9 (3D valid traces)
    valid_rows = (data(:, 9) == 9);
    
    % Only keep rows where we have at least 3 rows above (skip first 3 valid rows)
    valid_indices = find(valid_rows);
    valid_indices = valid_indices(valid_indices > 3);
    
    % Create traces array in one operation
    traces = data(valid_indices, :);
    
    % Replace column 6 with values from 3 rows above
    traces(:, 6) = data(valid_indices - 3, 6);
    
    % Keep only first 8 columns
    traces = traces(:, 1:8);
end