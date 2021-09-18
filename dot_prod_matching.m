function [quant_maps]=dot_prod_matching(dict, acquired_data)
% function [quant_maps]=dot_prod_matching(dict, acquired_data)
% Dot-product matching of MRF data
% Or Perlman 2021 (operlman@mgh.harvard.edu)
% Input:
% ----- %
% dict - MRF dictionry (struct):
%        dict.t1w - water T1
%        dict.t2w - water T2
%        dict.t1s - solute T1
%        dict.t2s - solute T2
%        dict.fs - solute volume fraction
%        dict.ksw - solute exchange rate
%        dict.sig - simulated signals
% raw_data: acquired raw data
%           double array of number_iterations x rows x columns
% Output
% ------ %
%  quant_maps - quantitative maps (struct):
%        quant_maps.t1w - water t1
%        quant_maps.t2w - water T2
%        quant_maps.t1s - solute T1
%        quant_maps.t2s - solute T2
%        quant_maps.fs  - solute volume fraction
%        quant_maps.ksw - solute exchange rate
%        quant_maps.dp  - dot product
%-------------------------------------------%

% Number of shceudule iterations and raw data dimensions
[n_iter, r_raw_data, c_raw_data] = size(acquired_data);

% Reshaping image data to voxel-associated columns
data = reshape(acquired_data, n_iter, r_raw_data * c_raw_data);

% Output quantitative maps, initially as 1D zero-vectors
dp  = zeros(1, r_raw_data * c_raw_data);
t1w = dp;
t2w = dp;
t1s = dp;
t2s = dp;
fs  = dp;
ksw = dp;

% 2-norm normalization
norm_dict = normc(dict.sig);
norm_data = normc(data);

% Matching in batches due to memory considerations
% (can be modified according to available RAM)

batch_size = 64;
assert(size(data, 2)/batch_size == round(size(data, 2)/batch_size), ...
    'The number of image pixel needs to be dividable by batch_size')

batch_indices = 1:batch_size:size(data, 2);

for ind = 1:length(batch_indices)
    
    % Dot-product for the current batch
    current_score = norm_data(:, batch_indices(ind):batch_indices(ind) + batch_size - 1)' *  norm_dict;
    
    % Finding maximum dot-product and storing the corresponding parameters
    [dp(batch_indices(ind):batch_indices(ind) + batch_size - 1), dp_ind] = max(current_score, [], 2);
    t1w(batch_indices(ind):batch_indices(ind) + batch_size - 1) = dict.t1w(dp_ind);
    t2w(batch_indices(ind):batch_indices(ind) + batch_size - 1) = dict.t2w(dp_ind);
    t1s(batch_indices(ind):batch_indices(ind) + batch_size - 1) = dict.t1s(dp_ind);
    t2s(batch_indices(ind):batch_indices(ind) + batch_size - 1) = dict.t2s(dp_ind);
    fs(batch_indices(ind):batch_indices(ind) + batch_size - 1) = dict.fs(dp_ind);
    ksw(batch_indices(ind):batch_indices(ind) + batch_size - 1) = dict.ksw(dp_ind);
    
    % Display progress every 10 percent
    if mod(ind/length(batch_indices).*100, 10) == 0
        disp(['Dot-product matching: ', num2str(ind/length(batch_indices).*100), '% done'])
    end
    
end

% Reshaping the output to the original image dimensions
quant_maps.dp = reshape(dp, r_raw_data, c_raw_data);
quant_maps.t1w = reshape(t1w, r_raw_data, c_raw_data);
quant_maps.t2w = reshape(t2w, r_raw_data, c_raw_data);
quant_maps.fs = reshape(fs, r_raw_data, c_raw_data);
quant_maps.ksw = reshape(ksw, r_raw_data, c_raw_data);

end