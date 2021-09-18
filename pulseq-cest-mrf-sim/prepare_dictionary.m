%% allocate space for dictionary variables and generate unique combination

% Input:  dict:    dictionary variable
%
% Output: dict:    dictionary variable, filled with all combinations
%         numComb: total number of all combinations

function [dict, numComb] = prepare_dictionary(dict)

% calculate  number of combinations
numComb = 1;
varNames = fieldnames(dict.variables);
numVars = numel(varNames);
C = cell(numVars,1);
for comb = 1:numVars
    cVar = dict.variables.(varNames{comb});
    numComb = numComb * numel(cVar);
    C{comb} = cVar;
end
disp(['Found ' num2str(numComb) ' different parameter combinations.']);

% reshape combinations
n = length(C);
[C{:}] = ndgrid(C{:});
combinations = reshape(cat(n+1,C{:}),[],n);

% put all combinations in dict
for comb = 1:numVars
    dict.(varNames{comb}) = combinations(:,comb);
end

