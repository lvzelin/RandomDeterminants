one_dict = createDictionary();
[col_dict,net_dict]=create_column_dict(one_dict);
keysList = keys(col_dict);
% save('column_splits_info.mat');

% s='121314252635364546'
col_dict('121314252635364546')
s='123456132546'
calculate_sign(s);
% compute_net_number(m)

for i = 1:length(keysList)
    if length(keysList{i})==18
        disp(['Key: ', keysList{i}, ', Value: ', num2str(net_dict(keysList{i}))]);
    end
    % disp(['Key: ', keysList{i}, ', Value: ', num2str(net_dict(keysList{i}))]);
end

function [col_dict,net_dict]=create_column_dict(one_dict)
keysList = keys(one_dict);
keysList{end+1}=' ';
col_dict = containers.Map('KeyType', 'char', 'ValueType', 'char');
net_dict= containers.Map('KeyType', 'char', 'ValueType', 'double');
for i=1:16
    for j=i:16
        for k=j:16
            for l=k:16
                for h=l:16
                    for q=h:16
                        current_string=strcat(keysList{i},keysList{j},keysList{k},keysList{l},keysList{h},keysList{q});
                        current_string = current_string(~isspace(current_string)); % Remove spaces from the key
                        if isempty(current_string)
                            break
                        end
                        matrix_form=stringToMatrix(current_string);
                        net_number=compute_net_number(matrix_form);
                        matrix_form=sortrows(matrix_form);

                        key = reshape(num2str(matrix_form'), 1, []);
                        key = key(~isspace(key));  % Remove spaces from the key

                        if ~isKey(col_dict, key)
                            col_dict(key) = current_string;
                            net_dict(key) = net_number;
                        else
                            col_dict(key) = append(current_string,col_dict(key));
                        end
                    end
                end
            end
        end
    end
end
end

function result = removeMatchingRows(cell_matrices, current_row)
% Initialize the result cell array with the same size as the input
result = cell(size(cell_matrices));

for i = 1:numel(cell_matrices)
    matrix = cell_matrices{i};
    matching_rows = ismember(matrix, current_row, 'rows');
    matrix(matching_rows, :) = [];
    result{i} = matrix;
end
end

% First we compute the number of different column placement, after placing
% all the columns, we can swap two elements if they have same pairing
function net_number = compute_net_number(M)
% Check the dimensions of M
[m, n] = size(M);
if mod(m, 3) ~= 0 || n ~= 2
    error('The input matrix M must have size (3k, 2)');
end

k = m / 3;

% Split the matrix M into k sub-matrices of size (3,2)
cell_matrices = cell(1, k);
for i = 1:k
    start_idx = (i-1)*3 + 1;
    end_idx = i*3;
    cell_matrices{i} = M(start_idx:end_idx, :);
end

% Reshape each (3,2) matrix into a row of 6 elements and concatenate
all_matrices_as_rows = cellfun(@(x) reshape(x', 1, []), cell_matrices, 'UniformOutput', false);
concatenated_rows = vertcat(all_matrices_as_rows{:});
% Sort the matrix to group identical rows together
m_sorted = sortrows(M);
row_m_sorted=sortrows(concatenated_rows);

% how many different columns we can have
current_row = row_m_sorted(1,:);
count = 1;
results = [];

for i = 2:size(row_m_sorted, 1)
    if isequal(row_m_sorted(i,:), current_row)
        % If current row is the same as the previous one, increase count
        count = count + 1;
    else
        % If current row is different, store the count in results array
        results = [results, count];
        % Move to the next row
        current_row = m_sorted(i,:);
        count = 1; % reset count
    end
end
% Check for the last group of rows
results = [results, count];

net_number = multinomial(results);

% Now, different elements placement
current_row = m_sorted(1,:);
count = 1;
results = [];

for i = 2:size(m_sorted, 1)
    if isequal(m_sorted(i,:), current_row)
        % If current row is the same as the previous one, increase count
        count = count + 1;
    else
        % If current row is different, store the count in results array
        results = [results, count];
        % Move to the next row
        current_row = m_sorted(i,:);
        count = 1; % reset count
    end
end
% Check for the last group of rows
results = [results, count];
for i=1:length(results)
    net_number=net_number*results(i);
end
end

function t = multinomial(M)
% Input:
% M - a vector of non-negative integers [m_1, m_2, ..., m_k]
%
% Output:
% t - the computed multinomial coefficient

% Check if M is a vector
if ~isvector(M)
    error('Input must be a vector');
end

% Check if elements of M are non-negative integers
if any(M < 0) || any(mod(M, 1) ~= 0)
    error('All elements of input vector must be non-negative integers');
end

% Calculate n, the sum of elements in M
n = sum(M);

% Calculate n!
n_factorial = factorial(n);

% Calculate the product of factorials of the elements in M
m_factorials_product = prod(factorial(M));

% Calculate the multinomial coefficient
t = n_factorial / m_factorials_product;
end

function result=countDuplicates(m)
% Sort the matrix to group identical rows together
m_sorted = sortrows(m);

% Initialize variables
current_row = m_sorted(1,:);
count = 1;
results = [];

for i = 2:size(m_sorted, 1)
    if isequal(m_sorted(i,:), current_row)
        % If current row is the same as the previous one, increase count
        count = count + 1;
    else
        % If current row is different, check if count is greater than 1
        % If yes, store the count in results array
        if count > 1
            results = [results, count];
        end
        current_row = m_sorted(i,:);
        count = 1;
    end
end
% Check for the last group of rows
if count > 1
    results = [results, count];
end

% Get unique counts and their frequencies
unique_counts = unique(results);
result = zeros(1, 6);
for i = 1:6
    result(i) = sum(results == i);
end

end

function dict = createDictionary()
% All possible combinations of 2 elements from 1 to 6
combList = nchoosek(1:6, 2);

% Create the dictionary
dict = containers.Map('KeyType', 'char', 'ValueType', 'double');

% Iterate over combinations and form matrices
for i = 1:size(combList, 1)
    for j = 1:size(combList, 1)
        for k = 1:size(combList, 1)
            if isempty(intersect(combList(i, :), combList(j, :))) && ...
                    isempty(intersect(combList(i, :), combList(k, :))) && ...
                    isempty(intersect(combList(j, :), combList(k, :)))

                matrix = cat(1, sort(combList(i, :)), sort(combList(j, :)), sort(combList(k, :)));
                matrix = sortrows(matrix, 1);

                % Convert matrix to a unique string key in row order
                key = reshape(num2str(matrix'), 1, []);
                key = key(~isspace(key));  % Remove spaces from the key

                % Add to dictionary with value 1 (only if key is new)
                if ~isKey(dict, key)
                    dict(key) = 1;
                end
            end
        end
    end
end

end

function matrix = stringToMatrix(s)
% if length(s) ~= 6
if mod(length(s),2)~=0
    error('Input string should have a length of 6');
end

% Convert string to double array
numArray = double(s) - double('0');  % Subtracting the ASCII value of '0' to get numeric values

% Reshape to x2 matrix
matrix = reshape(numArray, [2, length(s)/2])';
end


% Function to calculate the net number for a given column split S
function overall_perm_sign = calculate_sign(S)
% First construct the corresponding matrix M from S
index_info=stringToMatrix(S);
order_index_info=sortrows(index_info);
M=zeros(6,size(index_info,1)/3);
index_dict= containers.Map('KeyType', 'char', 'ValueType', 'char');
for i=1:size(index_info,1)
    key = reshape(num2str(order_index_info(i,:)), 1, []);
    key = key(~isspace(key));

    if ~isKey(index_dict, key)
        index_dict(key) = int2str(i);
    else
        to_store=strcat('|',int2str(i));
        index_dict(key) = append(index_dict(key),to_store);
    end
end

for i=1:size(index_info,1)/3
    for j=1:3
        key = extractBetween(S,(i-1)*6+(j-1)*2+1,(i-1)*6+(j-1)*2+2);
        key=key{1};
        curr_info=index_dict(key);

        separatorIndex = find(curr_info == '|', 1);
        % If a separator was found, extract the first number
        if ~isempty(separatorIndex)
            firstNumberStr = curr_info(1:separatorIndex-1);
            curr_info = curr_info(separatorIndex+1:end);
        else
            % If there is no separator, the whole string is the first number
            firstNumberStr = curr_info;
        end
        index_dict(key)=curr_info;

        M(index_info((i-1)*3+j,:),i)=str2num(firstNumberStr);

    end

end
    
M
% Initialize overall permutation sign
overall_perm_sign = 1;
% Calculate permutation sign for each row and update overall sign
for row = 1:size(M, 1)
    perm_sign = sign_of_permutation(M(row, :));
    overall_perm_sign = overall_perm_sign * perm_sign;
end
end

% Function to calculate the sign of a permutation
function sgn = sign_of_permutation(p)
sgn = 1;  % Initialize sign
n = length(p);
for i = 1:n
    for j = i+1:n
        sgn = sgn * sign(p(i) - p(j));
    end
end
end