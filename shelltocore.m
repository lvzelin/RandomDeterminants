load('column_splits_info.mat');
% current ver: sign are computed in this program, may make a dict before to
% save computational time
tic;

% input with a cell that has info about + - and the corresponding number
% input_info{1} =[1,-1]; %12
% input_info{2} = [2,-1]; %13
% input_info{3} = [3,-1]; %14
% input_info{4} = []; %15
% input_info{5} = []; %16
% input_info{6} = []; %23
% input_info{7} = [4,-1]; %24
% input_info{8} = [5,-1]; %25
% input_info{9} = []; %26
% input_info{10} = [6,-1]; %34
% input_info{11} = []; %35
% input_info{12} = [7,-1]; %36
% input_info{13} = []; %45
% input_info{14} = []; %46
% input_info{15} = [8,-1;9,-1]; %56

input_info{1} =[1,-1]; %12
input_info{2} = []; %13
input_info{3} = []; %14
input_info{4} = [4,0;5,-1]; %15
input_info{5} = []; %16
input_info{6} = []; %23
input_info{7} = []; %24
input_info{8} = []; %25
input_info{9} = []; %26
input_info{10} = [2,-1]; %34
input_info{11} = []; %35
input_info{12} = []; %36
input_info{13} = []; %45
input_info{14} = []; %46
input_info{15} = [3,-1]; %56

indices_info={'12','13','14','15','16','23','24','25','26','34','35','36','45','46','56'};

% assumed to start from 1 (to_do: get the min)
% convert the -1 (-) into possible 0 (+) places, and get the canonical forms
max_element=0;
for i = 1:15
    curr_matrix=input_info{i};
    if size(curr_matrix,1)==0
        continue
    end
    max_element=max(max(curr_matrix, [], 'all'),max_element);
end

pairing_info=cell(max_element,1);

for i = 1:15
    curr_matrix=input_info{i};
    if size(curr_matrix,1)==0
        continue
    end
    for j=1:size(curr_matrix,1)
        curr_row=curr_matrix(j,:);
        if curr_row(2)==0
            pairing_info{curr_row(1)}=indices_info(i);
        else
            curr_row(1);
            string_index=indices_info(i);
            curr_string=string_index{1};
            pairing_info{curr_row(1)}=getTwoTwoPartitions(curr_string(1),curr_string(2));
        end
    end
end

canonical_forms={''};
for i =1:max_element
    curr_strings=pairing_info{i};
    new_canonical_forms={};

    for j=1:size(curr_strings,2)
        curr_s=curr_strings{j};
        for k=1:size(canonical_forms,2)
            c_s=strcat(canonical_forms{k},curr_s);
            new_canonical_forms{end+1}=c_s;
        end
    end
    canonical_forms=new_canonical_forms;
end
p=0;
size(canonical_forms,2)
for i=1:size(canonical_forms,2)
    s=canonical_forms{i};
    m=stringToMatrix(s);
    m=sortrows(m);
    key = reshape(num2str(m'), 1, []);
    key = key(~isspace(key));
    if isKey(net_dict, key)
        all_cols=col_dict(key);
        num_col_splits=length(all_cols)/length(s);
        for j=1:num_col_splits
            start=(j-1)*length(s)+1;
            endp=j*length(s);
            curr_split=all_cols(start:endp);
            curr_sign=calculate_sign(curr_split);
            v=net_dict(key);
            p=p+net_dict(key)*curr_sign;
        end
        
    end
end

elapsedTime = toc;
fprintf('Total net number: %d.\n', p);
fprintf('Total running time: %.6f seconds.\n', elapsedTime);

function partitions = getTwoTwoPartitions(a, b)
    % Get the four numbers that are different from a and b
    all_nums = ['1','2','3','4','5','6'];
    remaining_nums = setdiff(all_nums, [a, b]);

    % Generate all 2-2 partitions of the remaining numbers without repetitions
    partitions = {};
    for i = 1:3
        s=strcat(remaining_nums(1),remaining_nums(i+1));
        reremaining_nums = setdiff(remaining_nums, [remaining_nums(1),remaining_nums(i+1)]);
        s=strcat(s,reremaining_nums(1),reremaining_nums(2));
        partitions{end+1}=s;
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
