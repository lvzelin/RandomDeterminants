function p=compute_net_pairing_core(sign_table)
load('column_splits_info.mat');
% input with a cell that has info about + - and no specific number
% sign_table{1} =[]; %12
% sign_table{2} = []; %13
% sign_table{3} = []; %14
% sign_table{4} = []; %15
% sign_table{5} = []; %16
% sign_table{6} = []; %23
% sign_table{7} = []; %24
% sign_table{8} = []; %25
% sign_table{9} = []; %26
% sign_table{10} = [0,-1]; %34
% sign_table{11} = []; %35
% sign_table{12} = []; %36
% sign_table{13} = []; %45
% sign_table{14} = []; %46
% sign_table{15} = [0,-1]; %56

% initialization
for i=1:15
    input_info{i} =[];
end

indices_info={'12','13','14','15','16','23','24','25','26','34','35','36','45','46','56'};

% first put placeholding elements in to the table
total_elements=0; % each - corresponding to 4
total_placeholders=0; % each + and - corresponding to 2

for i=1:15
    A=sign_table{i};
    num_neg = sum(A(:) == -1);
    total_elements=total_elements+num_neg*4;
    num_pos = sum(A(:) == 0);
    total_elements=total_elements+num_pos*2;
    total_placeholders=total_placeholders+size(A,2);
end
num_cols=total_elements/6;
placeholder_index=1;
ordering_size=zeros(total_placeholders,1);
for i=1:15
    A=sign_table{i};

    for j = 1:size(A,2)
        input_info{i}=[input_info{i};placeholder_index,A(j)];
        if A(j)==0
            ordering_size(placeholder_index,1)=2;
        elseif A(j)==-1
            ordering_size(placeholder_index,1)=4;
        end
        placeholder_index=placeholder_index+1;
    end
end

canonicalPTable=createCanonicalPositive(input_info);
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
% use this to track how the minus signs are splited
pairing_order={''};
for i =1:max_element
    curr_strings=pairing_info{i};
    new_canonical_forms={};
    new_order={};
    for j=1:size(curr_strings,2)
        curr_s=curr_strings{j};
        j_str=num2str(j);

        for k=1:size(canonical_forms,2)
            c_s=strcat(canonical_forms{k},curr_s);
            new_canonical_forms{end+1}=c_s;

            j_s=strcat(pairing_order{k},j_str);
            new_order{end+1}=j_s;
        end
    end
    canonical_forms=new_canonical_forms;
    pairing_order=new_order;
end

p=0;
size(canonical_forms,2);
for i=1:size(canonical_forms,2)
    s=canonical_forms{i};
    element_split=createElementSplit(s,ordering_size);
    curr_order=pairing_order{i};
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
            curr_sign=calculate_sign(curr_split,element_split,canonicalPTable);
            v=net_dict(key);
            p=p+net_dict(key)*curr_sign;
        end
        
    end
end

% p
end

% create canonical positive table
function T=createCanonicalPositive(info)
% first compute the overall size
total_elements=0;
for i=1:15
    A=info{i};
    num_neg = sum(A(:) == -1);
    total_elements=total_elements+num_neg*4;
    num_pos = sum(A(:) == 0);
    total_elements=total_elements+num_pos*2;
end
num_cols=total_elements/6;
T=zeros(6,num_cols);
indices_info=[1,2;1,3;1,4;1,5;1,6;2,3;2,4;2,5;2,6;3,4;3,5;3,6;4,5;4,6;5,6];
nsix=[1,2,3,4,5,6];
y=1;
% construct table
% track where to put
cur_index=ones(6,1);
for i=1:15
    A=info{i};
    for j=1:size(A,1)
        if A(j,2)==0
            cur_rows=indices_info(i,:);
            cur_col=cur_index(cur_rows(1));
            cur_index(cur_rows(1))=cur_index(cur_rows(1))+1;
            T(cur_rows(1),cur_col)=y;
            second_col=cur_index(cur_rows(2));
            cur_index(cur_rows(2))=cur_index(cur_rows(2))+1;
            T(cur_rows(2),second_col)=y;
            y=y+1;
        end
    end
end

for i=1:15
    A=info{i};
    for j=1:size(A,1)
        if A(j,2)==-1
            cur_rows=setdiff(nsix,indices_info(i,:));

            for k=1:4
                cur_col=cur_index(cur_rows(k));
                cur_index(cur_rows(k))=cur_index(cur_rows(k))+1;
                T(cur_rows(k),cur_col)=y;
            end
            y=y+1;

        end
    end
end
end

function split_string=createElementSplit(s, split)
split_string = cell(1, length(split));  % Initialize the cell array
startIdx = 1;  % Start index for substring extraction

for i = 1:length(split)
    endIdx = startIdx + split(i) - 1;  % End index for substring
    split_string{i} = s(startIdx:endIdx);  % Extract and store substring
    startIdx = endIdx + 1;  % Update the start index for next iteration
end
end

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

function s = matrixToString(M)
% Reshape the matrix into a single row
singleRow = reshape(M.', 1, []);

% Convert each element to a string and concatenate
s = '';
for i = 1:length(singleRow)
    s = strcat(s, num2str(singleRow(i)));
end
s=s(~isspace(s));
end

function su=find_and_delete(s,k)
su=s;
l=length(k);
for i=1:(length(s)/l)
    if strcmp(s((i-1)*l+1:i*l),k)
        su=strcat(s(1:(i-1)*l),s(i*l+1:end));
        break
    end
end
end

% Function to calculate the net number for a given column split S
function overall_perm_sign = calculate_sign(S,E,T)
positive_pairing='';
negative_pairing='';
for i=1:size(E,2)
    pairing=E{i};
    if length(pairing)==2
        positive_pairing=append(positive_pairing,pairing);
    else
        negative_pairing=append(negative_pairing,pairing);
    end
end

% first put positive pairing accordingly then negative pairing
index_dict= containers.Map('KeyType', 'char', 'ValueType', 'char');
for i=1:2:length(positive_pairing)
    key = positive_pairing(i:i+1);
    v=(i+1)/2;
    if ~isKey(index_dict, key)
        index_dict(key) = int2str(v);
    else
        to_store=strcat('|',int2str(v));
        index_dict(key) = append(index_dict(key),to_store);
    end
end

for i=1:4:length(negative_pairing)
    key1 = negative_pairing(i:i+1);
    key2 = negative_pairing(i+2:i+3);
    v=(i+3)/4+length(positive_pairing)/2;
    if ~isKey(index_dict, key1)
        index_dict(key1) = int2str(v);
    else
        to_store=strcat('|',int2str(v));
        index_dict(key1) = append(index_dict(key1),to_store);
    end

    if ~isKey(index_dict, key2)
        index_dict(key2) = int2str(v);
    else
        to_store=strcat('|',int2str(v));
        index_dict(key2) = append(index_dict(key2),to_store);
    end
end

% First construct the corresponding matrix M from S
index_info=stringToMatrix(S);
order_index_info=sortrows(index_info);
M=zeros(6,size(index_info,1)/3);


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
% M

T_sign=1;
for row = 1:size(T, 1)
    perm_sign = sign_of_permutation(T(row, :));
    T_sign = T_sign * perm_sign;
end

m_sign=1;
for row = 1:size(M, 1)
    perm_sign = sign_of_permutation(M(row, :));
    m_sign = m_sign * perm_sign;
end
overall_perm_sign=m_sign*T_sign;
end

% Function to calculate the net number for a given column split S
function overall_perm_sign = calculate_sign_old(S)
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
        if p(i) == p(j)
            continue
        else
            sgn = sgn * sign(p(i) - p(j));
        end
    end
end
end