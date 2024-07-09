% given a table of marks, generate all possible partial shells
% even some variables are named latex_something, it might not be in the
% latex formats
% for mu_5=10 mu_3

% clear
tic

% 3-mark
% table_input_list{1}=[1,1,1,0,0,0]';
% table_input_list{2}=[1,1,2,0,0,0]';
% table_input_list{3}=[1,2,3,0,0,0]';
% table_input_list{4}=[1,1,0,0,0,0;0,0,1,0,0,0]';
% table_input_list{5}=[1,1,0,0,0,0;0,0,2,0,0,0]';
% table_input_list{6}=[1,2,0,0,0,0;0,0,1,0,0,0]';
% table_input_list{7}=[1,2,0,0,0,0;0,0,3,0,0,0]';
% table_input_list{8}=[1,0,0,0,0,0;0,2,0,0,0,0;0,0,3,0,0,0]';

table_input_list{1}=[1,1,1,1,0,0]';
table_input_list{2}=[1,1,1,2,0,0]';
table_input_list{3}=[1,1,2,3,0,0]';
table_input_list{4}=[1,2,3,4,0,0]';
table_input_list{5}=[1,1,0,0,0,0;0,0,1,1,0,0]';
table_input_list{6}=[1,1,0,0,0,0;0,0,2,2,0,0]';
table_input_list{7}=[1,2,0,0,0,0;0,0,1,2,0,0]';
table_input_list{8}=[1,1,0,0,0,0;0,0,2,3,0,0]';
table_input_list{9}=[1,2,0,0,0,0;0,0,1,3,0,0]';
table_input_list{10}=[1,1,2,0,0,0;0,0,0,3,0,0]';
table_input_list{11}=[1,2,3,0,0,0;0,0,0,1,0,0]';
table_input_list{12}=[1,2,3,0,0,0;0,0,0,4,0,0]';
table_input_list{13}=[1,2,0,0,0,0;0,0,3,4,0,0]';

table_input_list{14}=[1,1,0,0,0,0;0,0,2,0,0,0;0,0,0,3,0,0]';
table_input_list{15}=[1,2,0,0,0,0;0,0,3,0,0,0;0,0,0,4,0,0]';
table_input_list{16}=[1,2,0,0,0,0;0,0,1,0,0,0;0,0,0,3,0,0]';
table_input_list{17}=[1,0,0,0,0,0;0,1,0,0,0,0;0,0,2,3,0,0]';
table_input_list{18}=[1,1,0,0,0,0;0,0,1,0,0,0;0,0,0,2,0,0]';
table_input_list{19}=[1,2,0,0,0,0;0,0,1,0,0,0;0,0,0,1,0,0]';
table_input_list{20}=[1,0,0,0,0,0;0,2,0,0,0,0;0,0,3,0,0,0;0,0,0,4,0,0]';

table_input_list{21}=[1,1,1,0,0,0;0,0,0,2,0,0]';
table_input_list{22}=[1,1,2,0,0,0;0,0,0,1,0,0]';
table_input_list{23}=[1,1,2,2,0,0]';
% togolist=[5,6,7,8,9,10,11,14,16,17,23]
togolist=[1,2,3,4]
for overall_i_index=1:size(togolist,2)
% for overall_i_index=1:23
overall_i=togolist(overall_i_index)
table_input=table_input_list{overall_i};
filePath = sprintf('./results/4_mark/case%d.txt',overall_i);

num_occurrences=zeros(max(table_input(:)),1);
max_element=max(table_input(:));
existing_rows = cell(max_element,1);
target_rows = cell(max_element,1);
for i=1: max_element
    num_occurrences(i)=sum(table_input(:) == i);
    % the rows with i
    existing_rows{i}=find(any(table_input == i, 2));
    target_rows{i}=setdiff([1,2,3,4,5,6],existing_rows{i});
end

odd_columns=countColumnsWithOddZeros(table_input);
num_triples=sum(num_occurrences(:)==3)+sum(num_occurrences(:)==1);
num_extra_columns=num_triples-odd_columns;
if mod(num_extra_columns,2)==1 || num_extra_columns<0
    'EXCEPTION: extra column'
else
    num_extra_columns=num_extra_columns/2;
end
extra_columns = zeros(6, num_extra_columns);
table_input = [table_input, extra_columns];

generated_second_stage_shell={table_input};

% first deal with elements covered by two marks
for i=1: max_element
    if num_occurrences(i)==2
        new_second_stage_shell={};
        for j=1:size(generated_second_stage_shell,2)
            curr_table=generated_second_stage_shell{j};
            new_table_1=curr_table;
            new_table_2=curr_table;
            new_column=zeros(6, 1);
            new_column(target_rows{i}) = -i;
            new_table_1 = [new_table_1, new_column];
            new_second_stage_shell{end+1}=new_table_1;

            zero_counts = sum(curr_table == 0);
            % Find the column indices where the count of zeros is exactly four
            column_indices = find(zero_counts == 4);
            rows_with_i = any(curr_table == i, 2);
            % Get the row indices where element i appears
            row_indices = find(rows_with_i);
            % we check if in each of the columns with four zeros, we can
            % place a known four-column
            for k=1:size(column_indices,2)
                curr_col=curr_table(:,column_indices(k));
                zero_rows=any(curr_col == 0, 2);
                zero_row_indices = find(zero_rows);
                if isequal(target_rows{i},zero_row_indices')
                    new_table_2(target_rows{i}, column_indices(k)) = -i;
                    new_second_stage_shell{end+1}=new_table_2;
                end
            end
        end
        generated_second_stage_shell=[generated_second_stage_shell,new_second_stage_shell];
    end
end

% put triples
for i=1:max_element
    new_second_stage_shell={};
    if num_occurrences(i)==1 || num_occurrences(i)==3
        % j represents the current partial shell
        for j=1:size(generated_second_stage_shell,2)
            curr_table=generated_second_stage_shell{j};

            % check columns one by one. If no columns can fit this element
            % i, we just delete this shell
            for k=1:size(curr_table,2)

                zero_rows=find(curr_table(:, k) == 0);
                if size(zero_rows,1)==1||size(zero_rows,1)==2 || size(zero_rows,1)==4
                    continue
                end
                possible_rows=intersect(target_rows{i},zero_rows);

                possible_three_rows=nchoosek(possible_rows, 3);

                for rows_idx = 1:size(possible_three_rows, 1)
                    subset = possible_three_rows(rows_idx, :);
                    new_table=curr_table;
                    new_table(subset, k) = -i;  % Update the k-th column for the rows in the subset
                    new_second_stage_shell{end+1}=new_table;
                end
            end
        end
        generated_second_stage_shell=new_second_stage_shell;
    end
end

shell_dict = containers.Map('KeyType', 'char', 'ValueType', 'char');
shell_weight_dict = containers.Map('KeyType', 'char', 'ValueType', 'double');
% display
for j=1:size(generated_second_stage_shell,2)
    generated_second_stage_shell{j};
end
size_second_stage_shells=size(generated_second_stage_shell,2);

for j=1:size(generated_second_stage_shell,2)
    % if mod(j,10)==0
    %     j
    % end
    
    % if ~ismember(j,want_list)
    %     continue
    % end

    curr_second_stage_table=generated_second_stage_shell{j};
    [mark_table_latex_string,GF_key]=latex_output(curr_second_stage_table);
    % now the labeled information can be deleted
    curr_second_stage_table=abs(curr_second_stage_table);
    % assign the pairing the generate partial shell, in the output negative
    % elements mean unlabeled
    partial_shells=generate_partial_shell(curr_second_stage_table);

    pssize=size(partial_shells,2);

    % first we create the shell_info and shell_table for each partial shell
    % in the list of partial shells
    for i=1:size(partial_shells,2)
        % sort first, then use the minus-plus table as its key
        string_canonical_shell=sortrows(partial_shells{i});
        input_info=getInputInfo(string_canonical_shell);
        key=append(GF_key,cell_to_string(input_info));
        key = strrep(key, '[', '');
        key = strrep(key, ']', '');
        if ~isKey(shell_dict, key)
            % input_info=getInputInfo(partial_shells{i});
            partial_shell_string=compute_factors_partial(input_info);
            output_string=append(mark_table_latex_string,partial_shell_string);
            shell_dict(key) = output_string;
            shell_weight_dict(key) = 1;
        else
            partial_shell_string=shell_dict(key);
            shell_weight_dict(key) = shell_weight_dict(key)+1;
        end

    end

end

% output
fileID = fopen(filePath, 'w');
    if fileID == -1
        error('Failed to open the file.');
    end
keysList = keys(shell_dict);
fprintf(fileID, '(');
for i =1:length(keysList)
    partial_shell_string=shell_dict(keysList{i});
    fprintf(fileID, '%d(%s)',shell_weight_dict(keysList{i}),partial_shell_string);
    if i==length(keysList)
        fprintf(fileID, ')\n');
    else
        fprintf(fileID, '+\n');
    end
end

fclose(fileID);
end

toc

function s=get_canonical_shell(shell)
    s=sortrows(shell);
end

% I think there is no need to check if the generated shells are actually
% partial shells. This property is guareented.

function partial_shells=generate_partial_shell(ss_table)
partial_shells={ss_table};
l=-1;

for i=1:size(ss_table,2)
    new_shells={};
    four_col_flag=false;
    two_col_flag=false;
    for j=1:size(partial_shells,2)
        curr_table=partial_shells{j};
        curr_col=curr_table(:,i);
        if sum(curr_col == 0)==4
            row_partitions=four_rows_partition(curr_col);
            for k=1:3
                new_col=curr_col;
                new_table=curr_table;
                curr_part=row_partitions{k};
                new_col(curr_part(1,:)')=l;
                new_col(curr_part(2,:)')=l-1;
                new_table(:, i) = new_col;
                new_shells{end+1}=new_table;
            end
            four_col_flag=true;
        elseif sum(curr_col==0)==0
            new_shells=partial_shells;
        elseif sum(curr_col==0)==2
            new_col=curr_col;
            new_table=curr_table;
            zero_indices = find(new_col == 0);
            new_col(zero_indices)=l;
            new_table(:, i) = new_col;
            new_shells{end+1}=new_table;
            two_col_flag=true;
        elseif sum(curr_col == 0)==1 || sum(curr_col == 0)==3 || sum(curr_col == 0)==5 || sum(curr_col == 0)==6
            'EXCEPTION: incorrect number of pairings'
        end
    end

    if four_col_flag==true
        l=l-2;
    elseif two_col_flag==true
        l=l-1;
    end

    partial_shells=new_shells;
end

end

% given a col with four zeros, generate the three possible parition of the indices of these four rows.
function partitions = four_rows_partition(col)
zero_indices = find(col == 0);

% Step 2: Generate all 2-element subsets of these indices
combinations = nchoosek(zero_indices, 2);

% Step 3: Form partitions by pairing the subsets
num_combinations = size(combinations, 1);
partitions = cell(num_combinations, 1);

k = 1;
for i = 1:num_combinations-1
    for j = i+1:num_combinations
        % Ensure no overlap between pairs
        if isempty(intersect(combinations(i,:), combinations(j,:)))
            partitions{k} = [combinations(i,:); combinations(j,:)];
            k = k + 1;
        end
    end
end
partitions = partitions(~cellfun('isempty', partitions)); % Remove empty cells
% disp('Possible partitions:');
% for i = 1:length(partitions)
%     disp(partitions{i});
% end
end

function input_info = getInputInfo(M)

input_info{1} =[]; %12
input_info{2} = []; %13
input_info{3} = []; %14
input_info{4} = []; %15
input_info{5} = []; %16
input_info{6} = []; %23
input_info{7} = []; %24
input_info{8} = []; %25
input_info{9} = []; %26
input_info{10} = []; %34
input_info{11} = []; %35
input_info{12} = []; %36
input_info{13} = []; %45
input_info{14} = []; %46
input_info{15} = []; %56
input_info{16}=size(M,2); % total columns
input_info{17}=0; % total elements
input_info{18}=1; % multiplicative factor
input_info{19}=1; % pairing number
input_info{20}=1; % shell sign
input_info{21}=1; % operation sign

rows_info=[0,1,2,3,4,5;1,0,6,7,8,9;2,6,0,10,11,12;3,7,10,0,13,14;4,8,11,13,0,15;5,9,12,14,15,0];

[shell_info, shell_table] = getShellInfo(M);

input_info{21}=compute_shell_sign(shell_info,shell_table);

unique_elements = unique(M);
% Count the number of unique elements
input_info{17} = numel(unique_elements);

min_element=min(M(:));
max_element=max(M(:));

for i=min_element:-1
    logical_matrix = (M == i);

    % Step 2: Find the row indices where i appears
    [row_indices, ~] = find(logical_matrix);

    % Get unique row indices
    unique_row_indices = unique(row_indices);

    if size(unique_row_indices,1)==2
        upd=rows_info(unique_row_indices(1),unique_row_indices(2));

        input_info{upd}=[input_info{upd};0,-1];
    else
        'EXCEPTION: info table 1'
    end
end

for i=1:max_element
    logical_matrix = (M == i);
    % Step 2: Find the row indices where i appears
    [row_indices, ~] = find(logical_matrix);
    % Get unique row indices
    unique_row_indices = unique(row_indices);

    if size(unique_row_indices,1)==2
        upd=rows_info(unique_row_indices(1),unique_row_indices(2));
        input_info{upd}=[input_info{upd};i,-1];

    elseif size(unique_row_indices,1)==4
        complement_rows=setdiff([1,2,3,4,5,6],unique_row_indices);
        upd=rows_info(complement_rows(1),complement_rows(2));
        input_info{upd}=[input_info{upd};i,0];
    elseif size(unique_row_indices,1)==6
    else
        'EXCEPTION: info table 2'
    end
end

% merge + - of same label 
for i =1:15
    curr_M=input_info{i};
    if isempty(curr_M)
        continue
    end
    firstColumn = curr_M(:,1); % Extract the first column
    distinctElements = unique(firstColumn);
    distinctElements=setdiff(distinctElements,0);
    for j=1:length(distinctElements)
        a=distinctElements(j);
        rowsWithA = curr_M(:,1) == a;
        % Extract the second column elements from these rows
        outputElements = curr_M(rowsWithA, 2);
        if length(unique(outputElements))>=2
            curr_M(rowsWithA, :) = [];
            input_info{i}=curr_M;
            if isempty(curr_M)
                input_info{i}=[];
            end
        end
    end
end


end


function [shell_info, shell_table] = getShellInfo(M)

shell_info{1} =[]; %12
shell_info{2} = []; %13
shell_info{3} = []; %14
shell_info{4} = []; %15
shell_info{5} = []; %16
shell_info{6} = []; %23
shell_info{7} = []; %24
shell_info{8} = []; %25
shell_info{9} = []; %26
shell_info{10} = []; %34
shell_info{11} = []; %35
shell_info{12} = []; %36
shell_info{13} = []; %45
shell_info{14} = []; %46
shell_info{15} = []; %56
rows_info=[0,1,2,3,4,5;1,0,6,7,8,9;2,6,0,10,11,12;3,7,10,0,13,14;4,8,11,13,0,15;5,9,12,14,15,0];
% to positive
M(M > 0) = M(M > 0) - 1;
min_element=min(M(:));
M=M-((min_element-1)*ones(size(M)));
max_element=max(M(:));

shell_table=M;

for i=1:max_element
    logical_matrix = (M == i);

    % Step 2: Find the row indices where i appears
    [row_indices, ~] = find(logical_matrix);

    % Get unique row indices
    unique_row_indices = unique(row_indices);

    if size(unique_row_indices,1)==2
        upd=rows_info(unique_row_indices(1),unique_row_indices(2));
        shell_info{upd}=[shell_info{upd};i,-1];
    elseif size(unique_row_indices,1)==4
        complement_rows=setdiff([1,2,3,4,5,6],unique_row_indices);
        upd=rows_info(complement_rows(1),complement_rows(2));
        shell_info{upd}=[shell_info{upd};i,0];
    elseif size(unique_row_indices,1)==6
    else
        'EXCEPTION: shell table'
    end
end


end

function count = countColumnsWithOddZeros(matrix)
% countColumnsWithFiveZeros counts the number of columns with exactly five zeros
% matrix: The input matrix

% Initialize the count
count = 0;

% Get the number of columns
[~, numCols] = size(matrix);

% Loop through each column
for col = 1:numCols
    % Count the number of zeros in the current column
    numZeros = sum(matrix(:, col) == 0);

    % Check if the current column has exactly five zeros
    if mod(numZeros,2)==1
        count = count + 1;
    end
end
end

% given a second-stage table (not choosing the pairing of zeros), generate
% the generating function of this table
function [latex_string,GF_key]=latex_output(M)
s='';
col_list=zeros(1,size(M,2));

for i=1:size(M,2)
    curr_col=M(:,i);
    [col_latex,k]=column_latex_output(curr_col);
    s=append(s,col_latex);
    col_list(i)=k;
end
latex_string=s;
col_list=sort(col_list);
GF_key=mat2str(col_list);
GF_key=GF_key(~isspace(GF_key));
end

% the input is a column, we first convert it into the canonical form
% up to three marks

% TODO: we may need a double check for the generating function of four
% columns
function [latex_string,key]=column_latex_output(col)
col=column_canonical(col);

s='';
key=0;
for i=1:6
    if col(i)<0
        s=append(s,'x');
    elseif col(i)==0
        s=append(s,'0');
    else
        s=append(s,char(96+col(i)));
    end
end

% switch based on s
switch s
    case 'aaabbb'
        latex_string='((Subscript[\[Mu], 3]^2 t) / (1 + Subscript[\[Mu], 3]^2 t))';
        key=1;
    case '00aaaa'
        latex_string='((Subscript[\[Mu], 4] - 3) t / (1 - (Subscript[\[Mu], 4] - 3) t))';
        key=2;
    case 'x00aaa'
        latex_string='(Subscript[\[Mu], 3] t / (1 + Subscript[\[Mu], 3]^2 t) / (1 - (Subscript[\[Mu], 4] - 3) t))';
        key=3;
    case 'xx0000'
        latex_string='(t(1/(1-(Subscript[\[Mu],4]-3)t)-2Subscript[\[Mu], 3]^2 t / (1 + Subscript[\[Mu], 3]^2 t)) / (1 - (Subscript[\[Mu], 4] - 3) t)^2)';
        key=4;
    case 'xxaaaa'
        latex_string='(((Subscript[\[Mu], 4] - 3) t) / (1 - (Subscript[\[Mu], 4] - 3) t))';
        key=5;
    case 'xxxaaa'
        latex_string='((Subscript[\[Mu],3]t/(1+Subscript[\[Mu],3]^2t))*((1+2(Subscript[\[Mu],4]-3)t)/(1 -(Subscript[\[Mu],4]-3)t)))';
        key=6;
    case 'xxxx00'
        latex_string='t(1+4(-Subscript[\[Mu],3]^2 t)/(1+Subscript[\[Mu], 3]^2 t))(1+7(Subscript[\[Mu],4]-3)t/(1-(Subscript[\[Mu],4]-3)t))';
        key=6;
    otherwise
        disp(s)
end


end

function s=column_canonical(col)
col=sort(-col);
% Step 1: Replace all negative elements with -1
col(col < 0) = -1;

% Step 2: Identify unique positive elements in descending order
positive_elements = unique(col(col > 0), 'stable');
ranks = 1:length(positive_elements);

% Step 3: Substitute the positive elements with their corresponding ranks
for i = 1:length(positive_elements)
    col(col == positive_elements(i)) = ranks(i);
end
s=col;
end

function s=cell_to_string(C)
s='';
for i=1:15
    s=append(s,sortrows(mat2str(C{i})));
    s=append(s,'|');
end
for i=16:21
    s=append(s,num2str(C{i}));
    s=append(s,'|');
end
end