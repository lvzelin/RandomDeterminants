% given a table of marks, generate all possible partial shells
% even some variables are named latex_something, it might not be in the
% latex formats
tic
table_input=[1,0,0,0,0,0;0,2,0,0,0,0;0,0,3,0,0,0]';
filePath = './results/three_mark_8.txt';

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
            for k=1:size(column_indices)
                curr_col=curr_table(:,k);
                zero_rows=any(curr_col == 0, 2);
                zero_row_indices = find(zero_rows);
                if isequal(target_rows{i},zero_row_indices')
                    new_table_2(target_rows{i}, k) = -i;
                    new_second_stage_shell{end+1}=new_table_2;
                end
            end
        end
        generated_second_stage_shell=[generated_second_stage_shell,new_second_stage_shell];
    end
end


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

% display
for j=1:size(generated_second_stage_shell,2)
    generated_second_stage_shell{j}
end
size(generated_second_stage_shell,2)

for j=1:size(generated_second_stage_shell,2)
    curr_second_stage_table=generated_second_stage_shell{j};
    mark_table_latex_string=latex_output(curr_second_stage_table);
    % now the labeled information can be deleted
    curr_second_stage_table=abs(curr_second_stage_table);
    % assign the pairing the generate partial shell, in the output negative
    % elements mean unlabeled
    partial_shells=generate_partial_shell(curr_second_stage_table);
    
    fileID = fopen(filePath, 'a');
    if fileID == -1
        error('Failed to open the file.');
    end

    % fprintf(fileID, '(*%d -th table*)%s(\n',j, mark_table_latex_string);
    fprintf(fileID, '%s(\n', mark_table_latex_string);

    % first we create the shell_info and shell_table for each partial shell
    % in the list of partial shells
    for i=1:size(partial_shells,2)
        input_info=getInputInfo(partial_shells{i});
        partial_shell_string=compute_factors_partial(input_info);
        fprintf(fileID, '+%s\n', partial_shell_string);
    end
    if j==size(generated_second_stage_shell,2)
        fprintf(fileID, ')');
        fclose(fileID);
        continue
    end
    fprintf(fileID, ')+\n');
    fclose(fileID);
end
toc

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

% for i=1:15
%     input_info{i}
% end
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
function latex_string=latex_output(M)
s='';
for i=1:size(M,2)
    curr_col=M(:,i);
    col_latex=column_latex_output(curr_col);
    s=append(s,col_latex);
end
latex_string=s;
end

% the input is a column, we first convert it into the canonical form
% up to three marks
function latex_string=column_latex_output(col)
col=column_canonical(col);

s='';
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
    case '00aaaa'
        latex_string='((Subscript[\[Mu], 4] - 3) t / (1 - (Subscript[\[Mu], 4] - 3) t))';
    case 'x00aaa'
        latex_string='(Subscript[\[Mu], 3] t / (1 + Subscript[\[Mu], 3]^2 t) / (1 - (Subscript[\[Mu], 4] - 3) t))';
    case 'xx0000'
        latex_string='(t (1 / (1 - (Subscript[\[Mu], 4] - 3) t) - 2 Subscript[\[Mu], 3]^2 t / (1 + Subscript[\[Mu], 3]^2 t)) / (1 - (Subscript[\[Mu], 4] - 3) t)^2)';
    case 'xxaaaa'
        latex_string='(((Subscript[\[Mu], 4] - 3) t) / (1 - (Subscript[\[Mu], 4] - 3) t))';
    case 'xxxaaa'
        latex_string='( (Subscript[\[Mu], 3] t / (1 + Subscript[\[Mu], 3]^2 t)) * ((1 + (2 Subscript[\[Mu], 4]- 3) t) / (1 - (Subscript[\[Mu], 4] - 3) t)))';
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