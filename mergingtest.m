% string_canonical_shell =[     
%      -1  ,  -2,3,4;
%      -1   , -2,3,4;
%     1   ,  -3,-5,4;
%     1   ,  -3,-5,4;
%     -4,2   , 3,-6;
%     -4,2   , 3,-6]

string_canonical_shell =[     
1,-1,-2;
1,-1,-7;
-3,-5,-2;
-3,-6,2;
-4,-5,2;
-4,-6,-7]

% string_canonical_shell =[     
% -1,-4,-6;
% -1,-4,-6;
% -2,1,-7;
% -2,1,-7;
% -3,-5,-8;
% -3,-5,-8]

input_info=getInputInfo(string_canonical_shell)
partial_shell_string=compute_factors_partial(input_info)

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

