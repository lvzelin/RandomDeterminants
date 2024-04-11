% the first input is used to indicate which element it is (but when merging, this is not implemented),
% the second input is used to indicate the sign, -1 for - and 0 for +
% e.g. [1,0;0,-1] means + with label 1 and - with no label
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
input_info{15} = [0,-1;0,-1;2,0;1,0]; %56
input_info{16}=2; % total columns
input_info{17}=4; % total elements
input_info{18}=1; % multiplicative factor
input_info{19}=1; % pairing number
filePath = './results/two_mark_c_11.txt';
filePath = './results/randomrandom.txt';

indices_info=[1,2;1,3;1,4;1,5;1,6;2,3;2,4;2,5;2,6;3,4;3,5;3,6;4,5;4,6;5,6];
rows_info=[0,1,2,3,4,5;1,0,6,7,8,9;2,6,0,10,11,12;3,7,10,0,13,14;4,8,11,13,0,15;5,9,12,14,15,0];

all_possible_table=cell(1);
% initilization
all_possible_table{1}=input_info;

% take care of + - of same label at first
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
all_possible_table{1}=input_info;

% minus_plus_case
for i=1:15
    if isempty(input_info{i})
        continue
    end
    merged_flag=false;
    curr_M=input_info{i};
    curr_signs=curr_M(:,2);
    num_minus=sum(curr_signs==-1);
    num_plus=sum(curr_signs==0);

    if min(num_minus,num_plus)<1
        continue
    end

    col2IsZero = curr_M(:, 2) == 0;
    col1IsPositive = curr_M(:, 1) > 0;
    rowsMeetConditions = col2IsZero & col1IsPositive;
    num_labeled_plus = sum(rowsMeetConditions);

    col2IsZero = curr_M(:, 2) == -1;
    col1IsPositive = curr_M(:, 1) > 0;
    rowsMeetConditions = col2IsZero & col1IsPositive;
    num_labeled_minus = sum(rowsMeetConditions);

    num_unlabeled_plus=num_plus-num_labeled_plus;
    num_unlabeled_minus=num_minus-num_labeled_minus;

    new_possible_table=cell(1);
    k=1;

    for j=1:min(num_unlabeled_minus,num_plus)
        merged_flag=true;
        mf=(factorial(num_unlabeled_minus)/factorial(num_unlabeled_minus-j))*(factorial(num_plus)/factorial(num_plus-j))/(factorial(j));
        % mf=mf*nchoosek(num_plus,j);
        for l=1:size(all_possible_table,2)
            curr_input_info=all_possible_table{l};
            curr_input_info{17}=curr_input_info{17}-j;
            curr_input_info{18}=curr_input_info{18}*mf;
            % merge j of + -
            curr_input_info{i}=mergeMinusPlus(curr_input_info{i},j);
            if (isempty(curr_input_info{i}))
                curr_input_info{i}=[];
            end
            new_possible_table{k}=curr_input_info;
            k=k+1;
        end
    end

    if merged_flag
        all_possible_table=[all_possible_table,new_possible_table];
    end
    
end

% 3 minuses
indices_to_row=[1,2;1,3;1,4;1,5;1,6;2,3;2,4;2,5;2,6;3,4;3,5;3,6;4,5;4,6;5,6];
table_index=1;
while table_index<=size(all_possible_table,2)
    curr_input_info=all_possible_table{table_index};
    for i=1:15
        for j=i:15
            for k=j:15
                if isempty(intersect(indices_to_row(i,:),indices_to_row(j,:))) & isempty(intersect(indices_to_row(k,:),indices_to_row(j,:)))&isempty(intersect(indices_to_row(k,:),indices_to_row(j,:)))
                    i_matrix=curr_input_info{i};
                    j_matrix=curr_input_info{j};
                    k_matrix=curr_input_info{k};
                    if isempty(i_matrix) | isempty(j_matrix) | isempty(k_matrix)
                        continue
                    end
                    % above are basic check, now check if they contain
                    % minus signs and at most one of these minus are
                    % labeled

                    num_i=sum(i_matrix(:,2) == -1);
                    num_j=sum(j_matrix(:,2) == -1);
                    num_k=sum(k_matrix(:,2) == -1);

                    col2IsZero = i_matrix(:, 2) == -1;
                    col1IsPositive = i_matrix(:, 1) > 0;
                    rowsMeetConditions = col2IsZero & col1IsPositive;
                    i_labeled_minus = sum(rowsMeetConditions);

                    col2IsZero = j_matrix(:, 2) == -1;
                    col1IsPositive = j_matrix(:, 1) > 0;
                    rowsMeetConditions = col2IsZero & col1IsPositive;
                    j_labeled_minus = sum(rowsMeetConditions);

                    col2IsZero = k_matrix(:, 2) == -1;
                    col1IsPositive = k_matrix(:, 1) > 0;
                    rowsMeetConditions = col2IsZero & col1IsPositive;
                    k_labeled_minus = sum(rowsMeetConditions);

                    i_unlabeled=num_i-i_labeled_minus;
                    j_unlabeled=num_j-j_labeled_minus;
                    k_unlabeled=num_k-k_labeled_minus;
                    new_possible_table=cell(1);
                    merged_flag=false;
                    r=1;
                    % all three to be unlabeled
                    if min([i_unlabeled,j_unlabeled,k_unlabeled])>=1
                        merged_flag=true;
                        % loop through the number of interactions across i j k
                        for l=1:min([i_unlabeled,j_unlabeled,k_unlabeled])
                            mf=(factorial(i_unlabeled)/factorial(i_unlabeled-l))*(factorial(j_unlabeled)/factorial(j_unlabeled-l))*(factorial(k_unlabeled)/factorial(k_unlabeled-l))/(factorial(l));
                            % 0 to indicate only merge unlabeled minus
                            [new_i,new_j,new_k]=mergeTripleMinus(i_matrix,j_matrix,k_matrix,l,0);
                            new_input_info=curr_input_info;

                            new_input_info{17}=curr_input_info{17}-2*l;
                            new_input_info{18}=curr_input_info{18}*mf;
                            new_input_info{i}=new_i;
                            new_input_info{j}=new_j;
                            new_input_info{k}=new_k;
                            if (isempty(new_input_info{i}))
                                new_input_info{i}=[];
                            end
                            if (isempty(new_input_info{j}))
                                new_input_info{j}=[];
                            end
                            if (isempty(new_input_info{k}))
                                new_input_info{k}=[];
                            end
                            new_possible_table{r}=new_input_info;
                            r=r+1;
                        end
                    end

                    % i is labeled and j,k are not labeled
                    if min([i_labeled_minus,j_unlabeled,k_unlabeled])>=1
                        merged_flag=true;
                        % loop through the number of interactions across i j k
                        for l=1:min([i_labeled_minus,j_unlabeled,k_unlabeled])
                            mf=(factorial(i_labeled_minus)/factorial(i_labeled_minus-l))*(factorial(j_unlabeled)/factorial(j_unlabeled-l))*(factorial(k_unlabeled)/factorial(k_unlabeled-l))/(factorial(l));
                            % 0 to indicate only merge unlabeled minus
                            [new_i,new_j,new_k]=mergeTripleMinus(i_matrix,j_matrix,k_matrix,l,1);
                            new_input_info=curr_input_info;

                            new_input_info{17}=curr_input_info{17}-2*l;
                            new_input_info{18}=curr_input_info{18}*mf;
                            new_input_info{i}=new_i;
                            new_input_info{j}=new_j;
                            new_input_info{k}=new_k;
                            if (isempty(new_input_info{i}))
                                new_input_info{i}=[];
                            end
                            if (isempty(new_input_info{j}))
                                new_input_info{j}=[];
                            end
                            if (isempty(new_input_info{k}))
                                new_input_info{k}=[];
                            end
                            new_possible_table{r}=new_input_info;

                            r=r+1;
                        end
                    end

                    % j is labeled and i,k are not labeled
                    if min([i_unlabeled,j_labeled_minus,k_unlabeled])>=1
                        merged_flag=true;
                        % loop through the number of interactions across i j k
                        for l=1:min([i_unlabeled,j_labeled_minus,k_unlabeled])
                            mf=(factorial(i_unlabeled)/factorial(i_unlabeled-l))*(factorial(j_labeled_minus)/factorial(j_labeled_minus-l))*(factorial(k_unlabeled)/factorial(k_unlabeled-l))/(factorial(l));
                            % 0 to indicate only merge unlabeled minus
                            [new_i,new_j,new_k]=mergeTripleMinus(i_matrix,j_matrix,k_matrix,l,2);
                            new_input_info=curr_input_info;

                            new_input_info{17}=curr_input_info{17}-2*l;
                            new_input_info{18}=curr_input_info{18}*mf;
                            new_input_info{i}=new_i;
                            new_input_info{j}=new_j;
                            new_input_info{k}=new_k;
                            if (isempty(new_input_info{i}))
                                new_input_info{i}=[];
                            end
                            if (isempty(new_input_info{j}))
                                new_input_info{j}=[];
                            end
                            if (isempty(new_input_info{k}))
                                new_input_info{k}=[];
                            end
                            new_possible_table{r}=new_input_info;
                            r=r+1;
                        end
                    end

                    % k is labeled and i,j are not labeled
                    if min([i_unlabeled,j_unlabeled,k_labeled_minus])>=1
                        merged_flag=true;
                        for l=1:min([i_unlabeled,j_unlabeled,k_labeled_minus])
                            mf=(factorial(i_unlabeled)/factorial(i_unlabeled-l))*(factorial(j_unlabeled)/factorial(j_unlabeled-l))*(factorial(k_labeled_minus)/factorial(k_labeled_minus-l))/(factorial(l));
                            % 0 to indicate only merge unlabeled minus
                            [new_i,new_j,new_k]=mergeTripleMinus(i_matrix,j_matrix,k_matrix,l,3);
                            new_input_info=curr_input_info;

                            new_input_info{17}=curr_input_info{17}-2*l;
                            new_input_info{18}=curr_input_info{18}*mf;
                            new_input_info{i}=new_i;
                            new_input_info{j}=new_j;
                            new_input_info{k}=new_k;
                            if (isempty(new_input_info{i}))
                                new_input_info{i}=[];
                            end
                            if (isempty(new_input_info{j}))
                                new_input_info{j}=[];
                            end
                            if (isempty(new_input_info{k}))
                                new_input_info{k}=[];
                            end
                            new_possible_table{r}=new_input_info;
                            r=r+1;
                        end
                    end
                    if merged_flag
                        all_possible_table=[all_possible_table,new_possible_table];
                    end
                end
            end
        end
    end
    table_index=table_index+1;
end

% 2 minus
table_index=1;
while table_index<=size(all_possible_table,2)

    curr_input_info=all_possible_table{table_index};
    for i=1:15
        for j=i:15
            if isempty(intersect(indices_to_row(i,:),indices_to_row(j,:)))
                    i_matrix=curr_input_info{i};
                    j_matrix=curr_input_info{j};
                    ij_union=union(indices_info(i,:),indices_info(j,:));
                    k_index=setdiff([1,2,3,4,5,6],ij_union);
                    k=rows_info(k_index(1),k_index(2));
                    k_matrix=curr_input_info{k};

                    if isempty(i_matrix) | isempty(j_matrix)
                        continue
                    end

                    % above are basic check, now check if they contain
                    % minus signs and at most one of these minus are
                    % labeled

                    num_i=sum(i_matrix(:,2) == -1);
                    num_j=sum(j_matrix(:,2) == -1);

                    col2IsZero = i_matrix(:, 2) == -1;
                    col1IsPositive = i_matrix(:, 1) > 0;
                    rowsMeetConditions = col2IsZero & col1IsPositive;
                    i_labeled_minus = sum(rowsMeetConditions);

                    col2IsZero = j_matrix(:, 2) == -1;
                    col1IsPositive = j_matrix(:, 1) > 0;
                    rowsMeetConditions = col2IsZero & col1IsPositive;
                    j_labeled_minus = sum(rowsMeetConditions);

                    i_unlabeled=num_i-i_labeled_minus;
                    j_unlabeled=num_j-j_labeled_minus;

                    new_possible_table=cell(1);
                    merged_flag=false;
                    r=1;
                    % all two are unlabeled
                    if min([i_unlabeled,j_unlabeled])>=1
                        merged_flag=true;
                        % loop through the number of interactions across i j k
                        for l=1:min([i_unlabeled,j_unlabeled])
                            mf=(factorial(i_unlabeled)/factorial(i_unlabeled-l))*(factorial(j_unlabeled)/factorial(j_unlabeled-l))/(factorial(l));
                            % 0 to indicate only merge unlabeled minus
                            [new_i,new_j,new_k]=mergeDoubleMinus(i_matrix,j_matrix,k_matrix,l,0);
                            new_input_info=curr_input_info;
                            new_input_info{17}=curr_input_info{17}-l;
                            new_input_info{18}=curr_input_info{18}*mf;
                            new_input_info{i}=new_i;
                            new_input_info{j}=new_j;
                            new_input_info{k}=new_k;

                            if (isempty(new_input_info{i}))
                                new_input_info{i}=[];
                            end
                            if (isempty(new_input_info{j}))
                                new_input_info{j}=[];
                            end
                            if (isempty(new_input_info{k}))
                                new_input_info{k}=[];
                            
                            end
                            new_possible_table{r}=new_input_info;
                            r=r+1;
                        end
                    end
                    
                    % i is labeled and j is not labeled
                   if min([i_labeled_minus,j_unlabeled])>=1
                        merged_flag=true;
                        % loop through the number of interactions across i j k
                        for l=1:min([i_labeled_minus,j_unlabeled])
                            mf=(factorial(i_labeled_minus)/factorial(i_labeled_minus-l))*(factorial(j_unlabeled)/factorial(j_unlabeled-l))/(factorial(l));
                            % 0 to indicate only merge unlabeled minus
                            [new_i,new_j,new_k]=mergeDoubleMinus(i_matrix,j_matrix,k_matrix,l,1);
                            new_input_info=curr_input_info;

                            new_input_info{17}=curr_input_info{17}-l;
                            new_input_info{18}=curr_input_info{18}*mf;
                            new_input_info{i}=new_i;
                            new_input_info{j}=new_j;
                            new_input_info{k}=new_k;

                            if (isempty(new_input_info{i}))
                                new_input_info{i}=[];
                            end
                            if (isempty(new_input_info{j}))
                                new_input_info{j}=[];
                            end
                            if (isempty(new_input_info{k}))
                                new_input_info{k}=[];
                            
                            end
                            new_possible_table{r}=new_input_info;
                            r=r+1;
                        end
                    end
                    % j is labeled and i is not labeled
                    % all two are unlabeled
                    if min([i_unlabeled,j_labeled_minus])>=1
                        merged_flag=true;
                        % loop through the number of interactions across i j k
                        for l=1:min([i_unlabeled,j_labeled_minus])
                            mf=(factorial(i_unlabeled)/factorial(i_unlabeled-l))*(factorial(j_labeled_minus)/factorial(j_labeled_minus-l))/(factorial(l));
                            % 0 to indicate only merge unlabeled minus
                            [new_i,new_j,new_k]=mergeDoubleMinus(i_matrix,j_matrix,k_matrix,l,2);
                            new_input_info=curr_input_info;

                            new_input_info{17}=curr_input_info{17}-l;
                            new_input_info{18}=curr_input_info{18}*mf;
                            new_input_info{i}=new_i;
                            new_input_info{j}=new_j;
                            new_input_info{k}=new_k;

                            if (isempty(new_input_info{i}))
                                new_input_info{i}=[];
                            end
                            if (isempty(new_input_info{j}))
                                new_input_info{j}=[];
                            end
                            if (isempty(new_input_info{k}))
                                new_input_info{k}=[];
                            
                            end
                            new_possible_table{r}=new_input_info;
                            r=r+1;
                        end
                    end
                    
                    if merged_flag
                        all_possible_table=[all_possible_table,new_possible_table];
                    end
            end
        end
    end
    table_index=table_index+1;
end

% use all possible table to obtain the generating function
for i=1:size(all_possible_table,2)
    curr_info=all_possible_table{i};
    curr_table=cell(1,15);
    output_table=cell(1,15);
    for j=1:15
        if isempty(curr_info{j})
            curr_table{j}=[];
            output_table{j}=[];
        else
            curr_table{j}=curr_info{j}(:,2)';
            output_table{j}=curr_info{j};
        end
    end

    if check_table_empty(curr_table)
        pairing_number=1;
    else
        pairing_number=compute_net_pairing_core(curr_table);
    end
    all_possible_table{i}{19}=pairing_number;

    s = formatMatrices(output_table);
    fileID = fopen(filePath, 'a');
    if fileID == -1
        error('Failed to open the file.');
    end
    
    fprintf(fileID, '%s\n', '\item');
    fprintf(fileID, '%s\n\n', 'table:');
    fprintf(fileID, '%s\n', s);
    k=curr_info{17}-curr_info{16};
    fprintf(fileID, 'value of k is: %d;\n', k);
    fprintf(fileID, 'net paring number is: %d;\n', pairing_number);
    denominator=1;
    for j=1:k
        denominator=j*(j+2)*(j+4)*denominator;
    end
    fprintf(fileID, 'denominator is: %d;\n \n', denominator);
    fprintf(fileID, 'factor for this term is: \n');
    fprintf(fileID, '$$\\frac{O(t)}{N(\\frac{t}{(1-(\\mu_4-3)t)^3})}(\\frac{t}{(1-(\\mu_4-3)t)^3})^{%d} \\frac{d^{%d}N}{dt^{%d}}(\\frac{t}{(1-(\\mu_4-3)t)^3})\\frac{%d}{%d}$$\n\n',k,k,k,pairing_number,denominator);
    % TO-DO: write the mathematica form
    % Close the file
    fclose(fileID);
end

% merge minus and plus in M for j times
function M_modified = mergeMinusPlus(M,j)
M_modified=M;
for i=1:j
    % Find rows where the second column is -1 and it's unlabeled
    rowsWithMinus = find(ismember(M_modified, [0, -1], 'rows'), 1);

    % Remove the first of those rows found
    if ~isempty(rowsWithMinus)
        M_modified(rowsWithMinus, :) = [];
    else
        % If no such row exists, return the original matrix
        'no minus to merge'
    end

    rowsWithPlus = find(M_modified(:,2) == 0, 1, 'first');
    % Remove the first of those rows found
    if ~isempty(rowsWithPlus)
        M_modified(rowsWithPlus, :) = [];
    end
end

end

function [i_M,j_M,k_M]=mergeTripleMinus(i_matrix,j_matrix,k_matrix,j,option_handle)
i_M=i_matrix;
j_M=j_matrix;
k_M=k_matrix;

for i=1:j
    % Find rows where the second column is -1 and unlabled
    i_rowsWithMinus = find(ismember(i_M, [0, -1], 'rows'), 1);
    j_rowsWithMinus = find(ismember(j_M, [0, -1], 'rows'), 1);
    k_rowsWithMinus = find(ismember(k_M, [0, -1], 'rows'), 1);

    if option_handle==1
        condition1 = i_M(:,1) > 0; % Finds rows where the first column is greater than 0
        condition2 = i_M(:,2) == -1; % Finds rows where the second column is -1
        combinedCondition = condition1 & condition2; % Combines the conditions using logical AND
        i_rowsWithMinus = find(combinedCondition, 1);
    elseif option_handle==2
        condition1 = j_M(:,1) > 0; % Finds rows where the first column is greater than 0
        condition2 = j_M(:,2) == -1; % Finds rows where the second column is -1
        combinedCondition = condition1 & condition2; % Combines the conditions using logical AND
        j_rowsWithMinus = find(combinedCondition, 1);
    elseif option_handle==3
        condition1 = k_M(:,1) > 0; % Finds rows where the first column is greater than 0
        condition2 = k_M(:,2) == -1; % Finds rows where the second column is -1
        combinedCondition = condition1 & condition2; % Combines the conditions using logical AND
        k_rowsWithMinus = find(combinedCondition, 1);
    end

    % Remove the first of those rows found
    if ~isempty(i_rowsWithMinus) & ~isempty(j_rowsWithMinus) & ~isempty(k_rowsWithMinus)
        i_M(i_rowsWithMinus, :) = [];
        j_M(j_rowsWithMinus, :) = [];
        k_M(k_rowsWithMinus, :) = [];
    else
        % If no such row exists, return the original matrix
        'no minus to merge'
    end
end

end

% merge 2 minus
function [i_M,j_M,k_M]=mergeDoubleMinus(i_matrix,j_matrix,k_matrix,j,option_handle)

i_M=i_matrix;
j_M=j_matrix;
k_M=k_matrix;

for i=1:j
    % Find rows where the second column is -1 and unlabled
    i_rowsWithMinus = find(ismember(i_M, [0, -1], 'rows'), 1);
    j_rowsWithMinus = find(ismember(j_M, [0, -1], 'rows'), 1);

    if option_handle==1
        condition1 = i_M(:,1) > 0; % Finds rows where the first column is greater than 0
        condition2 = i_M(:,2) == -1; % Finds rows where the second column is -1
        combinedCondition = condition1 & condition2; % Combines the conditions using logical AND
        i_rowsWithMinus = find(combinedCondition, 1);
    elseif option_handle==2
        condition1 = j_M(:,1) > 0; % Finds rows where the first column is greater than 0
        condition2 = j_M(:,2) == -1; % Finds rows where the second column is -1
        combinedCondition = condition1 & condition2; % Combines the conditions using logical AND
        j_rowsWithMinus = find(combinedCondition, 1);
    end

    % Remove the first of those rows found
    if ~isempty(i_rowsWithMinus) & ~isempty(j_rowsWithMinus)
        new_k_row_label=max(i_M(i_rowsWithMinus, 1),j_M(j_rowsWithMinus, 1));
        new_k_row=[new_k_row_label,0];
        k_M=[k_M;new_k_row];
        i_M(i_rowsWithMinus, :) = [];
        j_M(j_rowsWithMinus, :) = [];
        
    else
        % If no such row exists, return the original matrix
        'no minus to merge'
    end
end
end

function flag=check_table_empty(C)
flag=true;
for i =1:15
    if ~isempty(C{i})
        flag=false;
    end
end
end

function M = cell_to_matrix(C)
% Check if the input cell array is empty
if isempty(C)
    M = [];
    return;
end
% Preallocate a matrix to hold the combined data
% First, calculate the total number of rows required
totalRows = sum(cellfun(@(x) size(x, 1), C));
% Determine the number of columns based on the first element of the cell array
numCols = size(C{1}, 2);
% Preallocate the matrix with zeros
M = zeros(totalRows, numCols);
% Now, concatenate all matrices vertically
lastRow = 0;
for i = 1:length(C)
    currentMatrix = C{i};
    numRows = size(currentMatrix, 1);

    % Assign the matrix into the preallocated space
    M(lastRow+1:lastRow+numRows, :) = currentMatrix;

    % Update the last row index
    lastRow = lastRow + numRows;
end
end

function outputTable = formatMatrices(input_info)
    % Initialize the output table with the header row
    indices = ["12", "13", "14", "15", "16", "23", "24", "25", "26", "34", "35", "36", "45", "46", "56"];
    headerRow = strjoin(indices, " & ");
    
    % Start building the output as a string, including LaTeX table structure
    outputTable = "\begin{table}[H] \begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|l|l|l|l|} \hline ";
    outputTable = outputTable + headerRow + " \\ \hline ";
    
    % Determine the maximum number of rows among all matrices for row iteration
    maxRows = max(cellfun(@(c) size(c, 1), input_info));

    % Iterate through each possible row
    for rowIdx = 1:maxRows
        rowContent = strings(1, length(input_info)); % Initialize row content
        
        for matrixIdx = 1:length(input_info)
            matrix = input_info{matrixIdx};
            
            if rowIdx <= size(matrix, 1) % Check if the current matrix has this row
                % Convert first element to character
                charPart = char(matrix(rowIdx, 1) + 96);
                
                % Determine the sign based on the second element
                if matrix(rowIdx, 2) == -1
                    signPart = '-';
                else
                    signPart = '+';
                end
                
                rowContent(matrixIdx) = strcat(charPart, signPart);
            else
                rowContent(matrixIdx) = " ";
            end
        end
        
        % Add row to outputTable if it's not just empty strings
        if any(rowContent ~= " ")
            outputRow = strjoin(rowContent, " & ");
            outputTable = outputTable + outputRow + " \\ \hline ";
        end
    end
    
    % Close the LaTeX table structure
    outputTable = outputTable + "\end{tabular} \end{table}";
end
