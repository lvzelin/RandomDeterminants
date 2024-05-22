function p=compute_shell_sign(shell_info,shell)
T=createCanonicalPositive(shell_info);
if isempty(T)
    p=+1;
else
    p=sign_of_shell(T,shell);
end
end

function T=createCanonicalPositive(info)
% first compute the overall size
total_elements=0;
for i=1:15
    A=info{i};
    num_neg = sum(A(:) == -1);
    total_elements=total_elements+num_neg*2;
    num_pos = sum(A(:) == 0);
    total_elements=total_elements+num_pos*4;
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
        if A(j,2)==-1
            cur_rows=indices_info(i,:);
            cur_col=cur_index(cur_rows(1));
            cur_index(cur_rows(1))=cur_index(cur_rows(1))+1;
            T(cur_rows(1),cur_col)=A(j,1);
            second_col=cur_index(cur_rows(2));
            cur_index(cur_rows(2))=cur_index(cur_rows(2))+1;
            T(cur_rows(2),second_col)=A(j,1);
        end

        if A(j,2)==0
            cur_rows=setdiff(nsix,indices_info(i,:));

            for k=1:4
                cur_col=cur_index(cur_rows(k));
                cur_index(cur_rows(k))=cur_index(cur_rows(k))+1;
                T(cur_rows(k),cur_col)=A(j,1);
            end

        end
    end
end

end

% T is reference, S is the shell
function sgn = sign_of_shell(T,S)
sgn=1;
for i=1:6
    row_sign = sign_of_permutation(T(i, :))*sign_of_permutation(S(i, :));
    sgn = row_sign * sgn;
end
end

function sgn = sign_of_permutation(p)
sgn = 1;  % Initialize sign
n = length(p);
for i = 1:n
    for j = i+1:n
        sgn = sgn * sign(p(i) - p(j));
    end
end
end