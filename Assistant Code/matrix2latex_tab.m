function [str]=matrix2latex_tab(A,pre)
for i =1:size(A,1)
    str{i} = '';
    for j=1:size(A,2)
        str{i} = [str{i} ' & {' num2str(A(i,j),pre) '}'];
    end
    display(str{i});
end