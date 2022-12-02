function DispSymMatrixMatlab_file(AA,NAME,nameOfFunction)

if ndims(AA) == 2
    N1 = size(AA,1);
    N2 = size(AA,2);
    N3 = 1;
    N4 = 1;
elseif ndims(AA) == 3
    N1 = size(AA,1);
    N2 = size(AA,2);
    N3 = size(AA,3);
    N4 = 1;
elseif ndims(AA) == 4
    N1 = size(AA,1);
    N2 = size(AA,2);
    N3 = size(AA,3);
    N4 = size(AA,4);
end
text = ['% START DISPLAYING ' NAME];
dlmwrite(nameOfFunction,text,'-append','delimiter','');
for g=1:N4
    for k=1:N3
        for j=1:N2
            if any(AA(:,j,k,g) ~= 0)
                for i=1:N1
                    if AA(i,j,k,g) ~= 0
                        if ndims(AA) == 2
                            text = [NAME '(' num2str(i) ',' num2str(j) ')=' char(AA(i,j,k,g)) ';'];
                        elseif ndims(AA) == 3
                            text = [NAME '(' num2str(i) ',' num2str(j) ',' num2str(k) ')=' char(AA(i,j,k,g)) ';'];
                        elseif ndims(AA) == 4
                            text = [NAME '(' num2str(i) ',' num2str(j) ',' num2str(k) ',' num2str(g) ')=' char(AA(i,j,k,g)) ';'];
                        end
                        dlmwrite(nameOfFunction,text,'-append','delimiter','');
                    end
                end
            end
        end
    end
end
text = ['% END DISPLAYING ' NAME];
dlmwrite(nameOfFunction,text,'-append','delimiter','');
end
