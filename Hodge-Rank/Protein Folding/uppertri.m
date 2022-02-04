%This function takes the upper triangular values of matrix M such that
%there is an edge between ith and jth entry.
function V = uppertri(M,W)
    s = size(M,1);
    V = [];
    for i = 1:(s-1)
        for j = (i+1):(s)
            if W(i,j)>0
                V = [V;M(i,j)];
            end
        end
    end
end