

function [F]=v2m(def_G,sp)
   
    F=zeros(sp+1,sp+1);

    %Build matrix
    for i=1:sp
        for j=1:sp
            F(i,j)=def_G((i-1)*sp+j,1);
        end
    end
    F(i+1,j+1)=def_G(sp*sp+1,1);
end