
function [f]=m2v(F,sp)
   
    f=zeros(sp*sp+1,1);

    %Build vector
    for i=1:sp
        for j=1:sp
           f((i-1)*sp+j,1)=F(i,j);
        end
    end
    f(sp*sp+1,1)=F(i+1,j+1);
end