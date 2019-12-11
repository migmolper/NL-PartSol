
function [B]=b_m(elements,near,dp,sp,wrap,B)

    for i=1:elements
        if wrap(i)==1
            nb=near{i};
            sh=dp{i};
            n=length(nb);
            b=zeros(3,sp*n);
            for j=1:n
                b(1,sp*j-1)=sh(j,1);
                b(2,sp*j)  =sh(j,2);
                b(3,sp*j-1)=sh(j,2);
                b(3,sp*j)  =sh(j,1);
            end
            B(i)={b};
            clear t b;
        end
    end

end