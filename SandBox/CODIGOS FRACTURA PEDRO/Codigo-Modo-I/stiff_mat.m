
function [stiff_mtx]=stiff_mat(B,near,volume,i,stiff_mtx,T,A)

    global sp
    % ----------------------------
    % Geometric stiffness
    % ----------------------------
    nb=near{i};
    n=length(nb);
    sh=zeros(n*2,4);
    b=B{i};
    for j=1:n
        sh(j*2-1,1)=b(1,j*2-1);
        sh(j*2-1,2)=b(2,j*2);
        sh(j*2,3)=b(1,j*2-1);
        sh(j*2,4)=b(2,j*2);
    end

    Sigma=zeros(sp*2);
    Sigma(1:2,1:2)=T(1:2,1:2);
    Sigma(3:4,3:4)=T(1:2,1:2);
    
    K_geo=sh*Sigma*sh';
    
    % ----------------------------
    % Assembling of K
    % ----------------------------
    D=A(1:3,1:3);
    K_el=volume(i)*(b'*D*b+K_geo);

    for j=1:n
        for l=1:n
            for k=1:sp
                for r=1:sp
                    stiff_mtx(nb(j)*sp+1-k,nb(l)*sp+1-r)...
                    =stiff_mtx(nb(j)*sp+1-k,nb(l)*sp+1-r)...
                    +K_el(j*sp+1-k,l*sp+1-r);
                end
            end
        end
    end


end