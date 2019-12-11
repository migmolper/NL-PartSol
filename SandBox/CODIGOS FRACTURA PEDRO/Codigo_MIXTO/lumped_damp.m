function [C]=lumped_damp(near,p,volume,perm)

global elements nodes sp
C=zeros(nodes*sp,nodes*sp);
    %% Lumped Mass **********************
    for i=1:elements
        m=length(near{i});
        nd = near{i};
        sh= p{i};
        for t1=1:m
            for k=1:sp
                C(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                C(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                    +volume(i)*sh(t1)/perm(i);
            end
        end
    end
end