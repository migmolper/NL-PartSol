function [lmass]=lumped_mass(near,p,volume,dens)

    global elements sp nodes
    lmass=zeros(nodes*sp,nodes*sp);
    
    %% Lumped Mass **********************
    for i=1:elements
        m=length(near{i});
        nd = near{i};
        sh= p{i};
        for t1=1:m
            for k=1:sp
                lmass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)=...
                lmass(nd(t1)*sp+1-k,nd(t1)*sp+1-k)...
                    +dens(i)*volume(i)*sh(t1);
            end
        end
    end
end