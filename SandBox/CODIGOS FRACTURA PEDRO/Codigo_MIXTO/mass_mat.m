function [mass_mtx]=mass_mat(near,p,dens,volume)

    global elements nodes sp 
    mass_mtx=zeros(nodes*sp,nodes*sp);
    %% Mass matrix **********************
    for i=1:elements
        m=length(near{i});
        nd = near{i};
        sh= p{i};
        for t1=1:m
            for t2=1:m
                for k=1:sp
                    mass_mtx(nd(t1)*sp+1-k,nd(t2)*sp+1-k)=...
                        mass_mtx(nd(t1)*sp+1-k,nd(t2)*sp+1-k)+...
                        +dens(i)*volume(i)*sh(t1)*sh(t2);
                end
            end
        end
    end
end