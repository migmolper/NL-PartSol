
function comp_sup_energy

load FILE
load X0

[nodes,sp]=size(x_0);
AMP=1;

Sup_energy(1)=0;

for ss=2:ste_p
    D=0;
    D1=0;
    for i=1:nodes
        if x_0(i,2)==min(x_0(:,2))
            if (x_0(i,1) <= 72.0*AMP && x_0(i,1) >= 48.0*AMP) || ...
            (x_0(i,1) <= 372.0*AMP && x_0(i,1) >= 348.0*AMP)
                D=min(D,d(sp*i,ss));
                D1=min(D1,d(sp*i,ss-1));
            end
        end
    end
    Sup_energy(ss)=Sup_energy(ss-1)+(React(ss)+React(ss-1))/2*(D-D1);
end

save FILE