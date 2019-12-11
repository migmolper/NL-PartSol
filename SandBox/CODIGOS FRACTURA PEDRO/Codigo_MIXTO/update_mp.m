
function [xg]=update_mp(d,near,p,xg)
    
    [elements,sp]=size(xg);
        
    %% Update **********************
    for i=1:elements
        m=length(near{i});
        nd = near{i};
        sh= p{i};
        for t1=1:m
            for k=1:sp
                xg(i,k)=xg(i,k)...
                +sh(t1)*(d((nd(t1)-1)*sp+k,1)-d((nd(t1)-1)*sp+k,2));
            end
        end
    end  
end