function banthia(d,React,ste_p,tp,K_beam,W_beam,Sup_energy,D_energy)

    %INTEGER
    NODO=2110;
    sp=2;

    Frac(ste_p,1)=0;
    Frac2(ste_p,1)=0;
    for i=2:ste_p

        i_d=-(d(NODO*sp,i)-d(NODO*sp,i-1));
        F_m=(React(i)+React(i-1))/2;
        Frac(i)=Frac(i-1)+F_m*i_d/1000;
        Frac2(i)=Frac(i)-K_beam(i)/1000+W_beam(i)/1000;

    end

    figure;
    plot(tp,Frac,tp,K_beam/1000,tp,W_beam/1000,tp,-Sup_energy/1000, tp, D_energy/1000);
    figure;
    plot(tp,Frac,tp,K_beam/1000-Sup_energy/1000+D_energy/1000+W_beam/1000,...
    tp,K_beam/1000-Sup_energy/1000+W_beam/1000);


end

