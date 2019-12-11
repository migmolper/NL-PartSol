function post(velo,tp,ste_p,Sup_energy,K_hammer,K_beam,W_hammer,W_beam,D_energy,React,Impact,React2,Impact2)

if velo==2640
    load EXP2640;
elseif velo==1760
    load EXP1760;
elseif velo==881
    load EXP881;
end
del=0.1/velo;


figure;
plot(tp(1:ste_p),(K_hammer(1)-K_hammer(1:ste_p)-W_hammer(1:ste_p))/1000,...
    tp(1:ste_p),(K_beam(1:ste_p)+D_energy(1:ste_p)-Sup_energy(1:ste_p)+W_beam(1:ste_p))/1000);
xlabel('Time (s)');
ylabel('Energy [N·m]');
legend('Input energy','Total energy');

figure;
plot(tp(1:ste_p),K_beam(1:ste_p)/1000,tp(1:ste_p),...
    W_beam(1:ste_p)/1000,tp(1:ste_p),D_energy(1:ste_p)/1000,...
    tp(1:ste_p),-Sup_energy(1:ste_p)/1000);
xlabel('Time (s)');
ylabel('Beam Energy [N·m]');
legend('Kinetic en.','Strain en.','Dissipated en.','Support en.');


figure;
plot(tp(1:ste_p),(K_hammer(1)-K_hammer(1:ste_p))/1000,tp(1:ste_p),...
    W_hammer(1:ste_p)/1000);
xlabel('Time (s)');
ylabel('Hammer Energy [N·m]');
legend('Kinetic en.','Strain en.');

figure;
plot(tp(1:ste_p),React(1:ste_p)/1000,tp(1:ste_p),-Impact(1:ste_p)/1000);
xlabel('Time (s)');
ylabel('Forces [kN]');
legend('Reaction','Impact');

figure;
plot(tp(1:ste_p),React(1:ste_p)/1000,tp(1:ste_p),-Impact(1:ste_p)/1000,...
    tp(1:ste_p),React2(1:ste_p)/1000,tp(1:ste_p),Impact2(1:ste_p)/1000);
xlabel('Time (s)');
ylabel('Forces [kN]');
legend('Reaction','Impact','Reaction2','Impact2');

if velo==0
    figure1 = figure;
    axes1 = axes('Parent',figure1);
    plot(tp(1:ste_p),React(1:ste_p)/1000,tp(1:ste_p),-Impact(1:ste_p)/1000,...
        EXP(:,1)+del,EXP(:,3),EXP(:,1)+del,EXP(:,2));
    xlabel('Time (s)');
    ylabel('Forces [kN]');
    legend('Reaction','Impact','EXP Reaction','EXP Impact');
    xlim(axes1,[0 0.0008]);
end


