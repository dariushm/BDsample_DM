matX0 = dlmread('rate_at_x0_1_dt1.68885e-08.dat');
matXdx = dlmread('rate_at_x0_1.02513_dt1.68885e-08.dat');

matBDsim = (matX0+matXdx)/2;

ktonlyDtnoDr = dlmread('rate_for_Dtrans_only_dt1.68885e-08.dat');
ktreal = dlmread('rate_vs_time_calc_dt1.68885e-08.dat');

semilogx(matBDsim(:,1),matBDsim(:,2),'r','LineWidth',2)
hold on
semilogx(ktreal(:,1),ktreal(:,2),'b','LineWidth',2)
semilogx(ktonlyDtnoDr(:,1),ktonlyDtnoDr(:,2),'b--','LineWidth',2)


xlabel('Time (s.)','FontSize',16);
ylabel('k(t) (nm^3/{\mu}s)','FontSize',16);
set(gca,'FontSize',16);
h=legend('Simulation','Theory','Theory with D_{trans} only','Location','NorthEast');
set(h,'FontSize',12);
set(gcf,'Color','w')
saveas(gcf, 'kttest.fig', 'fig')
