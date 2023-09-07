close all
ftick=[0.25 0.5 1 2 5 10 20 40 65];
for j=1:3
    figure()
    plot((fc),mpmat(:,j,1),'-bs','LineWidth',1,'MarkerEdgeColor','b', 'MarkerSize',4);
    hold on
    plot((fc),mpmat(:,j,2),'-rd','LineWidth',1,'MarkerEdgeColor','r', 'MarkerSize',4);
    set(gcf,'position',[400,400,400,300])
    xlim([0.25 70]);
    set(gca, 'XScale', 'log')
    set(gca,'XTick',ftick)
    legend('VKOT','SVKOT');
    ylabel(append('MPAE',num2str(j)));xlabel('cutoff freq.')
end
figure()
plot((fc),sum(mpmat(:,:,1),2)/3,'-bs','LineWidth',1,'MarkerEdgeColor','b', 'MarkerSize',4);
hold on
plot((fc),sum(mpmat(:,:,2),2)/3,'-rd','LineWidth',1,'MarkerEdgeColor','r', 'MarkerSize',4);
set(gcf,'position',[400,400,400,300])
xlim([0.25 70]);
set(gca, 'XScale', 'log')
set(gca,'XTick',ftick)
legend('VKOT','SVKOT');
ylabel('Avg. MPAE');xlabel('cutoff freq.')