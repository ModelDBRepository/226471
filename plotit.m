load('time_tc');
plot(time_tc(:,1),time_tc(:,2))
xlim([2000 4000])
ylim([-58 -50])
title('PN LFP')
xlabel('Time [ms]')
ylabel('Voltage [mV]')
saveas(gca,'LFP.png','png')
