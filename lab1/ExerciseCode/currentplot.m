load ECM1
deffilem5
%subplot(211),
plot(Time/60,Current)

axis([0 40 -300 200])
xlabel('Time (min)')
ylabel('Current (A)')
text(25,-200,'Discharge')
text(25,130,'Charge')
%print -depsc currentplot
