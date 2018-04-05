y=[12357 12427; 11748 11777; 10864 10880; 9753 9755; 8724 8732;];
b=bar(y);
grid on;
ch = get(b,'children');
set(gca,'XTickLabel',{'0.8','0.9','1','1.1','1.2'})
set(ch{1},'FaceColor',[0 0 0])
set(ch{2},'FaceColor',[0.5 0.5 0.5])
legend([ch{1},ch{2}],'# of coverage holes found by RBA','# of coverage holes found by LBA');
xlabel('intensity \lambda ');
ylabel('number of coverage holes');