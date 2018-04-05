%y=[11777 11748 13469 13359; 10880 10864  13473 13405 ;  9755 9753  12999 12913;  8732 8724  12310 12260; 7622 7618  11647 11610];
y=[12427 12357; 11777 11748; 10880 10864; 9755 9753; 8732 8724;];
b=bar(y);
grid on;
ch = get(b,'children');
set(gca,'XTickLabel',{'0.008','0.009','0.010','0.011','0.012'})

xlabel('Node density \lambda ');
ylabel('Number of coverage holes');

legend('# of coverage holes/LBA','# of coverage holes/RBA');

% 
% data = [11748 11777 13359 13469; 10864 10880 13405 13473 ; 9753 9755 12913 12999; 8724 8732 12260 12310;7618 7622 11610 11647];;
%  bar(data,1)
%  %axis([0 6 0.0 100])
%  legend('# of coverage holes/RBA(k=1)','# of coverage holes/LBA(k=1)','# of coverage holes/RBA(k=2)','# of coverage holes/LBA(k=2)')
%  set(gca,'XTickLabel',{'0.009','0.010','0.011','0.012','0.013'})
%  applyhatch(gcf,'\+x.')

% x=[0.008 0.009 0.010 0.011 0.012 0.013];
% y1=[2982.18 3489.94 4223.99 5048.56 5688.13];
% y2=[1045.394 1395.022 1838.512 2327.598 2737.584];
% 
% y_all=[y1;y2]';
% bar(x,y_all)
% title('1995—2009产业总产值')
% xlabel('年份')
% ylabel('亿元')
% legend('y1','y2',2)