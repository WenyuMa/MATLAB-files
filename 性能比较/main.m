function[error]=main()
%clc;
clear;

% lambda is the intensity 
lambda=80;
nmb=poissrnd(lambda);

global node_coor node_x node_y

x = rand(1,nmb);
y = rand(1,nmb);
x = 10*x + 1;
y = 10*y + 1;

%coordinates of fence nodes
fence_x = [0,0,0,0,0, 0, 0,2, 2,4, 4,6, 6,8, 8,10,10,12,12,12,12,12,12,12];
fence_y = [0,2,4,6,8,10,12,0,12,0,12,0,12,0,12, 0,12, 0, 2, 4, 6, 8,10,12];

inner_x = [1,1,1,1,1, 1,3, 3,5, 5,7, 7,9, 9,11,11,11,11,11,11];
inner_y = [1,3,5,7,9,11,1,11,1,11,1,11,1,11, 1, 3, 5, 7, 9,11];

%coordinates of all the nodes
node_x = [x, inner_x, fence_x];
node_y = [y, inner_y, fence_y];


node_coor = [node_x; node_y];

plot(node_x, node_y, '.'); hold on;


% figure
Theta=[0:0.005:1]*2*pi;
Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
Ycircle=1.0*sin(Theta);
for i=1:length(node_coor)-24         % ����fence node֮�����нڵ��Rs
    Xc=Xcircle+node_x(i);
    Yc=Ycircle+node_y(i);
    plot(Xc,Yc,'k');
    fill(Xc,Yc,'g');
    
    axis square; %��������������ϵ
    xlim([0 12]); %x,y������������
    ylim([0 12]);
    hold on;
end

hold on;

%��ӡ�ڵ����
for i = 1: length(node_coor)   
    text(node_x(i)+0.1, node_y(i), num2str(i));
end


result1=boundary();
result2=boundary2();  %�ο��㷨���
error=0;
s=0;
if length(result1)==length(result2)   %�ն������Ƿ����
   for i=1:length(result2)
       for j=1:length(result1)          
           if length(result1{j})==length(result2{i})+1 
               hole1=unique(result1{j});
               hole2=unique(result2{i});
               s=isequal(sort(hole1),sort(hole2));
               if s==1
                  break;
               end 
           end
       end
       if s==0
           error=1;   %% breakpoint1
           break;
       end
   end
else
    error=1;    %% breakpoint2
end
