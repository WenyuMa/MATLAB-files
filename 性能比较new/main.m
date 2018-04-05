function[error,hole_num]=main()
%clc;
clear;

% lambda is the intensity 
lambda=130;
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
for i=1:length(node_coor)-24         % 画出fence node之外所有节点的Rs
    Xc=Xcircle+node_x(i);
    Yc=Ycircle+node_y(i);
    plot(Xc,Yc,'k');
    fill(Xc,Yc,'g');
    
    axis square; %产生正方形坐标系
    xlim([0 12]); %x,y轴上下限设置
    ylim([0 12]);
    hold on;
end

hold on;

%打印节点序号
for i = 1: length(node_coor)   
    text(node_x(i)+0.1, node_y(i), num2str(i));
end


result1=boundary();
result2=boundary2();  %参考算法结果

cycle_flag1=zeros(1,length(result1));
cycle_flag2=zeros(1,length(result2));
for i=1:length(result2)
    for j=1:length(result1)  
      if cycle_flag1(j)==0 && cycle_flag2(i)==0
        hole1=unique(result1{j});
        hole2=unique(result2{i});
        con_node=intersect(hole1,hole2);
        s1=isequal(sort(hole1),con_node);
        s2=isequal(sort(hole2),con_node);
        if s1==1||s2==1
            cycle_flag1(j)=1;    %匹配成功，置1
            cycle_flag2(i)=1;
            break;
        end
      end 
    end
end

error=sum(cycle_flag2(:)==0);
hole_num=length(result2);
% if error
%    disp('error');     %breakpoint
% end


