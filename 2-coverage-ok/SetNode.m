% This function is used for setting a WSN and calculate nodes' neighbors

function []=SetNode()

global k  node node_coor node_x node_y

%设置Rs=1，Rc=2*Rs，节点在12*12的正方形区域内分布

node_coor = [node_x; node_y];

% plot(node_x, node_y,'.','markersize',10,'color',[0 0 0]); hold on;


% figure
Theta=[0:0.005:1]*2*pi;
Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
Ycircle=1.0*sin(Theta);
for i=1:length(node_coor)-24*k         % 画出fence node之外所有节点的Rs
    Xc=Xcircle+node_x(i);
    Yc=Ycircle+node_y(i);
    plot(Xc,Yc,'k');
    fill(Xc,Yc,[0.5,0.5,0.5],'facealpha',0.8);  %半透明显示
    plot(node_x(i), node_y(i),'.','markersize',10,'color',[0 0 0]);
    
    hold on;
end

axis square; %产生正方形坐标系
xlim([0 12]); %x,y轴上下限设置
ylim([0 12]);
hold on;

% plot(node_x, node_y,'.','markersize',10,'color',[0 0 0]); hold on;

%打印节点序号
% for i = 1: length(node_coor)-24*k  
%     text(node_x(i)+0.1, node_y(i), num2str(i));
% end

% set the fence_flag for all nodes to indicate whether they are fence node
% or not
for i = 1: length(node_coor) 
    node(i).status = 1;   % sleep or active
    node(i).bn=0;         % a boundary node or not
    if i > length(node_coor) - 24*k 
        node(i).fence_flag = 1;
    else
        node(i).fence_flag = 0;    
    end
end


%相邻节点间的信息计算
for i=1: length(node_coor)         
    node(i).neighbour = [];
    node(i).dist = [];
    node(i).angle = [];
    node(i).msg=[];

    for j=1: length(node_coor)
        if (j==i) 
            continue;
        else
            distance = dist2(node_coor(:,i),node_coor(:,j));
            if distance <= 2 
                node(i).neighbour = [node(i).neighbour, j];
                node(i).dist = [node(i).dist, distance];
                node(i).angle = [node(i).angle, 2*acos(distance/2)];
            end
        end
    end
end
