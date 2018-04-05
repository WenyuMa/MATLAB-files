%% deploy nodes
clc;
clear;

% lambda is the intensity 
lambda=110;

nmb=poissrnd(lambda);

x = rand(1,nmb);
y = rand(1,nmb);
x = 10*x + 1;
y = 10*y + 1;

%coordinates of fence nodes
fence_x1 = [0,0,0,0,0, 0, 0,2, 2,4, 4,6, 6,8, 8,10,10,12,12,12,12,12,12,12];
fence_y1 = [0,2,4,6,8,10,12,0,12,0,12,0,12,0,12, 0,12, 0, 2, 4, 6, 8,10,12];

inner_x1 = [1,1,1,1,1, 1,3, 3,5, 5,7, 7,9, 9,11,11,11,11,11,11];
inner_y1 = [1,3,5,7,9,11,1,11,1,11,1,11,1,11, 1, 3, 5, 7, 9,11];

%coordinates of all the nodes
node_x = [x, inner_x1, fence_x1];
node_y = [y, inner_y1, fence_y1]; 

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
    fill(Xc,Yc,'g','facealpha',0.5);  %半透明显示
    
    axis square; %产生正方形坐标系
    xlim([0 12]); %x,y轴上下限设置
    ylim([0 12]);
    hold on;
end

hold on;

%打印节点序号
for i = 1: length(node_coor)-24  
    text(node_x(i)+0.1, node_y(i), num2str(i));
end

% set the fence_flag for all nodes to indicate whether they are fence node
% or not

for i=1: length(node_coor) 
    node(i).neighbor = [];
    node(i).weight=0; 
    node(i).simp1=[];
    node(i).status=1;
    
    if i > length(node_coor) - 24
        node(i).fence_flag = 1;
    else
        node(i).fence_flag = 0;    
    end
    
    for j=1: length(node_coor)
        if (j==i) continue;
        else
            distance = dist2(node_coor(:,i),node_coor(:,j));
            if distance <= 2 
                node(i).neighbor = [node(i).neighbor, j];%找到 node（i)的所有邻居节点  即构建出1-simplex
                 %line([node_x(i),node_x(j)],[node_y(i),node_y(j)]);
            end
        end
    end
end

%% find simplex
% find 1-simplex and 2-simplex
for i = 1 : length(node_coor)
  if node(i).status==1
      
      node(i).simp1 = node(i).neighbor;
      node(i).simp2 = cell(length(node(i).neighbor),1);    
    
      for j=1:length(node(i).neighbor)       
          u=node(i).neighbor(j);
          index1=ismember(node(i).neighbor,node(u).neighbor);
          co_neighbor=node(i).neighbor(index1);
          node(i).simp2{j}=co_neighbor;
      end
  end
end

%% calculate node weight
for i = 1 : length(node_coor)
    if node(i).fence_flag==1 
        node(i).weight=0;
    else
        flag1=0;  % whether exist one edge in E(v) that does not belong to any triangle
        for j=1:length(node(i).simp2)
            flag1=isempty(node(i).simp2{j});
            if flag1==1
                node(i).weight=0;
                break;
            end
        end
        if (flag1==0)
            flag2=0;  % There exist a triangle that does not have any neighbor;
            for m=1:length(node(i).simp2)
                for n=1:length(node(i).simp2{m})
                    v=node(i).simp2{m}(n);
                    t_neighbor=intersect(node(v).neighbor,node(i).simp2{m});
                    if isempty(t_neighbor)
                        flag2=1;
                        break;
                    end
                end
                if flag2==1
                    break;
                else
                    continue;
                end
            end
            if flag2==0
                node(i).weight=3;
            else
                node(i).weight=2;
            end
        end
    end
end


 

