%collapse nodes and edges based on neighbor set, neighbor set include node itself%

%% This is the main function
% 1. set nodes, bulid G=(V,E)
% 2. calculate 1 and 2-simplex to bulid the homology figure
% 3. strong collapse nodes, recompute the simplex 
% 4. strong collapse edges, recompute the simplex
% 5. repeat step 3 and 4 until there is no node and edge can be delete.

%% deploy nodes
close all; 
clc;
clear;

global node node_coor node_x node_y edge

% lambda is the intensity 
lambda=35;

nmb=poissrnd(lambda);

x = rand(1,nmb);
y = rand(1,nmb);
x = 6*x + 1;
y = 6*y + 1;

fence_x1 = [0,0,0,0,0,2,2,4,4,6,6,8,8,8,8,8];  %16
fence_y1 = [0,2,4,6,8,0,8,0,8,0,8,0,2,4,6,8];

inner_x1 = [1,1,1,1,3,3,5,5,7,7,7,7];  %12
inner_y1 = [1,3,5,7,1,7,1,7,1,3,5,7];

%coordinates of all the nodes
node_x = [x, inner_x1, fence_x1];
node_y = [y, inner_y1, fence_y1]; 

node_coor = [node_x; node_y];

subplot(2,2,1);
plot(node_x, node_y, '.'); 
title('The original network');
hold on;

subplot(2,2,2);
plot(node_x, node_y, '.'); 
title('Homology of the original network');
hold on;


% figure
Theta=[0:0.005:1]*2*pi;
Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
Ycircle=1.0*sin(Theta);
for i=1:length(node_coor)-16        % 画出fence node之外所有节点的Rs  /24
    Xc=Xcircle+node_x(i);
    Yc=Ycircle+node_y(i);
    subplot(2,2,1);
    plot(Xc,Yc,'k');
    fill(Xc,Yc,'g','facealpha',0.5);  %半透明显示
    text(node_x(i)+0.1, node_y(i), num2str(i));
    hold on;
end

axis square; %产生正方形坐标系
xlim([0 8]); %x,y轴上下限设置
ylim([0 8]);
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

edge.u=[];
edge.v=[];
edge.c_flag=[];

%% initial
for i=1: length(node_coor) 
    node(i).neighbor = [];
    node(i).weight=0; 
    node(i).simp=[];
    node(i).status=1;
    
    if i > length(node_coor)-16
        node(i).fence_flag = 1;
    else
        node(i).fence_flag = 0;    
    end
    
    if i>length(node_coor)-28
        node(i).c_flag=0;
    else
        node(i).c_flag=1;
    end
    
    for j=1: length(node_coor)
        %if (j==i) continue;
        %else
            distance = dist2(node_coor(:,i),node_coor(:,j));
            if distance <= 2 
                node(i).neighbor = [node(i).neighbor, j];%找到 node(i)的所有邻居节点  即构建出1-simplex
            end
        %end
    end
end

for i = 1: length(node_coor) - 16 
    if node(i).status==1
        no_neighbor = length(node(i).neighbor);
        if no_neighbor 
            for j = 1: no_neighbor
                if node(i).neighbor(j) <= i   % u>=v,continue || node(node(i).neighbor(j)).fence_flag == 1
                    continue;
                elseif node(node(i).neighbor(j)).status==1
                    edge.u = [edge.u,i];
                    edge.v=[edge.v,node(i).neighbor(j)];
                    edge.c_flag=[edge.c_flag,1];
                end 
            end 
        end 
    end
end

%% construct simplices
for i = 1 : length(node_coor)
   % construct 1-simplex (edge)
    neighbor_temp1=setdiff(node(i).neighbor,i);
    new_simp = struct('vert', [i], 'neighb', neighbor_temp1);
    node(i).simp{1} = new_simp;
    no_neighb = length(neighbor_temp1);
    
    if no_neighb > 0
        for j = 1 : no_neighb
            new_node = neighbor_temp1(j);
            vert_set_temp = sort([i, new_node]);
            neighbor_temp2=setdiff( node(new_node).neighbor,new_node);
            neighb_set_temp = intersect(neighbor_temp1, neighbor_temp2);
            simp_temp = struct('vert', vert_set_temp, 'neighb', neighb_set_temp);

            node(i).simp{2}(j) = simp_temp;
        end
    else
        node(i).simp{2} = [];
    end
    
   % construct 2-simplex (triangle)
    if ~isempty(node(i).simp{2})
        no_edge = size(node(i).simp{2}, 2);
        for j = 1 : no_edge
            vert_set = node(i).simp{2}(j).vert;
            neighb_set = node(i).simp{2}(j).neighb;
            if ~isempty(neighb_set)
                no_neighb = length(neighb_set);
                for k = 1 : no_neighb
                    new_node = neighb_set(k);
                    if new_node < max(setdiff(vert_set, i))
                        continue;
                    else
                        vert_set_temp = union(vert_set, new_node);
                        neighb_set_temp = intersect(neighb_set, setdiff(node(new_node).neighbor,new_node));
                        simp_temp = struct('vert', vert_set_temp, 'neighb', neighb_set_temp);

                        if size(node(i).simp, 2) == 2 %this is the first new higher simplex
                            node(i).simp{3}(1) = simp_temp;
                        else
                            no_new_simp = size(node(i).simp{3}, 2);
                            node(i).simp{3}(no_new_simp+1) = simp_temp;
                        end
                    end
                end
            end
        end
    end
end

%% figure
subplot(2,2,2);
for i = 1 : length(node_coor)
    dim_max = size(node(i).simp, 2);
    if dim_max == 3
        no_tri = size(node(i).simp{3}, 2);
        for j = 1 : no_tri
            vert_set = node(i).simp{3}(j).vert;
            if i == min(vert_set)
                node_ID1 = vert_set(1);
                node_ID2 = vert_set(2);
                node_ID3 = vert_set(3);
                X = [node_x(node_ID1), node_x(node_ID2), node_x(node_ID3)];
                Y = [node_y(node_ID1), node_y(node_ID2), node_y(node_ID3)];
                fill(X, Y, [204/255 204/255 204/255]); hold on;
            else 
                continue;
            end
        end
    end
end

for i=1: length(node_coor)
    text(node_x(i)+0.1, node_y(i), num2str(i));
    neigh_set = node(i).neighbor;
    for j = 1 : length(node(i).neighbor)
        neigh_ID = node(i).neighbor(j);
        if neigh_ID > i
            line([node_x(i),node_x(neigh_ID)],[node_y(i),node_y(neigh_ID)], 'Color', 'k');
        end
    end
end

axis square; %产生正方形坐标系
xlim([0 8]); %x,y轴上下限设置
ylim([0 8]);
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

flag1=NodeCollapse();
flag2=EdgeCollapse();
ConstructSimp();

%% figure

subplot(2,2,3);

for i = 1 : length(node_coor)
    dim_max = size(node(i).simp, 2);
    if dim_max > 2
        no_tri = size(node(i).simp{3}, 2);
        for j = 1 : no_tri
            vert_set = node(i).simp{3}(j).vert;
            if i == min(vert_set)
                node_ID1 = vert_set(1);
                node_ID2 = vert_set(2);
                node_ID3 = vert_set(3);
                X = [node_x(node_ID1), node_x(node_ID2), node_x(node_ID3)];
                Y = [node_y(node_ID1), node_y(node_ID2), node_y(node_ID3)];
                fill(X, Y, [204/255 204/255 204/255]); hold on; %204/255 is color grey
            else 
                continue;
            end
        end
    end
end

for i=1: length(node_coor)
  if node(i).status==1
      text(node_x(i)+0.1, node_y(i), num2str(i));
      neigh_set = node(i).neighbor;
      for j = 1 : length(node(i).neighbor)
          neigh_ID = node(i).neighbor(j);
          if neigh_ID > i
              line([node_x(i),node_x(neigh_ID)],[node_y(i),node_y(neigh_ID)], 'Color', 'k');
          end
      end 
  end 
end

axis square; %产生正方形坐标系
xlim([0 8]); %x,y轴上下限设置
ylim([0 8]);
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);
title('Homology of the network after first strong collapse');

while(flag1 || flag2)
    flag1=NodeCollapse();
    flag2=EdgeCollapse();    
end

%flag2=EdgeCollapse();
ConstructSimp();
subplot(2,2,4);  %%%%%

for i = 1 : length(node_coor)
    dim_max = size(node(i).simp, 2);
    if dim_max > 2
        no_tri = size(node(i).simp{3}, 2);
        for j = 1 : no_tri
            vert_set = node(i).simp{3}(j).vert;
            if i == min(vert_set)
                node_ID1 = vert_set(1);
                node_ID2 = vert_set(2);
                node_ID3 = vert_set(3);
                X = [node_x(node_ID1), node_x(node_ID2), node_x(node_ID3)];
                Y = [node_y(node_ID1), node_y(node_ID2), node_y(node_ID3)];
                fill(X, Y, [204/255 204/255 204/255]); hold on; %204/255 is color grey
            else 
                continue;
            end
        end
    end
end

for i=1: length(node_coor)
  if node(i).status==1
      text(node_x(i)+0.1, node_y(i), num2str(i));
      neigh_set = node(i).neighbor;
      for j = 1 : length(node(i).neighbor)
          neigh_ID = node(i).neighbor(j);
          if neigh_ID > i
              line([node_x(i),node_x(neigh_ID)],[node_y(i),node_y(neigh_ID)], 'Color', 'k');
          end
      end 
  end 
end

axis square; %产生正方形坐标系
xlim([0 8]); %x,y轴上下限设置
ylim([0 8]);
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);
title('Homology of the network after several strong collapse');

