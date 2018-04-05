%% deploy nodes
close all;
clc;
clear;

% lambda is the intensity 
lambda=15;

nmb=poissrnd(lambda);

x = rand(1,nmb);
y = rand(1,nmb);
x = 4*x + 1;
y = 4*y + 1;

%coordinates of fence nodes
%fence_x1 = [0,0,0,0,0, 0, 0,2, 2,4, 4,6, 6,8, 8,10,10,12,12,12,12,12,12,12];
%fence_y1 = [0,2,4,6,8,10,12,0,12,0,12,0,12,0,12, 0,12, 0, 2, 4, 6, 8,10,12];
fence_x1=[0,0,0,0,2,2,4,4,6,6,6,6];
fence_y1=[0,2,4,6,0,6,0,6,0,2,4,6];

%inner_x1 = [1,1,1,1,1, 1,3, 3,5, 5,7, 7,9, 9,11,11,11,11,11,11];
%inner_y1 = [1,3,5,7,9,11,1,11,1,11,1,11,1,11, 1, 3, 5, 7, 9,11];
inner_x1=[1,1,1,3,3,5,5,5];
inner_y1=[1,3,5,1,5,1,3,5];

%coordinates of all the nodes
node_x = [x, inner_x1, fence_x1];
node_y = [y, inner_y1, fence_y1]; 

node_coor = [node_x; node_y];
sleep_node=[];

subplot(2,2,1);
plot(node_x, node_y, '.'); hold on;

subplot(2,2,2);
plot(node_x, node_y, '.'); hold on;

% figure
Theta=[0:0.005:1]*2*pi;
Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
Ycircle=1.0*sin(Theta);
for i=1:length(node_coor)-12        % 画出fence node之外所有节点的Rs  /24
    Xc=Xcircle+node_x(i);
    Yc=Ycircle+node_y(i);
    subplot(2,2,1);
    plot(Xc,Yc,'k');
    fill(Xc,Yc,'g','facealpha',0.5);  %半透明显示
    
    axis square; %产生正方形坐标系
    xlim([0 6]); %x,y轴上下限设置
    ylim([0 6]);
    hold on;
end

hold on;

%打印节点序号
for i = 1: length(node_coor)-12 
    text(node_x(i)+0.1, node_y(i), num2str(i));
end

for i=1: length(node_coor) 
    node(i).neighbor = [];
    node(i).weight=0; 
    node(i).simp=[];
    node(i).status=1;
    
    if i > length(node_coor) - 12
        node(i).fence_flag = 1;
    else
        node(i).fence_flag = 0;    
    end
    
    if i>length(node_coor)-20
        node(i).c_flag=0;
    else
        node(i).c_flag=1;
    end
    
    for j=1: length(node_coor)
        if (j==i) continue;
        else
            distance = dist2(node_coor(:,i),node_coor(:,j));
            if distance <= 2 
                node(i).neighbor = [node(i).neighbor, j];%找到 node(i)的所有邻居节点  即构建出1-simplex
            end
        end
    end
end

%%
% construct simplices
for i = 1 : length(node_coor)
%     construct 1-simplex (edge)
    no_neighb = length(node(i).neighbor);
    if no_neighb > 0
        for j = 1 : no_neighb
            new_node = node(i).neighbor(j);
            vert_set_temp = sort([i, new_node]);
            neighb_set_temp = intersect(node(i).neighbor, node(new_node).neighbor);
            simp_temp = struct('vert', vert_set_temp, 'neighb', neighb_set_temp);

            node(i).simp{1}(j) = simp_temp;
        end
    else
        node(i).simp{1} = [];
    end
    
%     construct 2-simplex (triangle)
    if ~isempty(node(i).simp{1})
        no_edge = size(node(i).simp{1}, 2);
        for j = 1 : no_edge
            vert_set = node(i).simp{1}(j).vert;
            neighb_set = node(i).simp{1}(j).neighb;
            if ~isempty(neighb_set)
                no_neighb = length(neighb_set);
                for k = 1 : no_neighb
                    new_node = neighb_set(k);
                    if new_node < max(setdiff(vert_set, i))
                        continue;
                    else
                        vert_set_temp = union(vert_set, new_node);
                        neighb_set_temp = intersect(neighb_set, node(new_node).neighbor);
                        simp_temp = struct('vert', vert_set_temp, 'neighb', neighb_set_temp);

                        if size(node(i).simp, 2) == 1 %this is the first new higher simplex
                            node(i).simp{2}(1) = simp_temp;
                        else
                            no_new_simp = size(node(i).simp{2}, 2);
                            node(i).simp{2}(no_new_simp+1) = simp_temp;
                        end
                    end
                end
            end
        end
    end
end

for i = 1 : length(node_coor)
    dim_max = size(node(i).simp, 2);
    if dim_max == 2
        no_tri = size(node(i).simp{2}, 2);
        for j = 1 : no_tri
            vert_set = node(i).simp{2}(j).vert;
            if i == min(vert_set)
                node_ID1 = vert_set(1);
                node_ID2 = vert_set(2);
                node_ID3 = vert_set(3);
                X = [node_x(node_ID1), node_x(node_ID2), node_x(node_ID3)];
                Y = [node_y(node_ID1), node_y(node_ID2), node_y(node_ID3)];
                axis square; %产生正方形坐标系
                xlim([0 6]); %x,y轴上下限设置
                ylim([0 6]);
                subplot(2,2,2);
                fill(X, Y, [204/255 204/255 204/255]); hold on;
            else 
                continue;
            end
        end
    end
end

for i=1: length(node_coor)
    neigh_set = node(i).neighbor;
    for j = 1 : length(node(i).neighbor)
        neigh_ID = node(i).neighbor(j);
        if neigh_ID > i
            line([node_x(i),node_x(neigh_ID)],[node_y(i),node_y(neigh_ID)], 'Color', 'k');
        end
    end
end
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);


%%

for i=1:length(node_coor)   % 节点i的邻居w是否能够被collapse
    if node(i).status==1
       sleep_temp=[];
       for j=1:length(node(i).neighbor)
          w=node(i).neighbor(j);
          if node(w).status==1 && node(w).c_flag==1
              co_neighbor=intersect(node(i).neighbor,node(w).neighbor);
              diff_node=setdiff(node(w).neighbor,co_neighbor);
              if diff_node==i
                  node(w).status=0;
                  node(i).c_flag=0;
                  disp(['node ' num2str(w) ' is dominated by node ' num2str(i)]);
                  sleep_node=[sleep_node,w];
                  sleep_temp=[sleep_temp,w];
              end               
          end
       end
       for m=1:length(sleep_temp)
           for n=1:length(node(sleep_temp(m)).neighbor)
               node_temp=node(sleep_temp(m)).neighbor(n);
               node(node_temp).neighbor=setdiff(node(node_temp).neighbor,sleep_temp(m));
           end
       end
    end
end


subplot(2,2,3);
Theta=[0:0.005:1]*2*pi;
Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
Ycircle=1.0*sin(Theta);
for i=1:length(node_coor)-12        % 画出fence node之外所有节点的Rs  /24
    if node(i).status==1
        plot(node_x(i), node_y(i), '.'); hold on;
        
        Xc=Xcircle+node_x(i);
        Yc=Ycircle+node_y(i);
        plot(Xc,Yc,'k');
        fill(Xc,Yc,'g','facealpha',0.5);  %半透明显示
        text(node_x(i)+0.1, node_y(i), num2str(i));
    
        axis square; %产生正方形坐标系
        xlim([0 6]); %x,y轴上下限设置
        ylim([0 6]);
        hold on;
    end
end

%%
% construct new simplices
for i = 1 : length(node_coor)
%     construct 1-simplex (edge)
  node(i).simp=[];  
  if node(i).status==1
    no_neighb = length(node(i).neighbor);
    if no_neighb > 0
        for j = 1 : no_neighb
            new_node = node(i).neighbor(j);
            vert_set_temp = sort([i, new_node]);
            neighb_set_temp = intersect(node(i).neighbor, node(new_node).neighbor);
            simp_temp = struct('vert', vert_set_temp, 'neighb', neighb_set_temp);

            node(i).simp{1}(j) = simp_temp;
        end
    else
        node(i).simp{1} = [];
    end
    
%     construct 2-simplex (triangle)
    if ~isempty(node(i).simp{1})
        no_edge = size(node(i).simp{1}, 2);
        for j = 1 : no_edge
            vert_set = node(i).simp{1}(j).vert;
            neighb_set = node(i).simp{1}(j).neighb;
            if ~isempty(neighb_set)
                no_neighb = length(neighb_set);
                for k = 1 : no_neighb
                    new_node = neighb_set(k);
                    if new_node < max(setdiff(vert_set, i))
                        continue;
                    else
                        vert_set_temp = union(vert_set, new_node);
                        neighb_set_temp = intersect(neighb_set, node(new_node).neighbor);
                        simp_temp = struct('vert', vert_set_temp, 'neighb', neighb_set_temp);

                        if size(node(i).simp, 2) == 1 %this is the first new higher simplex
                            node(i).simp{2}(1) = simp_temp;
                        else
                            no_new_simp = size(node(i).simp{2}, 2);
                            node(i).simp{2}(no_new_simp+1) = simp_temp;
                        end
                    end
                end
            end
        end
    end
  end
end

for i = 1 : length(node_coor)
    dim_max = size(node(i).simp, 2);
    if dim_max == 2
        no_tri = size(node(i).simp{2}, 2);
        for j = 1 : no_tri
            vert_set = node(i).simp{2}(j).vert;
            if i == min(vert_set)
                node_ID1 = vert_set(1);
                node_ID2 = vert_set(2);
                node_ID3 = vert_set(3);
                X = [node_x(node_ID1), node_x(node_ID2), node_x(node_ID3)];
                Y = [node_y(node_ID1), node_y(node_ID2), node_y(node_ID3)];
                subplot(2,2,4);
                axis square; %产生正方形坐标系
                xlim([0 6]); %x,y轴上下限设置
                ylim([0 6]);
                fill(X, Y, [204/255 204/255 204/255]); hold on; %204/255 is color grey
            else 
                continue;
            end
        end
    end
end

for i=1: length(node_coor)
  if node(i).status==1
    neigh_set = node(i).neighbor;
    for j = 1 : length(node(i).neighbor)
        neigh_ID = node(i).neighbor(j);
        if neigh_ID > i
            line([node_x(i),node_x(neigh_ID)],[node_y(i),node_y(neigh_ID)], 'Color', 'k');
        end
    end
  end
end
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);