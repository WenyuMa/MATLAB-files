%% deploy nodes
close all;
clc;
clear;

% lambda is the intensity 
lambda=10;

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

plot(node_x, node_y, '.'); hold on;

% figure
Theta=[0:0.005:1]*2*pi;
Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
Ycircle=1.0*sin(Theta);
for i=1:length(node_coor)        % 画出fence node之外所有节点的Rs  /24
    Xc=Xcircle+node_x(i);
    Yc=Ycircle+node_y(i);
    plot(Xc,Yc,'k');
    fill(Xc,Yc,'g','facealpha',0.5);  %半透明显示
    
    axis square; %产生正方形坐标系
    xlim([0 6]); %x,y轴上下限设置
    ylim([0 6]);
    hold on;
end

hold on;

%打印节点序号
for i = 1: length(node_coor) 
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
    
    for j=1: length(node_coor)
        if (j==i) continue;
        else
            distance = dist2(node_coor(:,i),node_coor(:,j));
            if distance <= 2 
                node(i).neighbor = [node(i).neighbor, j];%找到 node（i)的所有邻居节点  即构建出1-simplex
            end
        end
    end
end

% ？？？？？？？？？
for i = 1 : length(node_coor)
    new_simp = struct('vert', [i], 'neighb', node(i).neighbor);
    node(i).simp{1} = new_simp;
end
%然后依次递归 找出每个顶点的k-simplex
for i = 1 : length(node_coor)
    while(1)
        index_max = size(node(i).simp, 2);%最大的单形,返回矩阵的列数

        no_max_simp = size(node(i).simp{index_max}, 2); %????
        %最大单形的最大数目
%考察每个顶点的最大维单形有no_max_simp个   对于每个最大维单形
% 对于
        for j = 1 : no_max_simp
            vert_set = node(i).simp{index_max}(j).vert;%第j个最大维单形的顶点的集合 计算j-1单形的顶点集合    (从0-simplex开始计算）
            neighb_set = node(i).simp{index_max}(j).neighb;%第j个最大维单形的邻居集合

            if ~isempty(neighb_set)
                
                no_neighb = length(neighb_set);
                for k = 1 : no_neighb     %考察第j个最大维单形 邻居集合的每个顶点
                    new_node = neighb_set(k);%对于构成该第j个单形的邻居集合的第k个顶点
                    if new_node < max(setdiff(vert_set, i))  %setdiff(vert_set, i)返回vert_set中存在而i中不存在的元素
                        continue;
                    else
                       
                        vert_set_temp = union(vert_set, new_node);
                        neighb_set_temp = intersect(neighb_set, node(new_node).neighbor);%当前i顶点的邻居集合  和  当前i顶点 第k个邻居顶点的邻居集合  的交集   
                        simp_temp = struct('vert', vert_set_temp, 'neighb', neighb_set_temp);%由第i个顶点和第k个顶点构成  以及它们共同的邻居顶点构成的higher simplex

                        if index_max == size(node(i).simp, 2) %this is the first index_max+1 new higher simplex 代表比index_max更高维的复形且属于该维第一个
                            node(i).simp{index_max+1}(1) = simp_temp; %第一个更高维的单形，直接添加
                        else                                    %this is the next index_max+1 higher simplex   高维的复形；；
                            no_new_simp = size(node(i).simp{index_max+1}, 2);%若不是第一个，则需要计算下标应该是多少
                            node(i).simp{index_max+1}(no_new_simp+1) = simp_temp;%index_max+1  代表比index_max更高维的复形且属于该维第no_new_simp+1个
                        end
                    end
                end
            end
        end

        if index_max == size(node(i).simp, 2) % there is no higher simplex
            break;
        end
    end
   
end

% %weight权重计算
% for i=1:nmb
%     node(i).weight=0;
%     index_max=size(node(i).simp,2);
%     if index_max>=2&&~isempty(node(i).simp{2}) 
%         edges_index=size(node(i).simp{2},2);
%         for j=1:edges_index
%             if isempty(node(i).simp{2}(j).neighb)
%                 node(i).weight=0;
%                 
%                 break;
%             else if index_max>=3&&~isempty(node(i).simp{3})
%                     tri_index=size(node(i).simp{3},2);
%                     for t=1:tri_index
%                         if isempty(node(i).simp{3}(t).neighb)
%                             node(i).weight=2;
%                             break;
%                         else
%                             node(i).weight=3;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% %计算weight_max
% weight_max=0;
% for i=1:nmb
%     weight(i)=node(i).weight;
% end