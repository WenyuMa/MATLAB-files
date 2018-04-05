function[error,hole_num]=main()
%clc;
clear;

global k node1 node2 node_coor node_x node_y edge 

% lambda is the intensity 
lambda=130;
k=2;

nmb=poissrnd(lambda);

x = rand(1,nmb);
y = rand(1,nmb);
x = 10*x + 1;
y = 10*y + 1;

%coordinates of fence nodes
fence_x1 = [0,0,0,0,0, 0, 0,2, 2,4, 4,6, 6,8, 8,10,10,12,12,12,12,12,12,12];
fence_y1 = [0,2,4,6,8,10,12,0,12,0,12,0,12,0,12, 0,12, 0, 2, 4, 6, 8,10,12];

fence_x2 = [0,0,0,0,0, 0,1, 1,3, 3,5, 5,7, 7,9, 9,11,11,12,12,12,12,12,12];
fence_y2 = [1,3,5,7,9,11,0,12,0,12,0,12,0,12,0,12, 0,12, 1, 3, 5, 7, 9,11];

inner_x1 = [1,1,1,1,1, 1,3, 3,5, 5,7, 7,9, 9,11,11,11,11,11,11];
inner_y1 = [1,3,5,7,9,11,1,11,1,11,1,11,1,11, 1, 3, 5, 7, 9,11];

inner_x2 = [1,1,1,1,1, 2, 2,4, 4,6, 6,8, 8,10,10,11,11,11,11,11];
inner_y2 = [2,4,6,8,10,1,11,1,11,1,11,1,11, 1,11, 2, 4, 6, 8,10];

%coordinates of all the nodes
node_x = [x, inner_x1, inner_x2, fence_x1, fence_x2];
node_y = [y, inner_y1, inner_y2, fence_y1, fence_y2];


node_coor = [node_x; node_y];

plot(node_x, node_y, '.'); hold on;


% figure
Theta=[0:0.005:1]*2*pi;
Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
Ycircle=1.0*sin(Theta);
for i=1:length(node_coor)-24*k         % 画出fence node之外所有节点的Rs
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
for i = 1: length(node_coor)-24*k  
    text(node_x(i)+0.1, node_y(i), num2str(i));
end

% set the fence_flag for all nodes to indicate whether they are fence node
% or not
for i = 1: length(node_coor) 
    node1(i).status = 1;   % sleep or active
    node2(i).status = 1;   % sleep or active
    node1(i).bn=0;         % a boundary node or not
    node2(i).bn=0;         % a boundary node or not
    if i > length(node_coor) - 24*k 
        node1(i).fence_flag = 1;
        node2(i).fence_flag = 1;
    else
        node1(i).fence_flag = 0;    
        node2(i).fence_flag = 0;
    end
end

%--------------------------------------1.相邻节点间的信息计算
for i=1: length(node_coor)         
    node1(i).neighbour = [];
    node1(i).dist = [];
    node1(i).angle = [];
    node1(i).msg=[];

    for j=1: length(node_coor)
        if (j==i) 
            continue;
        else
            distance = dist2(node_coor(:,i),node_coor(:,j));
            if distance <= 2 
                node1(i).neighbour = [node1(i).neighbour, j];
                node1(i).dist = [node1(i).dist, distance];
                node1(i).angle = [node1(i).angle, 2*acos(distance/2)];
            end
        end
    end
end

edge.u=[];
edge.v=[];


for i = 1: length(node_coor) - 24*k 
    if node1(i).status==1
        no_neighbor = length(node1(i).neighbour);
        if no_neighbor 
            for j = 1: no_neighbor
                if node1(i).neighbour(j) < i || node1(node1(i).neighbour(j)).fence_flag == 1  % u>v
                    continue;
                elseif node1(node1(i).neighbour(j)).status==1
                    edge.u = [edge.u,i];
                    edge.v=[edge.v,node1(i).neighbour(j)];
                end 
            end 
        end 
    end
end

%---------------------------------RBA-------------------------------------
boundary1();
S1=[];
if k>1
    S=zeros(1,length(node_coor)-44*k);  %%
    for i=1:length(node_coor)-44*k
        S(i)=i;
        if node1(i).bn==1 
            S1=[S1,i];
        else        
          flag_sleep=SleepNode1(i);        
          if flag_sleep==1
            node1(i).status=0;
            S1=[S1,i];
          end
       end
    end
    S2=setdiff(S,S1);
    for m=1:length(S1)
        node1(S1(m)).status=1;
    end
    for n=1:length(S2)
        node1(S2(n)).status=0; 
        index1=find(edge.u==S2(n));
        index2=find(edge.v==S2(n));
        if index1
          for j=1:length(index1)
            edge.u(index1(j))=0;
            edge.v(index1(j))=0;
          end
        end
        if index2
          for j=1:length(index2)
            edge.u(index2(j))=0;
            edge.v(index2(j))=0;
          end         
        end 
    end 
    
%     %    ------------------------test----------------------%
%     figure;
%     Theta=[0:0.005:1]*2*pi;
%     Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
%     Ycircle=1.0*sin(Theta);
%     for i=1:length(node_coor)-24*k         % 画出fence node之外所有节点的Rs      
%     if node1(i).status==1
%         Xc=Xcircle+node_x(i);
%         Yc=Ycircle+node_y(i);
%         plot(Xc,Yc,'k');
%         fill(Xc,Yc,'g','facealpha',1);  %半透明显示
%         text(node_x(i)+0.1, node_y(i), num2str(i));
%         hold on;
%     end
%     end 
% % 
% %     ------------------------test----------------------%
% 
    boundary1();
    result1=find_cycle1();   
end

%%---------------------------------LBA-------------------------------------
boundary2();  %参考算法结果
find_cycle2();
S3=[];
if k>1
    S=zeros(1,length(node_coor)-44*k);  %%
    for i=1:length(node_coor)-44*k
        S(i)=i;
        if node2(i).bn==1 
            S3=[S3,i];
        else        
          flag_sleep=SleepNode2(i);        
          if flag_sleep==1
            node2(i).status=0;
            S3=[S3,i];
          end
       end
    end
    S4=setdiff(S,S3);
    for m=1:length(S3)
        node2(S3(m)).status=1;
    end
    for n=1:length(S4)
        node2(S4(n)).status=0; 
    end
%         %    ------------------------test----------------------%
%     figure;
%     Theta=[0:0.005:1]*2*pi;
%     Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
%     Ycircle=1.0*sin(Theta);
%     for i=1:length(node_coor)-24*k         % 画出fence node之外所有节点的Rs      
%     if node2(i).status==1
%         Xc=Xcircle+node_x(i);
%         Yc=Ycircle+node_y(i);
%         plot(Xc,Yc,'k');
%         fill(Xc,Yc,'g','facealpha',1);  %半透明显示
%         text(node_x(i)+0.1, node_y(i), num2str(i));
%         hold on;
%     end
%     end 
% 
%     ------------------------test----------------------%
    boundary2();
    result2=find_cycle2();   
end


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
            cycle_flag1(j)=1;   % 匹配成功，置1
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


