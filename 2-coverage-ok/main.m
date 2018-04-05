% This is the main function
clc;
clear;

global node node_coor boundary_edge cycle_time k node_x node_y edge 

% lambda is the intensity 
lambda=140;
k=2;
S1=[];

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

SetNode();

edge.u=[];
edge.v=[];

for i = 1: length(node_coor) - 24*k 
    if node(i).status==1
        no_neighbor = length(node(i).neighbour);
        if no_neighbor 
            for j = 1: no_neighbor
                if node(i).neighbour(j) < i || node(node(i).neighbour(j)).fence_flag == 1  % u>v
                    continue;
                elseif node(node(i).neighbour(j)).status==1
                    edge.u = [edge.u,i];
                    edge.v=[edge.v,node(i).neighbour(j)];
                end 
            end 
        end 
    end
end
cycle_time=[];
  for i=1:length(edge.u)
     line([node_x(edge.u(i)),node_x(edge.v(i))],[node_y(edge.u(i)),node_y(edge.v(i))],'Color', [0 0 0]);
  end

boundary_edge=[];

FindEdge(1);
%%----------------6/17-----------
index1=length(boundary_edge);
if index1
    disp('1-coverage holes：');
    FindCycle();
else
    disp('No 1-coverage holes exists!');
end

if k>1
    S=zeros(1,length(node_coor)-44*k);  %%
    for i=1:length(node_coor)-44*k
        S(i)=i;
        if node(i).bn==1 
            S1=[S1,i];
        else        
          flag_sleep=SleepNode(i);        
          if flag_sleep==1
            node(i).status=0;
            S1=[S1,i];
          end
       end
    end
    S2=setdiff(S,S1);
    for m=1:length(S1)
        node(S1(m)).status=1;
    end
    for n=1:length(S2)
        node(S2(n)).status=0; 
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
%    ------------------------test----------------------%
    figure;
    Theta=[0:0.005:1]*2*pi;
    Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
    Ycircle=1.0*sin(Theta);
    for i=1:length(node_coor)-24*k         % 画出fence node之外所有节点的Rs      
      if node(i).status==1
        Xc=Xcircle+node_x(i);
        Yc=Ycircle+node_y(i);
        plot(Xc,Yc,'k');
        fill(Xc,Yc,[0.5,0.5,0.5],'facealpha',0.8);  %半透明显示
        plot(node_x(i), node_y(i),'.','markersize',10,'color',[0 0 0]);
        hold on;
      end
    end 
% 
%     ------------------------test----------------------%
    FindEdge(2);  %2-coverage
    index1=length(boundary_edge);
    if index1
        disp('2-coverage holes：');
        FindCycle();
    else 
        disp('No 2-coverage holes exists!');
    end     
end

