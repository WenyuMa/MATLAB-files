% This is the main function
clc;
clear;

global node node_coor boundary_edge k node_x node_y 

% lambda is the intensity 
lambda=110;
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
FindEdge();

index1=length(boundary_edge);
if index1
    disp('1-coverage holes：');
    FindCycle();
elseif k>1 %
    S=zeros(1,length(node_coor)-44*k);
    for i=1:length(node_coor)-44*k
        S(i)=i;
        flag_sleep=SleepNode(i);
        if flag_sleep==1
            node(i).status=0;
            S1=[S1,i];
        end
    end
    S2=setdiff(S,S1);
    for m=1:length(S1)
        node(S1(m)).status=1;
    end
    for n=1:length(S2)
        node(S2(n)).status=0;
    end
    %------------------------test----------------------%
    % figure
    Theta=[0:0.005:1]*2*pi;
    Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
    Ycircle=1.0*sin(Theta);
    for i=1:length(node_coor)-24*k         % 画出fence node之外所有节点的Rs      
    if node(i).status==1
        Xc=Xcircle+node_x(i);
        Yc=Ycircle+node_y(i);
        plot(Xc,Yc,'k');
        fill(Xc,Yc,'g','facealpha',0.7);  %半透明显示
   
        hold on;
    end
    end 

    %------------------------test----------------------%
    
    
    FindEdge();
    index2=length(boundary_edge);
    if index2
        disp('2-coverage holes：');
        FindCycle();
    end  
end