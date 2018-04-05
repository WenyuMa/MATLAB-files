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
