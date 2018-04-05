function [hole]=boundary()
%设置Rs=1，Rc=2*Rs，节点在12*12的正方形区域内分布
clc;
clear;

% lambda is the intensity 
lambda=200;
k=2;
nmb=poissrnd(lambda);


x = rand(1,nmb);
y = rand(1,nmb);
x = 10*x + 1;
y = 10*y + 1;

%coordinates of fence nodes
fence_x = [0,0,0,0,0, 0, 0,2, 2,4, 4,6, 6,8, 8,10,10,12,12,12,12,12,12,12];
fence_y = [0,2,4,6,8,10,12,0,12,0,12,0,12,0,12, 0,12, 0, 2, 4, 6, 8,10,12];

fence_x2 = [0,0,0,0,0, 0,1, 1,3, 3,5, 5,7, 7,9, 9,11,11,12,12,12,12,12,12];
fence_y2 = [1,3,5,7,9,11,0,12,0,12,0,12,0,12,0,12, 0,12, 1, 3, 5, 7, 9,11];

inner_x = [1,1,1,1,1, 1,3, 3,5, 5,7, 7,9, 9,11,11,11,11,11,11];
inner_y = [1,3,5,7,9,11,1,11,1,11,1,11,1,11, 1, 3, 5, 7, 9,11];

inner_x2 = [1,1,1,1,1, 2, 2,4, 4,6, 6,8, 8,10,10,11,11,11,11,11];
inner_y2 = [2,4,6,8,10,1,11,1,11,1,11,1,11, 1,11, 2, 4, 6, 8,10];

%coordinates of all the nodes
% node_x = [x, inner_x, fence_x];
% node_y = [y, inner_y, fence_y];

node_x = [x, inner_x, inner_x2, fence_x, fence_x2];
node_y = [y, inner_y, inner_y2, fence_y, fence_y2];


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
    xlim([1 11]); %x,y轴上下限设置
    ylim([1 11]);
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
    if i > length(node_coor) - 24*k 
        node(i).fence_flag = 1;
    else
        node(i).fence_flag = 0;
    end
end


%相邻节点间的信息计算
for i=1: length(node_coor)-24  %*k  
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

% all the edges
edge.u=[];
edge.v=[];
angle=[];
for i = 1: length(node_coor) - 24*k 
    no_neighbor = length(node(i).neighbour);
    if no_neighbor 
        for j = 1: no_neighbor
            if node(i).neighbour(j) < i || node(node(i).neighbour(j)).fence_flag == 1  % u>v
                continue;
            else
                edge.u = [edge.u,i];
                edge.v=[edge.v,node(i).neighbour(j)];
                angle = [angle,node(i).angle(j)];
            end
        end
    end
end
 for i=1:length(edge.u)
   % line([node_x(edge.u(i)),node_x(edge.v(i))],[node_y(edge.u(i)),node_y(edge.v(i))],'Color', [0.5 0.5 0.5]);
 end
    
cycle_time=[];
boundary_edge=[];
edge_cover=[];

%find public neighbour
for m=1:length(edge.u)
    
    flag_cover=0;
    flag_stop=0;
    edge_weight=0;
    S1=[];
    index1=0;
    
    for i=1:length(node(edge.u(m)).neighbour)
        x=node(edge.u(m)).neighbour(i);
        index1=find(node(edge.v(m)).neighbour==x);  %公共相邻节点在v邻居中的位置 

        if index1
            y=node(edge.v(m)).neighbour(index1);  % 具体是哪一个节点w
            index2=i;
            index3=find(node(edge.u(m)).neighbour==edge.v(m));  %v在u邻居中的位置
            uv=node(edge.u(m)).dist(index3); %uv 
            uw=node(edge.u(m)).dist(index2); %uw
            wv=node(edge.v(m)).dist(index1); %wv
            wuv=acos((uv^2+uw^2-wv^2)/(2*uv*uw));
            
            angle1=node(edge.u(m)).angle(index3);%uv
            angle2=node(edge.u(m)).angle(index2);%uw

            if angle2>=angle1+2*wuv
               flag_cover=flag_cover+1;
            elseif  wuv > abs(angle2 - angle1)/2 && wuv <= (angle2 + angle1)/2
               S1=[S1,y];
            end        
        end  
    end  
%%------------------------------------  
    if flag_cover>=k 
       edge_weight=0;
    else 
        l=length(S1);
        if l<k- flag_cover
           edge_weight=2;
        else 
            if flag_cover==k-1 && l==1
               edge_weight=1;
            else
                flag_pair=zeros(1,l);

                for f=1:l-1
                    if flag_stop==0
                        if flag_pair(f)==0
                            for g=f+1:l
                                if  flag_pair(f)==0 && flag_stop==0
                                    if flag_pair(g)==0
                                      %  if g<=f
                                       %    continue; 
                                       % else 
                                            neibr=find(node(S1(f)).neighbour==S1(g));
                                            if neibr
                                                index4=find(node(edge.u(m)).neighbour==S1(g));%ux
                                                index5=find(node(edge.v(m)).neighbour==S1(g));%vx
                                                index6=find(node(edge.u(m)).neighbour==S1(f));%uw
                                                index7=find(node(edge.v(m)).neighbour==S1(f));%vw
                                                index8=neibr;  %wx
                                 
                                                ux=node(edge.u(m)).dist(index4); %ux 
                                                vx=node(edge.v(m)).dist(index5); %vx
                                                uw=node(edge.u(m)).dist(index6); %uw 
                                                wv=node(edge.v(m)).dist(index7); %wv
                                                wx=node(S1(f)).dist(index8);%wx
                                 
                                                wux=acos((uw^2+ux^2-wx^2)/(2*uw*ux));                      
                                                vux=acos((uv^2+ux^2-vx^2)/(2*uv*ux));
                                                wuv=acos((uw^2+uv^2-wv^2)/(2*uw*uv));
                                                if  abs(wux - vux - wuv) < 1e-5 || abs(2*pi-wux - vux - wuv) < 1e-5
                                                    flag_pair(f)=1;
                                                    flag_pair(g)=1;
                                                    flag_cover = flag_cover+1;
                                                end                                                                                      
                                            else                                    
                                                     flag_pair(f)=1;
                                                     flag_pair(g)=1; 
                                                     flag_cover = flag_cover+1;
                                            end
                                        %end  
                                    end
                                else
                                    break;
                                end
                                if flag_cover==k
                                   flag_stop=1;
                                   edge_weight=0;
                                   break;
                                end  
                            end
                        else
                            continue;
                        end
                    else
                        break;
                    end
                end
                remain=sum(flag_pair(:)==0);
                if flag_cover<k
                    if flag_cover+remain>=k
                       edge_weight=1;
                    else
                       edge_weight=2;
                    end
                end                          
            end
        end
    end
    
    if edge.u(m)>length(node_coor)-44*k && edge.v(m)>length(node_coor)-44*k && edge_weight>1   %
            edge_weight=1;
    end
    if edge_weight>0
       cycle_time=[cycle_time,edge_weight];
       boundary_edge=[boundary_edge,edge.u(m),edge.v(m)];
       edge_cover=[edge_cover,flag_cover];
    end;
       if edge_weight>0
         if flag_cover==0
            line([node_x(edge.u(m)),node_x(edge.v(m))],[node_y(edge.u(m)),node_y(edge.v(m))], 'Color', 'r', 'LineWidth', 3); 
         elseif flag_cover==1
            line([node_x(edge.u(m)),node_x(edge.v(m))],[node_y(edge.u(m)),node_y(edge.v(m))], 'Color', 'b', 'LineWidth',2); 
         end 
       end           
end
%% delete isolated edges
boundary_edge0=unique(boundary_edge);
for i=1:length(boundary_edge0)
   c=sum(boundary_edge==boundary_edge0(i));
   c1=0;
   if c>1
       continue;
   else
       i1=find(boundary_edge==boundary_edge0(i));
       j=ceil(i1/2);
       cycle_time(j)=0;
       if j==i1/2
          i1=i1-1; 
       elseif  2*j-1==i1
          i1=i1+1;
       end
       index=find(boundary_edge==boundary_edge(i1));
       if length(index)<3
              for i2=1:length(index)
                  j1=ceil(index(i2)/2);
                  if cycle_time(j1)>0
                     c1=c1+1;
                     break;
                  end  
              end  
              while(c1==1)
                   c1=0;
                   cycle_time(j1)=0;
                   if j1==index(i2)/2
                       i1=index(i2)-1;
                   elseif 2*j1-1==index(i2)
                       i1=index(i2)+1;
                   end   
                   index=find(boundary_edge==boundary_edge(i1));
                   if length(index)<3
                      for i2=1:length(index)
                          j1=ceil(index(i2)/2);
                          if cycle_time(j1)>0
                             c1=c1+1;
                             break; 
                          end     
                      end    
                   end  
              end  
       end 
   end 
end

%% boundary cycle dectect

%-----------------find initial edge----------------------%
for hole_type=1:k
   disp([num2str(hole_type),'-coverage holes：']);
   x=find(edge_cover==hole_type-1);  
   if x
       t=length(x);
       cycle_time_temp=zeros(1,t);
       boundary_edge_temp=zeros(1,2*t);
       for h=1:t
           cycle_time_temp(h)=cycle_time(x(h));         %all k-cover edges and their weights
           boundary_edge_temp(2*h-1)=boundary_edge(2*x(h)-1);
           boundary_edge_temp(2*h)=boundary_edge(2*x(h));
       end
   else
       continue;
   end
   next=[];
   last=[];
   ini_edge=[0,0];
   over=find(cycle_time_temp>0);
   flag1=length(over); 
   
   % find initial edge
  while flag1
    ini_old=ini_edge;
    k1=0;
    k2=0;
    k3=0;
    a=0;
    x=find(cycle_time_temp==2);
    if x
       for i=1:length(x)
          index1=find(boundary_edge_temp==boundary_edge_temp(2*x(i)-1));
          index2=find(boundary_edge_temp==boundary_edge_temp(2*x(i)));
          
          %find weight 2 edge that connects two other edges of weight 1
          if length(index1)==3 
             for i1=1:3
                 if index1(i1)==2*x(i)-1
                    continue;
                 else
                     num1=ceil(index1(i1)/2);
                     if cycle_time_temp(num1)==1
                         a=a+1;
                     end
                 end
             end
             if a==2
                ini_edge=[boundary_edge_temp(2*x(i)-1),boundary_edge_temp(2*x(i))] ;
                ini_weight=cycle_time_temp(x(i));
                k1=1;
                break;
             end
          end
          if length(index2)==3 
              a=0;
              for i1=1:3
                 if index2(i1)==2*x(i)
                    continue;
                 else
                     num1=ceil(index2(i1)/2);
                     if cycle_time_temp(num1)==1
                         a=a+1;
                     end
                 end
              end
              if a==2
                 ini_edge=[boundary_edge_temp(2*x(i)),boundary_edge_temp(2*x(i)-1)] ;
                 ini_weight=cycle_time_temp(x(i));
                 k1=1;
                 break;
              else
                  continue;
              end
          else
              continue;
          end
       end
       % find weight 2 edge that connect 2 other edges of any weight
       if k1==0
              for i=1:length(x)
                  index1=find(boundary_edge_temp==boundary_edge_temp(2*x(i)-1));
                  index2=find(boundary_edge_temp==boundary_edge_temp(2*x(i)));
                  if length(index1)==3 
                      ini_edge=[boundary_edge_temp(2*x(i)-1),boundary_edge_temp(2*x(i))] ;
                      ini_weight=cycle_time_temp(x(i));
                      k2=1;
                      break;
                  elseif length(index2)==3
                      ini_edge=[boundary_edge_temp(2*x(i)),boundary_edge_temp(2*x(i)-1)] ;
                      ini_weight=cycle_time_temp(x(i));
                      k2=1;
                      break;
                  else 
                      continue;
                  end 
              end
              % find any edge of weight 2
              if k2==0
                 number=unidrnd(length(x));
                 y=x(number);
                 ini_edge=[boundary_edge_temp(2*y-1),boundary_edge_temp(2*y)];
                 ini_weight=cycle_time_temp(y);          
              end 
       end
       elseif find(cycle_time_temp==1)  
              x=find(cycle_time_temp==1);
              for i=1:length(x)
                  index1=find(boundary_edge_temp==boundary_edge_temp(2*x(i)-1));
                  index2=find(boundary_edge_temp==boundary_edge_temp(2*x(i)));
                  % find weight 1 edge that connects 2 other edges
                  if length(index1)==2 && length(index2)==2
                     ini_edge=[boundary_edge_temp(2*x(i)-1),boundary_edge_temp(2*x(i))] ;
                     ini_weight=cycle_time_temp(x(i));
                     k3=1;
                     break;
                  else 
                      continue;
                  end 
              end 
              % find any edge of weight 1
              if k3==0
                 number=unidrnd(length(x));
                 y=x(number);
                 ini_edge=[boundary_edge_temp(2*y-1),boundary_edge_temp(2*y)];
                 ini_weight=cycle_time_temp(y);          
              end     
    end 
    % same initial edge as last time, choose another edge by random
    if ini_edge==ini_old
        y=find(cycle_time_temp>0);
        number=unidrnd(length(y));
        x=y(number);
        ini_edge=[boundary_edge_temp(2*x-1),boundary_edge_temp(2*x)];
        ini_weight=cycle_time_temp(x);
    end 
    
   %-------------------------broadcast--------------------------%
    ini=ini_edge(1);
    node(ini).msg=ini;
    next=ini_edge(2);
    start=ini_edge(2);
    node(start).msg=[ini,start];
    flag2=length(next);
    while flag2
       last=next;
       next=[];  %
       for i=1:length(last)
           index=find(boundary_edge_temp==last(i));
           flag_find=0;      %独立边标志位
           for j=1:length(index)
                   no_edge=ceil(index(j)/2);
               if cycle_time_temp(no_edge)>0
                   if mod(index(j),2)==1          %奇数偶数判断
                      next_jump=boundary_edge_temp(index(j)+1);
                   else  
                      next_jump=boundary_edge_temp(index(j)-1); 
                   end
                   l=length(node(last(i)).msg);
                   if l>1 && next_jump==node(last(i)).msg(l-1)  %回到原来的边，continue
                       flag_re=1;
                       flag_ok=0;
                   else
                       flag_find=1;
                       if ini==next_jump    % return to the initial node?
                          flag_ok=1;
                       else  
                          flag_ok=0;
                       end
                       
                       if node(next_jump).msg
                          if flag_ok==0 && node(last(i)).msg(1)== node(next_jump).msg(1) && node(last(i)).msg(2)== node(next_jump).msg(2)  % old message?
                             flag_re=1;
                          elseif flag_ok==1 
                              if l>2 && length(node(next_jump).msg)>1
                                 if node(last(i)).msg(1)== node(next_jump).msg(1) && node(last(i)).msg(2)== node(next_jump).msg(2)&&node(last(i)).msg(3)==node(next_jump).msg(3)
                                    flag_re=1;
                                 else 
                                    flag_re=0;
                                 end
                              else
                                  flag_re=0;
                              end
                          else  
                             flag_re=0;
                          end 
                       else 
                           flag_re=0;
                       end 
                   end
                       if flag_ok==0 && flag_re==0   
                          next=[next,next_jump];
                          node(next_jump).msg=[node(last(i)).msg,next_jump];
                       elseif flag_ok==0 && flag_re==1
                              continue; 
                       elseif flag_ok==1 && flag_re==1
                              continue;
                       elseif flag_ok==1 && ini_weight>0
                           node(next_jump).msg=[node(last(i)).msg,next_jump];   
                           hole=[node(last(i)).msg,next_jump];
                              disp(hole);
                              for n=1:length(hole)-1
                                  le=find(boundary_edge_temp==hole(n));
                                  for m=1:length(le)
                                      if hole(n)>hole(n+1) 
                                         if hole(n+1)==boundary_edge_temp(le(m)-1)
                                            num=le(m)/2;%
                                            break;
                                         else   
                                            continue;
                                         end  
                                      else 
                                          if hole(n+1)==boundary_edge_temp(le(m)+1)
                                             num=(le(m)+1)/2;%
                                             break;
                                          else  
                                             continue;
                                          end  
                                      end 
                                  end
                                  cycle_time_temp(num)=cycle_time_temp(num)-1;
                              end
                              ini_weight=ini_weight-1;
                       end  
               end 
           end
           if ini_weight==0
               next=[];
               break;
           end
           if flag_find==0               %找不到下一条边
               error=node(last(i)).msg;
               for n1=1:length(error)-1
                   le=find(boundary_edge_temp==error(n1));
                   for m1=1:length(le)
                       if error(n1)>error(n1+1) 
                          if error(n1+1)==boundary_edge_temp(le(m1)-1)
                             num2=le(m1)/2;%
                             break;
                          else    
                             continue;
                          end   
                       else  
                           if error(n1+1)==boundary_edge_temp(le(m1)+1)
                              num2=(le(m1)+1)/2;%
                              break;
                           else   
                              continue;
                           end   
                       end   
                   end 
                      cycle_time_temp(num2)=cycle_time_temp(num2)-1; 
               end
           end
       end
       flag2=length(next);
    end
    
   % find next initial edge
   over=find(cycle_time_temp>0);
   flag1=length(over);
  end
end