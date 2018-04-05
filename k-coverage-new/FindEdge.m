% This function is used for finding boundary edges

function []= FindEdge()
% all the edges

global k node node_coor node_x node_y boundary_edge cycle_time

edge.u=[];
edge.v=[];
angle=[];
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
                    angle = [angle,node(i).angle(j)];
                end 
            end 
        end 
    end
end
%  for i=1:length(edge.u)
%     line([node_x(edge.u(i)),node_x(edge.v(i))],[node_y(edge.u(i)),node_y(edge.v(i))],'Color', [0.5 0.5 0.5]);
%  end
 
cycle_time=[];
boundary_edge=[];

%find public neighbour
for m=1:length(edge.u)
    
    flag_cover=0;
    edge_weight=0;
    S1=[];
    index1=0;
    
    for i=1:length(node(edge.u(m)).neighbour)
        x=node(edge.u(m)).neighbour(i);
        if node(x).status==1
            index1=find(node(edge.v(m)).neighbour==x);  %公共相邻节点在v邻居中的位置
            if index1
                %y=node(edge.v(m)).neighbour(index1);  % 具体是哪一个节点w
                index2=i;
                index3=find(node(edge.u(m)).neighbour==edge.v(m));  %v在u邻居中的位置
                uv=node(edge.u(m)).dist(index3); %uv 
                uw=node(edge.u(m)).dist(index2); %uw
                wv=node(edge.v(m)).dist(index1); %wv
                wuv=acos((uv^2+uw^2-wv^2)/(2*uv*uw));
                
                angle1=node(edge.u(m)).angle(index3);%uv
                angle2=node(edge.u(m)).angle(index2);%uw
                
                if angle2>=angle1+2*wuv
                      flag_cover=1;
                      edge_weight=0;
                      break;      
                elseif wuv > abs(angle2 - angle1)/2 && wuv <= (angle2 + angle1)/2
                      S1=[S1,x];
                end 
            end 
        end
    end
    
    if flag_cover==0 
        l=length(S1);
           if l==0
              edge_weight=2;
           elseif l==1
              edge_weight=1;
           else
              flag_pair=0;
              for f=1:l-1
                  if flag_pair==0
                     for g=f+1:l 
                         if flag_pair==1
                             break;
                         else
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
                                 
                                 if abs(wux - vux - wuv) < 1e-5 || abs(2*pi-wux - vux - wuv) < 1e-5
                                    flag_pair=1;
                                 end 
                             else
                                 flag_pair=1; 
                             end 
                         end 
                     end
                  else
                      break;
                  end
              end
              
              if flag_pair==1
                  edge_weight=0;
              elseif flag_pair==0
                  edge_weight=1;
              end
           end
    end
    if edge.u(m)>length(node_coor)-44*k && edge.v(m)>length(node_coor)-44*k && edge_weight>1
            edge_weight=1;
    end
    if edge_weight>0
       cycle_time=[cycle_time,edge_weight];
       boundary_edge=[boundary_edge,edge.u(m),edge.v(m)];
    end;
       if edge_weight>0
           line([node_x(edge.u(m)),node_x(edge.v(m))],[node_y(edge.u(m)),node_y(edge.v(m))], 'Color', 'r', 'LineWidth', 3); 
       end           
end

