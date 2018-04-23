function [flag]=EdgeCollapse()

global node node_coor edge node_x node_y
flag=0;
sleep_edge=[];

%% edge collapse
for i=1:length(edge.u)
    if node(edge.u(i)).status==1 && node(edge.v(i)).status==1 && (edge.u(i)<=length(node_coor-28))
        u=edge.u(i);
        v=edge.v(i);
%         if node(u).c_flag==0 && node(v).c_flag==0
%             edge.c_flag(i)=0;
        if edge.c_flag(i)==1
            co_neighb1=intersect(node(u).neighbor,node(v).neighbor);
            co_neighb1=setdiff(co_neighb1,[u,v]);  %两条边一定不能有公共节点吗？

            if length(co_neighb1)>1
                sleep_temp=[];
                temp_index=[];
                for m=1:length(co_neighb1)-1
                    w=co_neighb1(m);
                    if node(w).fence_flag==0  % exclude edges with fence nodes
                        index1=find(edge.u==w);
                        if index1
                            for n=m+1:length(co_neighb1)
                                x=co_neighb1(n);
                                index2=find(edge.v(index1)==x);
                                if index2   %edge wx locates in edge.u(index1(index2))
                                    if node(x).fence_flag==0
                                        co_neighb2=intersect(node(w).neighbor,node(x).neighbor);
                                        co_neighb3=intersect(co_neighb1,co_neighb2);
                                        diff=setdiff(co_neighb1,co_neighb3);  %belong to 1 but not 3
                                        if isempty(diff)
                                            % ADD SOME JUDGE HERE                                             
                                            edge_index=EdgeCheck(i,index1(index2));
                                            sleep_temp=[sleep_temp,index1(index2)];
                                            temp_index=[temp_index,edge_index];
                                        end                                      
                                    end
                                end
                            end
                        end
                    end
                end
                if sleep_temp
                    m=min(temp_index);
                    min_index=find(temp_index==m);
                    if m<3
                        node(u).neighbor=setdiff(node(u).neighbor,v);
                        node(v).neighbor=setdiff(node(v).neighbor,u);
                        sleep_edge=[sleep_edge,i];
                        w=edge.u(sleep_temp(min_index(1)));
                        x=edge.v(sleep_temp(min_index(1)));
                        edge.c_flag(sleep_temp(min_index(1)))=0;
                        disp(['edge (' num2str(u) ',' num2str(v) ')' ' is dominated by edge (' num2str(w) ',' num2str(x) ')']);
                        %line([node_x(u),node_x(v)],[node_y(u),node_y(v)], 'Color', 'r', 'linewidth', 1);
                        %line([node_x(w),node_x(x)],[node_y(w),node_y(x)], 'Color', 'b', 'linewidth', 1);
                    end
                end
            end
        end
    end
end

if ~isempty(sleep_edge)
    for i=1:length(sleep_edge)
        edge.u(sleep_edge(i))=0;
        edge.v(sleep_edge(i))=0;
        edge.c_flag(sleep_edge(i))=-1;
    end
    edge.u=edge.u(edge.u~=0);
    edge.v=edge.v(edge.v~=0);
    edge.c_flag=edge.c_flag(edge.c_flag~=-1);
    flag=1;
end
