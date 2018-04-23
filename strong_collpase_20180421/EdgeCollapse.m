function [flag]=EdgeCollapse()

global node node_coor node_x node_y
flag=0;
%% edge collapse
for i=1:length(node_coor)-28
    if node(i).status==1
        for j=1:length(node(i).simp{1}.neighb)
            if i<node(i).simp{1}.neighb(j)
                v=node(i).simp{1}.neighb(j);
                co_neighb1=intersect(node(i).neighbor,node(v).neighbor);
                if length(node(i).simp{2}(j).neighb)>1
                    neighb_temp=node(i).simp{2}(j).neighb;
                    for m=1:length(neighb_temp)-1
                        w=neighb_temp(m);
                        for n=m+1:length(neighb_temp)
                            x=neighb_temp(n);
                            co_neighb2=intersect(node(w).neighbor,node(x).neighbor);
                            co_neighb3=intersect(co_neighb1,co_neighb2);
                            diff=setdiff(co_neighb1,co_neighb3);  %belong to 1 but not 3
                            if isempty(diff)
                                flag=1;
                                disp(['edge (' num2str(i) ',' num2str(v) ')' ' is dominated by edge (' num2str(w) ',' num2str(x) ')']);
                                line([node_x(i),node_x(v)],[node_y(i),node_y(v)], 'Color', 'r', 'linewidth', 1);
                                node(i).neighbor=setdiff(node(i).neighbor,v);
                                node(v).neighbor=setdiff(node(v).neighbor,i);
                            end                            
                        end
                    end
                end
            end           
        end
    end
end