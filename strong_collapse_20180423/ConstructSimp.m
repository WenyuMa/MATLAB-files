
function []=ConstructSimp()

global node node_coor

for i = 1 : length(node_coor)
    node(i).simp=[]; 
    if node(i).status==1
        new_simp = struct('vert', [i], 'neighb', setdiff(node(i).neighbor,i));
        node(i).simp{1} = new_simp;
    end
end
%然后依次递归 找出每个顶点的k-simplex
for i = 1 : length(node_coor)
    while(1)
      index_max = size(node(i).simp, 2);%最大的单形,返回矩阵的列数

      if index_max>0
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
                        neighb_set_temp = intersect(neighb_set, setdiff(node(new_node).neighbor,new_node));%当前i顶点的邻居集合  和  当前i顶点 第k个邻居顶点的邻居集合  的交集   
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
      end 

        if index_max == size(node(i).simp, 2) % there is no higher simplex
            break;
        end
    end
   
end