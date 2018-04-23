
function []=ConstructSimp()

global node node_coor

for i = 1 : length(node_coor)
    node(i).simp=[]; 
    if node(i).status==1
        new_simp = struct('vert', [i], 'neighb', setdiff(node(i).neighbor,i));
        node(i).simp{1} = new_simp;
    end
end
%Ȼ�����εݹ� �ҳ�ÿ�������k-simplex
for i = 1 : length(node_coor)
    while(1)
      index_max = size(node(i).simp, 2);%���ĵ���,���ؾ��������

      if index_max>0
        no_max_simp = size(node(i).simp{index_max}, 2); %????
        %����ε������Ŀ
%����ÿ����������ά������no_max_simp��   ����ÿ�����ά����
% ����
        for j = 1 : no_max_simp
            vert_set = node(i).simp{index_max}(j).vert;%��j�����ά���εĶ���ļ��� ����j-1���εĶ��㼯��    (��0-simplex��ʼ���㣩
            neighb_set = node(i).simp{index_max}(j).neighb;%��j�����ά���ε��ھӼ���

            if ~isempty(neighb_set)
                
                no_neighb = length(neighb_set);
                for k = 1 : no_neighb     %�����j�����ά���� �ھӼ��ϵ�ÿ������
                    new_node = neighb_set(k);%���ڹ��ɸõ�j�����ε��ھӼ��ϵĵ�k������
                    if new_node < max(setdiff(vert_set, i))  %setdiff(vert_set, i)����vert_set�д��ڶ�i�в����ڵ�Ԫ��
                        continue;
                    else
                       
                        vert_set_temp = union(vert_set, new_node);
                        neighb_set_temp = intersect(neighb_set, setdiff(node(new_node).neighbor,new_node));%��ǰi������ھӼ���  ��  ��ǰi���� ��k���ھӶ�����ھӼ���  �Ľ���   
                        simp_temp = struct('vert', vert_set_temp, 'neighb', neighb_set_temp);%�ɵ�i������͵�k�����㹹��  �Լ����ǹ�ͬ���ھӶ��㹹�ɵ�higher simplex

                        if index_max == size(node(i).simp, 2) %this is the first index_max+1 new higher simplex �����index_max����ά�ĸ��������ڸ�ά��һ��
                            node(i).simp{index_max+1}(1) = simp_temp; %��һ������ά�ĵ��Σ�ֱ�����
                        else                                    %this is the next index_max+1 higher simplex   ��ά�ĸ��Σ���
                            no_new_simp = size(node(i).simp{index_max+1}, 2);%�����ǵ�һ��������Ҫ�����±�Ӧ���Ƕ���
                            node(i).simp{index_max+1}(no_new_simp+1) = simp_temp;%index_max+1  �����index_max����ά�ĸ��������ڸ�ά��no_new_simp+1��
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