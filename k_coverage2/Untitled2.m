%% boundary cycle dectect

%-----------------find initial edge----------------------%
for hole_type=1:k
   x=find(edge_cover==hole_type-1);  
   if x
       t=length(x);
       cycle_time_temp=zeros(1,t);
       boundary_edge_temp=zeros(1,2*t);
       for h=1:t
           cycle_time_temp(h)=cycle_time_temp(x(h));         %all k-cover edges and their weights
           boundary_edge_temp(2*h-1)=boundary_edge_temp(2*x(h)-1);
           boundary_edge_temp(2*h)=boundary_edge_temp(2*x(h));
       end 
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