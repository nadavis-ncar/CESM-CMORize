%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Enforces dimension rules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var_out=stored_direction(var_out)

for i=1:length(var_out.dim)
   if length(var_out.dim{i}.out.value)>1
      flip_dim=0;
      dim_increasing=(var_out.dim{i}.out.value(2)-var_out.dim{i}.out.value(1))>0;
      if dim_increasing==1
         if strcmp(var_out.dim{i}.out.info.stored_direction,'decreasing')
            flip_dim=1;
         end
      else
         if strcmp(var_out.dim{i}.out.info.stored_direction,'increasing')
            flip_dim=1;
         end
      end
      if flip_dim==1
         var_out.native.value=flip(var_out.native.value,i);
         var_out.dim{i}.out.value=flip(var_out.dim{i}.out.value,i);
         var_out.dim{i}.flip=1;
      else
         var_out.dim{i}.flip=0;
      end
   end
end
