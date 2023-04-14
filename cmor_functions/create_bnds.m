%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create bounding variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bnds=create_bnds(dim_in,varargin)

if isstruct(dim_in)
   dim=dim_in.value;
else
   dim=dim_in;
end
bnds=zeros(2,length(dim));

%Passed interface field (lev, a, b, etc.)
if nargin>1
   dim_interface=varargin{1};
   for i=1:size(bnds,2)
      bnds(1,i)=dim_interface(i);
      bnds(2,i)=dim_interface(i+1);
   end
else
   %Valid dim bounds exist (lat, lon)
   if ~isempty(dim_in.info.valid_max)
      maxbnd=str2num(dim_in.info.valid_max);
      minbnd=str2num(dim_in.info.valid_min);

      if xor(dim(1)==minbnd,dim(end)==maxbnd)
         if dim(1)==minbnd
            loop_start=1;
            loop_end=size(bnds,2)-1;
            bnds(1,size(bnds,2))=dim(end);
            bnds(2,size(bnds,2))=maxbnd;
         else
            loop_start=2;
            loop_end=size(bnds,2);
            bnds(1,1)=minbnd;
            bnds(2,1)=dim(1);
         end
         for j=loop_start:loop_end
            bnds(1,j)=dim(j);
            bnds(2,j)=dim(j+1);
         end
      else
         bnds(1,1)=minbnd;
         bnds(2,end)=maxbnd;

         bnds(2,1)=bnds(1,1)+(dim(2)-dim(1))/2;
         bnds(1,end)=bnds(2,end)-(dim(end)-dim(end-1))/2;
         for j=2:size(bnds,2)-1
            bnds(1,j)=dim(j)-(dim(j)-dim(j-1))/2;
            bnds(2,j)=dim(j)+(dim(j+1)-dim(j))/2;
         end
      end
   else
      for j=2:size(bnds,2)-1
         bnds(1,j)=dim(j)-(dim(j)-dim(j-1))/2;
         bnds(2,j)=dim(j)+(dim(j+1)-dim(j))/2;
      end
      bnds(:,1)=[-(dim(2)-dim(1))/2 (dim(2)-dim(1))/2]+dim(1);
      bnds(:,end)=[-(dim(end)-dim(end-1))/2 (dim(end)-dim(end-1))/2]+dim(end);
   end
end
