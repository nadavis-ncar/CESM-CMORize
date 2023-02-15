%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Expand surface pressure to full pressure field with hybrid coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pres=calculate_pressure(ps,a,b,varargin)

disp('we have arrived')

if nargin>3
   column=varargin{1};
else
   column=0;
end

if column==0
size(ps)
   ps_mat=permute(repmat(ps,[1 1 1 length(a)]),[1 2 4 3]);
 size(ps_mat)
   a=permute(repmat(a(:),[1 size(ps)]),[2 3 1 4]);
   b=permute(repmat(b(:),[1 size(ps)]),[2 3 1 4]);
 size(a)
size(b)
size(ps_mat)
   pres=a+b.*ps_mat;
else
   pres=a(:)+b(:)*ps;
end
