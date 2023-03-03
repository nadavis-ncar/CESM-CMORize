%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates TEM output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var_out=calculate_tem(var_out,missing_value)

g=9.81;
a=6371e3;
omega=7.292e-5;
H=6800;

lat=var_out.dim{1}.native.value;
lev=var_out.dim{2}.native.value;
time=var_out.dim{3}.native.value;

lev=repmat(lev(:)',[length(lat) 1]);
rho=(1.5/1e5)*lev;
z=log(lev/1e5)*-H;
lat=repmat(lat(:),[1 size(lev,2)]);
f=2*omega*sind(lat);
cosmat=cosd(lat);
phi=deg2rad(lat);

uzm=var_out.var{1};
vzm=var_out.var{2};
wzm=var_out.var{3};
vthzm=var_out.var{4};
uvzm=var_out.var{5};
uwzm=var_out.var{6};
thzm=var_out.var{7};

var_out.tem.utendepfd=zeros(size(uzm));
var_out.tem.vtem=zeros(size(uzm));
var_out.tem.wtem=zeros(size(uzm));
var_out.tem.epfy=zeros(size(uzm));
var_out.tem.epfz=zeros(size(uzm));

%This operation is so fast, the overhead of parallel is probably longer than just
%doing this with one core unless the time dimension is VERY large
for t=1:length(time)
   [dudz,dummy]=gradient(squeeze(uzm(:,:,t)),squeeze(z(1,:,1)),squeeze(phi(:,1,1)));
   [dummy,dudy]=gradient(squeeze(uzm(:,:,t)).*cosmat,squeeze(z(1,:,1)),squeeze(phi(:,1,1)));
   dudy=dudy.*(1./(a*cosmat));
   [dtheta_dz,dummy]=gradient(squeeze(thzm(:,:,t)),squeeze(z(1,:,1)),squeeze(phi(:,1,1)));

   epfy=rho.*a.*cosmat.*(dudz.*squeeze(vthzm(:,:,t))./dtheta_dz-squeeze(uvzm(:,:,t)));
   epfz=rho.*a.*cosmat.*((f-dudy).*squeeze(vthzm(:,:,t))./dtheta_dz-squeeze(uwzm(:,:,t)));

   [dummy,dFydy]=gradient(epfy.*cosmat,squeeze(z(1,:,1)),squeeze(phi(:,1,1)));
   dFydy=dFydy.*(1./(a*cosmat));
   [dFzdz,dummy]=gradient(epfz,squeeze(z(1,:,1)),squeeze(phi(:,1,1)));
   var_out.tem.utendepfd(:,:,t)=(1./(rho.*a.*cosmat)).*(dFydy+dFzdz);

   %residual circulation
   [dTdz,dummy]=gradient(rho.*squeeze(vthzm(:,:,t))./dtheta_dz,squeeze(z(1,:,1)),squeeze(phi(:,1,1)));
   [dummy,dTdy]=gradient(cosmat.*squeeze(vthzm(:,:,t))./dtheta_dz,squeeze(z(1,:,1)),squeeze(phi(:,1,1)));
   var_out.tem.vtem(:,:,t)=squeeze(vzm(:,:,t))-(dTdz./rho);
   var_out.tem.wtem(:,:,t)=squeeze(wzm(:,:,t))+(1./(a*cosmat)).*dTdy;

   var_out.tem.epfy(:,:,t)=epfy./rho;
   var_out.tem.epfz(:,:,t)=epfz./rho;
end
