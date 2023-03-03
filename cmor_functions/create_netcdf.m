%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create NetCDF output file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_netcdf(var_out,filename,globals,varargin)

lev_ind=[];
do_ps=0;
if nargin>3
   ps=varargin{1};
   ps_dim=[];
   do_ps=1;
end

if exist(filename,'file')==2
   delete(filename);
   disp('Overwriting existing file')
end
ncid=netcdf.create(filename,'NETCDF4');
varid = netcdf.getConstant('GLOBAL');
for i=1:length(globals)
   netcdf.putAtt(ncid,varid,globals(i).name,globals(i).value);
end

for i=1:length(var_out.dim)
   ndim(i)=length(var_out.dim{i}.out.value);
   if strcmp(var_out.dim{i}.out.name,'time')
      dim(i)=netcdf.defDim(ncid,var_out.dim{i}.out.info.out_name,netcdf.getConstant('NC_UNLIMITED'));
   else
      dim(i)=netcdf.defDim(ncid,var_out.dim{i}.out.info.out_name,ndim(i));
   end
   if do_ps==1
      if ~strcmp(var_out.dim{i}.out.name,'lev')
         ps_dim=cat(1,ps_dim,dim(i));
      else
         lev_ind=i;
      end
   end
end
dim(i+1)=netcdf.defDim(ncid,'bnds',2);

var=netcdf.defVar(ncid,var_out.info.out_name,'NC_FLOAT',dim(1:end-1));
attnames=fieldnames(var_out.info);
atts=struct2cell(var_out.info);
for i=1:length(atts)
   if ~isempty(atts{i}) & sum(strcmp(attnames{i},{'frequency';'realm';'dimensions';'type';'out_name'}))==0
      netcdf.putAtt(ncid,var,attnames{i},atts{i});
   end
end

if do_ps==1
   var_ps=netcdf.defVar(ncid,'ps','NC_FLOAT',ps_dim);
   attnames=fieldnames(ps.info);
   atts=struct2cell(ps.info);
   for i=1:length(atts)
      if ~isempty(atts{i}) & sum(strcmp(attnames{i},{'frequency';'realm';'dimensions';'type';'out_name'}))==0
         netcdf.putAtt(ncid,var_ps,attnames{i},atts{i});
      end
   end
   formula_names=fieldnames(var_out.dim{lev_ind}.out.formula);
   formula=struct2cell(var_out.dim{lev_ind}.out.formula);
   for i=1:length(formula_names)
      if sum(size(formula{i}.value)==1)==0
         formula_dims{i}=[dim(end) dim(lev_ind)];
      else
         formula_dims{i}=dim(lev_ind);
      end
      var_formula(i)=netcdf.defVar(ncid,formula_names{i},'NC_FLOAT',formula_dims{i});
      attnames=fieldnames(formula{i}.info{1});
      atts=struct2cell(formula{i}.info{1});
      for j=1:length(atts)
         if ~isempty(atts{j}) & sum(strcmp(attnames{j},{'dimensions';'type';'out_name'}))==0
            netcdf.putAtt(ncid,var_formula(i),attnames{j},atts{j})
         end;
      end
   end
end
for i=1:length(var_out.dim)
   var_dim(i)=netcdf.defVar(ncid,var_out.dim{i}.out.info.out_name,'NC_FLOAT',dim(i));
   if length(var_out.dim{i}.out.value)>1
      var_dim_bnds(i)=netcdf.defVar(ncid,[var_out.dim{i}.out.info.out_name,'_bnds'],'NC_FLOAT',[dim(end) dim(i)]);
   end
   attnames=fieldnames(var_out.dim{i}.out.info);
   atts=struct2cell(var_out.dim{i}.out.info);
   for j=1:length(atts)
       if ~isempty(atts{j}) & sum(strcmp(attnames{j},{'out_name';'must_have_bounds';'stored_direction';'type';'requested'}))==0
          netcdf.putAtt(ncid,var_dim(i),attnames{j},atts{j});
          if length(var_out.dim{i}.out.value)>1
             netcdf.putAtt(ncid,var_dim_bnds(i),attnames{j},atts{j});
          end
       end
   end
end

netcdf.endDef(ncid);
start=zeros(length(var_out.dim),1);
for i=1:length(var_out.dim)
   ending(i)=length(var_out.dim{i}.out.value);
end
netcdf.putVar(ncid,var,start,ending,var_out.native.value)
if do_ps==1
   netcdf.putVar(ncid,var_ps,ps.value);
   for i=1:length(var_formula)
      netcdf.putVar(ncid,var_formula(i),formula{i}.value);
   end
end
for i=1:length(var_out.dim)
   netcdf.putVar(ncid,var_dim(i),var_out.dim{i}.out.value)
   if length(var_out.dim{i}.out.value)>1
      netcdf.putVar(ncid,var_dim_bnds(i),var_out.dim{i}.out.bnds);
   end
end
netcdf.close(ncid);
