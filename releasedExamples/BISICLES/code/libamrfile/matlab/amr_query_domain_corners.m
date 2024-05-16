function [ lo,hi ] = amr_query_domain_corners( amrID, level )
%extract the domain corners for a given level grid
load_libamrfile();

%Check that level is valid - matlab crashes if invalid
[nlevel]=amr_query_n_level(amrID);
if level > (nlevel-1)
    level=nlevel-1;
    error(['Input level too high. Max value is ' num2str(level)]);
end

status = libpointer('int32Ptr',-1);
levelp=libpointer('int32Ptr',level);
lop=libpointer('int32Ptr',[0,0]);
hip=libpointer('int32Ptr',[1,1]);
calllib('libamrfile','amr_query_domain_corners',status,lop,hip,amrID,levelp);
amr_error(status);
lo = lop.value;
hi = hip.value;
end

