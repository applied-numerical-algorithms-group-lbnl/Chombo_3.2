function [ n_level ] = amr_query_n_level( amrID )
%get the level of finest mesh in the amr file
load_libamrfile();
status = libpointer('int32Ptr',-1);
levelp=libpointer('int32Ptr',0);
calllib('libamrfile','amr_query_n_level',status,levelp,amrID);
amr_error(status);
n_level = levelp.value;
end

