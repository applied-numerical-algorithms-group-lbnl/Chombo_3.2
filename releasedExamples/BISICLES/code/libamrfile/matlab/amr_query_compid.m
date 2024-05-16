function [ compID ] = amr_query_compid( amrID, compname )
%lookup the component number (ID) given a name
load_libamrfile();
status = libpointer('int32Ptr',-1);
namelen=libpointer('int32Ptr',length(compname));
compID=libpointer('int32Ptr',0);
calllib('libamrfile','amr_query_comp_id',status,compID,amrID,compname,namelen);
amr_error(status);
end

