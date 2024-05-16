function amr_free( amrID )
%free memory holding Chombo AMR data
load_libamrfile();
status = libpointer('int32Ptr',-1);
calllib('libamrfile','amr_free',status,amrID);
amr_error(status)
end

