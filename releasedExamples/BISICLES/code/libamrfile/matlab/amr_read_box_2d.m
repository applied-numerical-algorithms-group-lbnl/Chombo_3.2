function [ x,y,v ] = amr_read_box_2d( amrID, level, lo, hi, comp, interp_order )
%read data from a Chombo AMR file (amrID) over a uniform mesh.
%Defines a cell centred box with corner cells at lo and hi, 
%on AMR level amr_level, then fills an array with data that
%is copied, averaged, or interploated as needed
load_libamrfile();
compid = libpointer('int32Ptr',-1);
if (ischar(comp))
    compid = amr_query_compid(amrID, comp);
elseif (isinteger(comp))
    compid.value = comp;
end

lop=libpointer('int32Ptr',lo);
hip=libpointer('int32Ptr',hi);
levelp=libpointer('int32Ptr',level);
nx = hi(1)-lo(1)+1; ny = hi(2)-lo(2)+1;
x = zeros(1,nx); xp=libpointer('doublePtr',x);
y = zeros(1,ny); yp=libpointer('doublePtr',y);
v = zeros(nx,ny); vp=libpointer('doublePtr',v);
interp=libpointer('int32Ptr',interp_order);
status = libpointer('int32Ptr',-1);
calllib('libamrfile','amr_read_box_data_2d',status,vp,xp,yp,amrID,levelp,lop,hip,compid,interp);
amr_error(status);
x = xp.value;
y = yp.value;
v = (vp.value)';

end

