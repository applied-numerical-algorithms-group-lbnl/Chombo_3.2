amrID = amr_load('plot.amundsen.2d.hdf5');
thkname = 'thickness'; % name of the ice thickness data
thkrange = [0,4000.0]; %sensible range for thickness data
%read data at the coarsest (level 0 resolution)
level = 0;

%work out the domain corners for  level 0
[ lo,hi ] = amr_query_domain_corners(amrID, level);

interp_order = 0; %0 for piecewise constant interpolation, 1 for linear
[ x0,y0,thk0 ] = amr_read_box_2d( amrID, level, lo, hi, thkname, interp_order  );

hold off;
imagesc(x0,y0,thk0,thkrange); colorbar();
axis image % make the pixel aspect ratio 1:1
set(gca,'ydir','normal'); %put thk(1,1) at the bottom left
hold on; 

thkc = [500.0,1000.0,1500.0,2000.0];

%read data at a finer (level 1 resolution)
level = 1;
lo = [50,50]; hi = [150,150]; %box corners
interp_order = 0; %0 for piecewise constant interpolation, 1 for linear
[ x1,y1,thk1 ] = amr_read_box_2d( amrID, level, lo, hi, thkname, interp_order  );
%imagesc(x1,y1,thk1,thkrange);
%draw a rectange around the high res data
dx = x1(2)- x1(1);
w = max(x1)-min(x1)+dx;
h = max(y1)-min(y1)+dx;
rectangle('Position',[min(x1)-dx/2.0,min(y1)-dx/2.0,w,h]);

time=amr_query_time(amrID);

amr_free(amrID);
