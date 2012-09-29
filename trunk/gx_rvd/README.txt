This code is basically taken from gx_cvt3d_bis and only contains the parts made
to build, maintain and display a restricted voronoi diagram.These classes
contain :

RVDApp deriving from GeexApp:
  all the necessary GUI associated with a RVD, and handling of the command line
  arguments

RVDObj deriving from GeexOb to let Geex consider it as the scene
  deriving from RVDGeometry and RVDGraphics

RVDGeometry builds the RVD, initialises the samples

RVDGraphics displays the RVD
