
find_path(KNITRO_INCLUDE_DIR 
 NAMES
  knitro.h
  PATHS
  $ENV{KNITRODIR32}
    PATH_SUFFIXES
    include)

find_library(KNITRO_LIBRARY knitro800 $ENV{KNITRO_LIB_PATH})


