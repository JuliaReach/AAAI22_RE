# Try to find the Gurobi libraries and the includes


find_path(GUROBI_header_PATH
          NAMES gurobi_c++.h
          PATHS /home/christian/programs/gurobi-8.0.1/linux64/include
          #/home/ubuntu/gurobi702/linux64/include  /Library/gurobi702/mac64/include
          )



find_library(GUROBI_C_LIBRARY
	     NAMES gurobi80
	     PATHS /home/christian/programs/gurobi-8.0.1/linux64/lib
	     #/Library/gurobi702/mac64 /Library/gurobi702/mac64/lib/ /home/ubuntu/gurobi701/linux64/lib /home/srirams/gurobi701/linux64/lib
	     )

find_library(GUROBI_CPP_LIBRARY
             NAMES gurobi_c++
             PATHS /home/christian/programs/gurobi-8.0.1/linux64/lib
             #/Library/gurobi702/mac64 /Library/gurobi702/mac64/lib/ /home/ubuntu/gurobi701/linux64/lib /home/srirams/gurobi701/linux64/lib
             )

if (GUROBI_header_PATH)
   message(1)
endif(GUROBI_header_PATH)
if ( GUROBI_C_LIBRARY)
   message(2)
endif( GUROBI_C_LIBRARY)
if (GUROBI_CPP_LIBRARY)
   message(3)
endif(GUROBI_CPP_LIBRARY)

if (GUROBI_header_PATH AND GUROBI_C_LIBRARY AND GUROBI_CPP_LIBRARY)
   set(GUROBI_FOUND "YES")
   set(GUROBI_INCLUDEDIR ${GUROBI_header_PATH})
   set(GUROBI_LIBRARIES ${GUROBI_CPP_LIBRARY} ${GUROBI_C_LIBRARY})
endif(GUROBI_header_PATH AND GUROBI_C_LIBRARY AND GUROBI_CPP_LIBRARY)
