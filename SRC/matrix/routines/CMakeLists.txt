#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Matrix
  PRIVATE
    eigSY3.c
    invGL2.c
    invGL3.c
    invGL4.c
    invGL5.c
    invGL6.c
  PUBLIC
    cmx.h
    SY3.h
)

target_include_directories(OPS_Matrix PUBLIC ${CMAKE_CURRENT_LIST_DIR})
