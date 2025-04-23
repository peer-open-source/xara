target_sources(OPS_Element
  PRIVATE
    IGASurfacePatch.cpp
    IGAFollowerLoad.cpp
    IGAKLShell.cpp
    IGAKLShell_BendingStrip.cpp
  PUBLIC
    IGASurfacePatch.h
    IGAFollowerLoad.h
    IGAKLShell.h
    IGAKLShell_BendingStrip.h
)
target_sources(OPS_Runtime PRIVATE dispatch.cpp)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
