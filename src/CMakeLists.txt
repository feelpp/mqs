
if (EXISTS ${FEELPP_INSTALL_DIR}/include/feelpp/toolboxes)
  MESSAGE(STATUS "Loading toolboxes hearders from ${FEELPP_INSTALL_DIR}/include/feelpp/toolboxes")
  include_directories(${FEELPP_INSTALL_DIR}/include/feelpp/toolboxes)
endif()

feelpp_add_application(blockform SRCS mqs-blockform.cpp INSTALL)
feelpp_add_application(form SRCS mqs-form.cpp
  LINK_LIBRARIES ${FEELPP_LIBRARIES} feelpp_modelcore feelpp_modelalg feelpp_toolbox_thermoelectric_lib_3dP1G1_3dP1G1 libfeelpp_toolbox_maxwell_lib_3dP1G1.so
  INSTALL)

feelpp_add_testcase( cases PREFIX mqs CATEGORY mqs )
