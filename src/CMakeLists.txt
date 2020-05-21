SET(SRC
  test2
  )
foreach(app ${SRC})
  if ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${app}.cfg AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${app}.geo  )
    feelpp_add_application( ${app} SRCS ${app}.cpp  CFG ${app}.cfg GEO ${app}.geo TESTS)
  elseif ( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${app}.cfg )
    feelpp_add_application( ${app} SRCS ${app}.cpp  CFG ${app}.cfg TESTS )
  else()
    feelpp_add_application( ${app} SRCS ${app}.cpp TESTS )
  endif()
endforeach()

feelpp_add_testcase(cases)
