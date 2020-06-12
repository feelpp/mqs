
feelpp_add_application(blockform SRCS mqs-blockform.cpp INSTALL)
feelpp_add_application(form SRCS mqs-form.cpp INSTALL)

feelpp_add_testcase( cases PREFIX mqs CATEGORY mqs )
