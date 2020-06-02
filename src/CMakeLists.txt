
feelpp_add_application(blockform SRCS mqs-blockform.cpp INSTALL TESTS)
feelpp_add_application(form SRCS mqs-form.cpp INSTALL TESTS)

feelpp_add_testcase( cases PREFIX mqs CATEGORY mqs )
