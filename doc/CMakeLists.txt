
find_package(LATEX)
if (LATEX_FOUND)
include(UseLATEX)
# set(LATEX_OUTPUT_PATH doc)

add_subdirectory(SimpleManual/manual2.5)

add_subdirectory(SimpleManual/DevelGuide)

add_subdirectory(SimpleTutorials2017)

endif()
