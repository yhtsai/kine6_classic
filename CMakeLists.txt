project(kine_classic_6dof)

add_executable(kine_classic_6dof main.c Classic6dofKine.c)

target_link_libraries(kine_classic_6dof m)

install(TARGETS kine_classic_6dof RUNTIME DESTINATION bin)
