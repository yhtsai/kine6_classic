project(kine6_classic)

add_executable(kine6_classic main.c Classic6dofKine.c)

target_link_libraries(kine6_classic m)

#install(TARGETS kine6_classic RUNTIME DESTINATION bin) #kdevelop
