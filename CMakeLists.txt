cmake_minimum_required(VERSION 3.15)
project(proof2)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS}")
#set(CMAKE_BUILD_TYPE  Release)

add_executable(proof2 main.c bn256.c bn256.h util.h aztec.c aztec.h util.c)
include_directories(/usr/local/include )
target_link_libraries(proof2 "-L /usr/local/lib -L /usr/local/opt/openssl/lib -lgmp -lgmpxx -lcrypto")