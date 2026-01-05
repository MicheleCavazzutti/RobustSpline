message(STATUS "Configuring QDLDL solver from local sources")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

# Invece di FetchContent_Declare con Git, puntiamo alla cartella locale
# Assicurati che il path sia corretto rispetto a dove si trova questo file
set(qdldl_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/qdldl_local")

# Definiamo i tipi come farebbe OSQP
set(QDLDL_FLOAT ${OSQP_USE_FLOAT} CACHE BOOL "QDLDL Float type")
set(QDLDL_LONG ${OSQP_USE_LONG} CACHE BOOL "QDLDL Integer type")
set(QDLDL_BUILD_STATIC_LIB OFF CACHE BOOL "Build QDLDL static library")
set(QDLDL_BUILD_SHARED_LIB OFF CACHE BOOL "Build QDLDL shared library")

# Aggiungiamo i sorgenti di QDLDL direttamente
add_subdirectory(${qdldl_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/qdldl_build EXCLUDE_FROM_ALL)

list(POP_BACK CMAKE_MESSAGE_INDENT)

# Il resto del file rimane simile, ma puntiamo ai path locali
set(qdldl_include "${qdldl_SOURCE_DIR}/include")

file(
    GLOB
    AMD_SRC_FILES
    CONFIGURE_DEPENDS
    ${OSQP_ALGEBRA_ROOT}/_common/lin_sys/qdldl/amd/src/*.c
    ${OSQP_ALGEBRA_ROOT}/_common/lin_sys/qdldl/amd/include/*.h
    )

set( LIN_SYS_QDLDL_NON_EMBEDDED_SRC_FILES
     ${AMD_SRC_FILES}
     )

set( LIN_SYS_QDLDL_EMBEDDED_SRC_FILES
     ${OSQP_ALGEBRA_ROOT}/_common/kkt.h
     ${OSQP_ALGEBRA_ROOT}/_common/kkt.c
     ${OSQP_ALGEBRA_ROOT}/_common/lin_sys/qdldl/qdldl_interface.h
     ${OSQP_ALGEBRA_ROOT}/_common/lin_sys/qdldl/qdldl_interface.c
     )

set( LIN_SYS_QDLDL_SRC_FILES
     ${LIN_SYS_QDLDL_EMBEDDED_SRC_FILES}
     ${LIN_SYS_QDLDL_NON_EMBEDDED_SRC_FILES}
     )

set( LIN_SYS_QDLDL_INC_PATHS
     ${qdldl_include}
     ${OSQP_ALGEBRA_ROOT}/_common/
     ${OSQP_ALGEBRA_ROOT}/_common/lin_sys/qdldl/
     ${OSQP_ALGEBRA_ROOT}/_common/lin_sys/qdldl/amd/include
     ${qdldl_SOURCE_DIR}/include
     ${CMAKE_CURRENT_BINARY_DIR}/qdldl_build/include
     )