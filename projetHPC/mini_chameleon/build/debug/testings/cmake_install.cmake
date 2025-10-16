# Install script for directory: /home/cisd-lahlah/mini-chameleon/testings

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/gnu/store/wvl8j6s9633m4wf5w3y0d2ld27wiz8xm-profile/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_ddot" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_ddot")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_ddot"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/cisd-lahlah/mini-chameleon/build/debug/testings/check_ddot")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_ddot" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_ddot")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_ddot"
         OLD_RPATH "/gnu/store/krxcx4inxn2p23zbygyfx5cn6vq5s8i7-starpu-1.4.8/lib:/gnu/store/vjdc7vss30gkrjnbwhhba6rzddwi5qmv-hwloc-2.12.2-lib/lib:/home/cisd-lahlah/mini-chameleon/build/debug/myblas:/home/cisd-lahlah/mini-chameleon/build/debug/algonum:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/gnu/store/wvl8j6s9633m4wf5w3y0d2ld27wiz8xm-profile/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_ddot")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgemm" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgemm")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgemm"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/cisd-lahlah/mini-chameleon/build/debug/testings/check_dgemm")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgemm" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgemm")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgemm"
         OLD_RPATH "/gnu/store/krxcx4inxn2p23zbygyfx5cn6vq5s8i7-starpu-1.4.8/lib:/gnu/store/vjdc7vss30gkrjnbwhhba6rzddwi5qmv-hwloc-2.12.2-lib/lib:/home/cisd-lahlah/mini-chameleon/build/debug/myblas:/home/cisd-lahlah/mini-chameleon/build/debug/algonum:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/gnu/store/wvl8j6s9633m4wf5w3y0d2ld27wiz8xm-profile/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgemm")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgetrf" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgetrf")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgetrf"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/cisd-lahlah/mini-chameleon/build/debug/testings/check_dgetrf")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgetrf" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgetrf")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgetrf"
         OLD_RPATH "/gnu/store/krxcx4inxn2p23zbygyfx5cn6vq5s8i7-starpu-1.4.8/lib:/gnu/store/vjdc7vss30gkrjnbwhhba6rzddwi5qmv-hwloc-2.12.2-lib/lib:/home/cisd-lahlah/mini-chameleon/build/debug/myblas:/home/cisd-lahlah/mini-chameleon/build/debug/algonum:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/gnu/store/wvl8j6s9633m4wf5w3y0d2ld27wiz8xm-profile/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/check_dgetrf")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_ddot" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_ddot")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_ddot"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/cisd-lahlah/mini-chameleon/build/debug/testings/perf_ddot")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_ddot" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_ddot")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_ddot"
         OLD_RPATH "/gnu/store/krxcx4inxn2p23zbygyfx5cn6vq5s8i7-starpu-1.4.8/lib:/gnu/store/vjdc7vss30gkrjnbwhhba6rzddwi5qmv-hwloc-2.12.2-lib/lib:/home/cisd-lahlah/mini-chameleon/build/debug/myblas:/home/cisd-lahlah/mini-chameleon/build/debug/algonum:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/gnu/store/wvl8j6s9633m4wf5w3y0d2ld27wiz8xm-profile/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_ddot")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgemm" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgemm")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgemm"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/cisd-lahlah/mini-chameleon/build/debug/testings/perf_dgemm")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgemm" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgemm")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgemm"
         OLD_RPATH "/gnu/store/krxcx4inxn2p23zbygyfx5cn6vq5s8i7-starpu-1.4.8/lib:/gnu/store/vjdc7vss30gkrjnbwhhba6rzddwi5qmv-hwloc-2.12.2-lib/lib:/home/cisd-lahlah/mini-chameleon/build/debug/myblas:/home/cisd-lahlah/mini-chameleon/build/debug/algonum:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/gnu/store/wvl8j6s9633m4wf5w3y0d2ld27wiz8xm-profile/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgemm")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgetrf" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgetrf")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgetrf"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/cisd-lahlah/mini-chameleon/build/debug/testings/perf_dgetrf")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgetrf" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgetrf")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgetrf"
         OLD_RPATH "/gnu/store/krxcx4inxn2p23zbygyfx5cn6vq5s8i7-starpu-1.4.8/lib:/gnu/store/vjdc7vss30gkrjnbwhhba6rzddwi5qmv-hwloc-2.12.2-lib/lib:/home/cisd-lahlah/mini-chameleon/build/debug/myblas:/home/cisd-lahlah/mini-chameleon/build/debug/algonum:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/gnu/store/wvl8j6s9633m4wf5w3y0d2ld27wiz8xm-profile/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/perf_dgetrf")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_ddot" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_ddot")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_ddot"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/cisd-lahlah/mini-chameleon/build/debug/testings/warm_ddot")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_ddot" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_ddot")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_ddot"
         OLD_RPATH "/gnu/store/krxcx4inxn2p23zbygyfx5cn6vq5s8i7-starpu-1.4.8/lib:/gnu/store/vjdc7vss30gkrjnbwhhba6rzddwi5qmv-hwloc-2.12.2-lib/lib:/home/cisd-lahlah/mini-chameleon/build/debug/myblas:/home/cisd-lahlah/mini-chameleon/build/debug/algonum:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/gnu/store/wvl8j6s9633m4wf5w3y0d2ld27wiz8xm-profile/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_ddot")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgemm" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgemm")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgemm"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/cisd-lahlah/mini-chameleon/build/debug/testings/warm_dgemm")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgemm" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgemm")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgemm"
         OLD_RPATH "/gnu/store/krxcx4inxn2p23zbygyfx5cn6vq5s8i7-starpu-1.4.8/lib:/gnu/store/vjdc7vss30gkrjnbwhhba6rzddwi5qmv-hwloc-2.12.2-lib/lib:/home/cisd-lahlah/mini-chameleon/build/debug/myblas:/home/cisd-lahlah/mini-chameleon/build/debug/algonum:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/gnu/store/wvl8j6s9633m4wf5w3y0d2ld27wiz8xm-profile/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgemm")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgetrf" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgetrf")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgetrf"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/cisd-lahlah/mini-chameleon/build/debug/testings/warm_dgetrf")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgetrf" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgetrf")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgetrf"
         OLD_RPATH "/gnu/store/krxcx4inxn2p23zbygyfx5cn6vq5s8i7-starpu-1.4.8/lib:/gnu/store/vjdc7vss30gkrjnbwhhba6rzddwi5qmv-hwloc-2.12.2-lib/lib:/home/cisd-lahlah/mini-chameleon/build/debug/myblas:/home/cisd-lahlah/mini-chameleon/build/debug/algonum:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/gnu/store/wvl8j6s9633m4wf5w3y0d2ld27wiz8xm-profile/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/warm_dgetrf")
    endif()
  endif()
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/cisd-lahlah/mini-chameleon/build/debug/testings/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
