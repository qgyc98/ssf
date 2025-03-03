# Install script for directory: D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/install/windows-default")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/GEFEL/cmake_install.cmake")
  include("D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/HPARAL/cmake_install.cmake")
  include("D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/MAN/cmake_install.cmake")
  include("D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/MEFEL/cmake_install.cmake")
  include("D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/METR/cmake_install.cmake")
  include("D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/PARGEF/cmake_install.cmake")
  include("D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/PARMEF/cmake_install.cmake")
  include("D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/PARMETR/cmake_install.cmake")
  include("D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/PARTRF/cmake_install.cmake")
  include("D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/PREP/cmake_install.cmake")
  include("D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/TRFEL/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "D:/work-bbefem/gongzuowenjian/SIFEL/svn-SIFEL/websvn2/ceshi/hh/SIFEL/out/build/windows-default/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
