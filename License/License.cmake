#--////////////////////////////////////////////////////////////////////////////
#-- Copyright (c) 2009, Michael A. Jackson. BlueQuartz Software
#-- All rights reserved.
#-- BSD License: http://www.opensource.org/licenses/bsd-license.html
#-- This code was partly written under US Air Force Contract FA8650-07-D-5800
#--////////////////////////////////////////////////////////////////////////////
# --------------------------------------------------------------------
# create Resource files for the various license files that are used and 
# also create a header file that lists all the License Files
set(LICENSE_FILES ${PROJECT_BINARY_DIR}/License.txt
                  ${PROJECT_SOURCE_DIR}/License/Boost.license
                  ${PROJECT_SOURCE_DIR}/License/MXA.license
                  ${PROJECT_SOURCE_DIR}/License/Qt.license
                  ${PROJECT_SOURCE_DIR}/License/tiff.license )
set(QRC_LICENSE_FILES "")
file(WRITE ${PROJECT_BINARY_DIR}/License/LicenseFiles.h "#ifndef _LICENSE_FILES_H_\n")
file(APPEND ${PROJECT_BINARY_DIR}/License/LicenseFiles.h "#define _LICENSE_FILES_H_\n")
file(APPEND ${PROJECT_BINARY_DIR}/License/LicenseFiles.h "namespace ${PROJECT_PREFIX} {\n")
file(APPEND ${PROJECT_BINARY_DIR}/License/LicenseFiles.h "  QStringList LicenseList = (QStringList()  ")
foreach(lf ${LICENSE_FILES})
    get_filename_component(cmp_text_file_name ${lf} NAME)
    get_filename_component(lf_fn ${lf} NAME_WE)
    # Copy the text file into the Build Directory
    configure_file("${lf}" ${PROJECT_BINARY_DIR}/License/${cmp_text_file_name} COPYONLY )
        
    # create the Qt Resource File
    configure_file(${CMP_CONFIGURED_FILES_SOURCE_DIR}/QtResourceFile.qrc.in 
                   ${PROJECT_BINARY_DIR}/License/${lf_fn}.qrc)
                   
    set(QRC_LICENSE_FILES ${QRC_LICENSE_FILES} ${PROJECT_BINARY_DIR}/License/${lf_fn}.qrc)
    file(APPEND ${PROJECT_BINARY_DIR}/License/LicenseFiles.h " << \":/${cmp_text_file_name}\"")
    
endforeach(lf ${LICENSE_FILES})
file(APPEND ${PROJECT_BINARY_DIR}/License/LicenseFiles.h ");\n")
file(APPEND ${PROJECT_BINARY_DIR}/License/LicenseFiles.h "}\n#endif /* _LICENSE_FILES_H_ */ \n")

