cmake_minimum_required(VERSION 3.11)
project(Practical3)

include(FetchContent)

# libigl
FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG        e699d257fe860d4fd82709835e72438bd0cd1533
)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_TETGEN            "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI  "Use OpenGL viewer with ImGui"  ON)
FetchContent_MakeAvailable(libigl)

# Add your project files
file(GLOB SRCFILES *.cpp *.h *.c)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::tetgen igl::opengl_glfw_imgui)
