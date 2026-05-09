
target_sources(OPS_Runtime PRIVATE
    element.cpp
    section/fiber.cpp
)

add_subdirectory(material)
add_subdirectory(element)
add_subdirectory(invoking)
