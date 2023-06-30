#pragma once

#define SET_FIELD_CONFIGURATIONS_GRID_METHODS(pyFieldType)                \
    .def("set", &pyFieldType::setFieldConfiguration<NullField>,           \
        py::arg("field_configuration"))                                   \
    .def("set", &pyFieldType::setFieldConfiguration<TightFocusingField>,  \
        py::arg("field_configuration")) 


#define SET_COMPUTATIONAL_GRID_METHODS(pyFieldType)                        \
    .def("set_J", &pyFieldType::setJ)                                      \
    .def("set_E", &pyFieldType::setE)                                      \
    .def("set_B", &pyFieldType::setB)                                      \
    .def("set_J", &pyFieldType::pySetJ)                                    \
    .def("set_E", &pyFieldType::pySetE)                                    \
    .def("set_B", &pyFieldType::pySetB)                                    \
    .def("set_J", &pyFieldType::setJxyz,                                   \
        py::arg("Jx"), py::arg("Jy"), py::arg("Jz"))                       \
    .def("set_E", &pyFieldType::setExyz,                                   \
        py::arg("Ex"), py::arg("Ey"), py::arg("Ez"))                       \
    .def("set_B", &pyFieldType::setBxyz,                                   \
        py::arg("Bx"), py::arg("By"), py::arg("Bz"))                       \
    .def("set_J", &pyFieldType::pySetJxyz,                                 \
        py::arg("Jx"), py::arg("Jy"), py::arg("Jz"))                       \
    .def("set_E", &pyFieldType::pySetExyz,                                 \
        py::arg("Ex"), py::arg("Ey"), py::arg("Ez"))                       \
    .def("set_B", &pyFieldType::pySetBxyz,                                 \
        py::arg("Bx"), py::arg("By"), py::arg("Bz"))                       \
    .def("set_J", &pyFieldType::setJxyzt,                                  \
        py::arg("Jx"), py::arg("Jy"), py::arg("Jz"), py::arg("t"))         \
    .def("set_E", &pyFieldType::setExyzt,                                  \
        py::arg("Ex"), py::arg("Ey"), py::arg("Ez"), py::arg("t"))         \
    .def("set_B", &pyFieldType::setBxyzt,                                  \
        py::arg("Bx"), py::arg("By"), py::arg("Bz"), py::arg("t"))         \
    SET_FIELD_CONFIGURATIONS_GRID_METHODS(pyFieldType)


#define SET_SCALAR_FIELD_METHODS(pyFieldType)                              \
    .def("get_Jx_array", &pyFieldType::getJxArray)                         \
    .def("get_Jy_array", &pyFieldType::getJyArray)                         \
    .def("get_Jz_array", &pyFieldType::getJzArray)                         \
    .def("get_Ex_array", &pyFieldType::getExArray)                         \
    .def("get_Ey_array", &pyFieldType::getEyArray)                         \
    .def("get_Ez_array", &pyFieldType::getEzArray)                         \
    .def("get_Bx_array", &pyFieldType::getBxArray)                         \
    .def("get_By_array", &pyFieldType::getByArray)                         \
    .def("get_Bz_array", &pyFieldType::getBzArray)


#define SET_COMMON_FIELD_METHODS(pyFieldType)                             \
    .def("change_time_step", &pyFieldType::changeTimeStep,                \
        py::arg("time_step"))                                             \
    .def("refresh", &pyFieldType::refresh)                                \
    .def("set_time", &pyFieldType::setTime, py::arg("time"))              \
    .def("get_time", &pyFieldType::getTime)


#define SET_SUM_AND_MAP_FIELD_METHODS(pyFieldType)                        \
    .def("apply_mapping", [](std::shared_ptr<pyFieldType> self,           \
        std::shared_ptr<Mapping> mapping) {                               \
        return self->applyMapping(                                        \
            std::static_pointer_cast<pyFieldBase>(self), mapping          \
            );                                                            \
    }, py::arg("mapping"))                                                \
    .def("__add__", [](std::shared_ptr<pyFieldType> self,                 \
        std::shared_ptr<pyFieldBase> other) {                             \
        return std::make_shared<pySumField>(                              \
            std::static_pointer_cast<pyFieldBase>(self), other            \
            );                                                            \
    }, py::is_operator())                                                 \
    .def("__mul__", [](std::shared_ptr<pyFieldType> self, FP factor) {    \
        return std::make_shared<pyMulField>(                              \
            std::static_pointer_cast<pyFieldBase>(self), factor           \
            );                                                            \
    }, py::is_operator())                                                 \
    .def("__rmul__", [](std::shared_ptr<pyFieldType> self, FP factor) {   \
        return std::make_shared<pyMulField>(                              \
            std::static_pointer_cast<pyFieldBase>(self), factor           \
            );                                                            \
    }, py::is_operator())


// ------------------- Methods for getting field sections -------------------


#define SET_PYFIELDBASE_XY_PLANE_METHOD(name, method)                                  \
    .def(name, [](std::shared_ptr<pyFieldBase> self, FP zpos,                          \
        FP xmin, FP xmax, size_t nx, FP ymin, FP ymax, size_t ny) {                    \
            return self->getSlice2d(CoordinateEnum::z, zpos, CoordinateEnum::x,        \
                xmin, xmax, nx, CoordinateEnum::y, ymin, ymax, ny, method);            \
        }, py::arg("z_pos"), py::arg("x_min"), py::arg("x_max"), py::arg("x_size"),    \
        py::arg("y_min"), py::arg("y_max"), py::arg("y_size"))

#define SET_PYFIELDBASE_XZ_PLANE_METHOD(name, method)                                  \
    .def(name, [](std::shared_ptr<pyFieldBase> self, FP ypos,                          \
        FP xmin, FP xmax, size_t nx, FP zmin, FP zmax, size_t nz) {                    \
            return self->getSlice2d(CoordinateEnum::y, ypos, CoordinateEnum::x,        \
                xmin, xmax, nx, CoordinateEnum::z, zmin, zmax, nz, method);            \
        }, py::arg("y_pos"), py::arg("x_min"), py::arg("x_max"), py::arg("x_size"),    \
        py::arg("z_min"), py::arg("z_max"), py::arg("z_size"))

#define SET_PYFIELDBASE_YZ_PLANE_METHOD(name, method)                                  \
    .def(name, [](std::shared_ptr<pyFieldBase> self, FP xpos,                          \
        FP ymin, FP ymax, size_t ny, FP zmin, FP zmax, size_t nz) {                    \
            return self->getSlice2d(CoordinateEnum::x, xpos, CoordinateEnum::y,        \
                ymin, ymax, ny, CoordinateEnum::z, zmin, zmax, nz, method);            \
        }, py::arg("x_pos"), py::arg("y_min"), py::arg("y_max"), py::arg("y_size"),    \
        py::arg("z_min"), py::arg("z_max"), py::arg("z_size"))


#define SET_PYFIELDBASE_X_LINE_METHOD(name, method)                                    \
    .def(name, [](std::shared_ptr<pyFieldBase> self, FP ypos, FP zpos,                 \
        FP xmin, FP xmax, size_t nx) {                                                 \
            return self->getSlice1d(CoordinateEnum::y, ypos,                           \
                CoordinateEnum::z, zpos, CoordinateEnum::x, xmin, xmax, nx, method);   \
        }, py::arg("y_pos"), py::arg("z_pos"),                                         \
        py::arg("x_min"), py::arg("x_max"), py::arg("x_size"))

#define SET_PYFIELDBASE_Y_LINE_METHOD(name, method)                                    \
    .def(name, [](std::shared_ptr<pyFieldBase> self, FP xpos, FP zpos,                 \
        FP ymin, FP ymax, size_t ny) {                                                 \
            return self->getSlice1d(CoordinateEnum::x, xpos,                           \
                CoordinateEnum::z, zpos, CoordinateEnum::y, ymin, ymax, ny, method);   \
        }, py::arg("x_pos"), py::arg("z_pos"),                                         \
        py::arg("y_min"), py::arg("y_max"), py::arg("y_size"))

#define SET_PYFIELDBASE_Z_LINE_METHOD(name, method)                                    \
    .def(name, [](std::shared_ptr<pyFieldBase> self, FP xpos, FP ypos,                 \
        FP zmin, FP zmax, size_t nz) {                                                 \
            return self->getSlice1d(CoordinateEnum::x, xpos,                           \
                CoordinateEnum::y, ypos, CoordinateEnum::z, zmin, zmax, nz, method);   \
        }, py::arg("x_pos"), py::arg("y_pos"),                                         \
        py::arg("z_min"), py::arg("z_max"), py::arg("z_size"))


#define SET_PYFIELDBASE_3D_AREA_METHOD(name, method)                                   \
    .def(name, [](std::shared_ptr<pyFieldBase> self, FP xmin, FP xmax, size_t nx,      \
        FP ymin, FP ymax, size_t ny, FP zmin, FP zmax, size_t nz) {                    \
            return self->getSlice3d(CoordinateEnum::x, xmin, xmax, nx,                 \
                CoordinateEnum::y, ymin, ymax, ny, CoordinateEnum::z, zmin, zmax, nz,  \
                method);                                                               \
        }, py::arg("x_min"), py::arg("x_max"), py::arg("x_size"),                      \
        py::arg("y_min"), py::arg("y_max"), py::arg("y_size"),                         \
        py::arg("z_min"), py::arg("z_max"), py::arg("z_size"))


#define SET_PYFIELDBASE_FIELD_SECTION_METHODS(macros, name)                            \
    macros("get_Ex" name, &pyFieldBase::getEx)                                         \
    macros("get_Ey" name, &pyFieldBase::getEy)                                         \
    macros("get_Ez" name, &pyFieldBase::getEz)                                         \
    macros("get_Bx" name, &pyFieldBase::getBx)                                         \
    macros("get_By" name, &pyFieldBase::getBy)                                         \
    macros("get_Bz" name, &pyFieldBase::getBz)                                         \
    macros("get_Jx" name, &pyFieldBase::getJx)                                         \
    macros("get_Jy" name, &pyFieldBase::getJy)                                         \
    macros("get_Jz" name, &pyFieldBase::getJz)

#define SET_ALL_PYFIELDBASE_SECTION_METHODS()                                          \
    SET_PYFIELDBASE_FIELD_SECTION_METHODS(SET_PYFIELDBASE_XY_PLANE_METHOD, "_xy_plane")\
    SET_PYFIELDBASE_FIELD_SECTION_METHODS(SET_PYFIELDBASE_XZ_PLANE_METHOD, "_xz_plane")\
    SET_PYFIELDBASE_FIELD_SECTION_METHODS(SET_PYFIELDBASE_YZ_PLANE_METHOD, "_yz_plane")\
    SET_PYFIELDBASE_FIELD_SECTION_METHODS(SET_PYFIELDBASE_X_LINE_METHOD, "_x_line")    \
    SET_PYFIELDBASE_FIELD_SECTION_METHODS(SET_PYFIELDBASE_Y_LINE_METHOD, "_y_line")    \
    SET_PYFIELDBASE_FIELD_SECTION_METHODS(SET_PYFIELDBASE_Z_LINE_METHOD, "_z_line")    \
    SET_PYFIELDBASE_FIELD_SECTION_METHODS(SET_PYFIELDBASE_3D_AREA_METHOD, "")


// ------------------- Field generator setter methods -------------------


#define SET_FIELD_GENERATOR_METHODS()                                                                 \
    .def("set_field_generator", &pyYeeField::setFieldGenerator,                                       \
        py::arg("left_index"), py::arg("right_index"),                                                \
        py::arg("bx_func"), py::arg("by_func"), py::arg("bz_func"),                                   \
        py::arg("ex_func"), py::arg("ey_func"), py::arg("ez_func"),                                   \
        py::arg("is_left_x_border_enabled") = true, py::arg("is_left_y_border_enabled") = true,       \
        py::arg("is_left_z_border_enabled") = true, py::arg("is_right_x_border_enabled") = true,      \
        py::arg("is_right_y_border_enabled") = true, py::arg("is_right_z_border_enabled") = true)     \
    .def("set_field_generator", &pyYeeField::setFieldGeneratorAllFunctions,                           \
        py::arg("left_index"), py::arg("right_index"),                                                \
        py::arg("left_x_bx_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_x_by_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_x_bz_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("right_x_bx_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_x_by_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_x_bz_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("left_y_bx_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_y_by_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_y_bz_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("right_y_bx_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_y_by_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_y_bz_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("left_z_bx_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_z_by_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_z_bz_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("right_z_bx_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_z_by_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_z_bz_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("left_x_ex_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_x_ey_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_x_ez_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("right_x_ex_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_x_ey_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_x_ez_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("left_y_ex_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_y_ey_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_y_ez_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("right_y_ex_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_y_ey_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_y_ez_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("left_z_ex_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_z_ey_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("left_z_ez_func") = (CFunctionPointer)field_generator::defaultFieldFunction,    \
        py::arg("right_z_ex_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_z_ey_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("right_z_ez_func") = (CFunctionPointer)field_generator::defaultFieldFunction,   \
        py::arg("is_left_x_border_enabled") = true, py::arg("is_left_y_border_enabled") = true,       \
        py::arg("is_left_z_border_enabled") = true, py::arg("is_right_x_border_enabled") = true,      \
        py::arg("is_right_y_border_enabled") = true, py::arg("is_right_z_border_enabled") = true)
