// some macros to organize similar grid index methods

#define ZERO_SHIFT FP3(0.0, 0.0, 0.0)

// method to return coordinates by grid index
#define GRID_GET_POSITION_IMPL(funcname, shift)                    \
forceinline const FP3 funcname(int x, int y, int z) const          \
{                                                                  \
    return baseCoords(x, y, z) + shift;                            \
}

// method to return grid index by coordinates
#define GRID_GET_INDEX_IMPL(funcname, shift)                       \
forceinline Int3 funcname(const FP3& coords) const {               \
    FP3 internalCoords;                                            \
    Int3 idx;                                                      \
    getGridCoords(coords, shift, idx, internalCoords);             \
    return idx;                                                    \
}

// method to return closest (left or right) grid index for given physical coords
#define GRID_GET_CLOSEST_INDEX_IMPL(funcname, shift)               \
forceinline Int3 funcname(const FP3& coords) const {               \
    FP3 internalCoords;                                            \
    Int3 idx;                                                      \
    getClosestGridCoords(coords, shift, idx, internalCoords);      \
    return idx;                                                    \
}

// method to return true if coords is inside of the area that grid defines (without external cells)
#define GRID_IS_INSIDE_IMPL(funcname, shift)                       \
forceinline bool funcname(const FP3& coords) const {               \
    return isInside(coords, shift);                                \
}

// method to return interpolated field value in arbitrary coords
#define GRID_GET_FIELD_IMPL(funcname, interpolation)               \
forceinline FP funcname(const FP3& coords) const                   \
{                                                                  \
    return (this->*interpolation)(coords);                         \
}

// method to return CIC-interpolated field value in arbitrary coords
#define GRID_GET_FIELD_CIC_IMPL(funcname, field, shift)            \
forceinline FP funcname(const FP3& coords) const                   \
{                                                                  \
    return getFieldCIC(coords, field, shift);                      \
}

// method to return TSC-interpolated field value in arbitrary coords
#define GRID_GET_FIELD_TSC_IMPL(funcname, field, shift)            \
forceinline FP funcname(const FP3& coords) const                   \
{                                                                  \
    return getFieldTSC(coords, field, shift);                      \
}

// method to return PCS-interpolated field value in arbitrary coords
#define GRID_GET_FIELD_PCS_IMPL(funcname, field, shift)            \
forceinline FP funcname(const FP3& coords) const                   \
{                                                                  \
    return getFieldPCS(coords, field, shift);                      \
}

// method to return second order interpolated field value in arbitrary coords
#define GRID_GET_FIELD_SECOND_ORDER_IMPL(funcname, field, shift)   \
forceinline FP funcname(const FP3& coords) const                   \
{                                                                  \
    return getFieldSecondOrder(coords, field, shift);              \
}

// method to return fourth order interpolated field value in arbitrary coords
#define GRID_GET_FIELD_FOURTH_ORDER_IMPL(funcname, field, shift)   \
forceinline FP funcname(const FP3& coords) const                   \
{                                                                  \
    return getFieldFourthOrder(coords, field, shift);              \
}
