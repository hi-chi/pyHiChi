#pragma once

#include "Dimension.h"
#include "FP.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>


namespace pfc {

// Vector of 1 component of arithmetic type T (T = int, double, etc.)
// with access by .x, and [], provides basic arithmetic operations
// This is essentially a wrapper around T with interface compatible with Vector2 and Vector3
template <typename T>
struct Vector1
{
    T x;

    Vector1() :
        x(0) {}

    Vector1(T _x) :
        x(_x) {}

    template<typename U>
    Vector1(const Vector1<U>& other) :
        x(other.x) {}

    inline T operator[](int ) const
    {
        return x;
    }

    inline T& operator[](int )
    {
        return x;
    }

    inline T volume() const
    {
        return x;
    }

    inline T norm() const
    {
        return x * ((x > 0) ? static_cast<T>(1) : static_cast<T>(-1));
    }

    inline T norm2() const
    {
        return x * x;
    }

    inline std::string toString()
    {
        std::stringstream s;
        s << *this;
        return s.str();
    }
};

template<typename T>
inline const Vector1<T> operator + (const Vector1<T>& v1, const Vector1<T>& v2)
{
    return Vector1<T>(v1.x + v2.x);
}

template<typename T>
inline Vector1<T>& operator += (Vector1<T>& v1, const Vector1<T>& v2)
{
    v1.x += v2.x;
    return v1;
}

template<typename T>
inline const Vector1<T> operator - (const Vector1<T>& v1, const Vector1<T>& v2)
{
    return Vector1<T>(v1.x - v2.x);
}

template<typename T>
inline Vector1<T>& operator -= (Vector1<T>& v1, const Vector1<T>& v2)
{
    v1.x -= v2.x;
    return v1;
}

template<typename T>
inline const Vector1<T> operator * (const Vector1<T>& v1, const Vector1<T>& v2)
{
    return Vector1<T>(v1.x * v2.x);
}

template<typename T>
inline const Vector1<T> operator * (const Vector1<T>& v, T a)
{
    return Vector1<T>(v.x * a);
}

template<typename T>
inline const Vector1<T> operator * (T a, const Vector1<T>& v)
{
    return Vector1<T>(v.x * a);
}

template<typename T>
inline Vector1<T>& operator *= (Vector1<T>& v1, const Vector1<T>& v2)
{
    v1.x *= v2.x;
    return v1;
}

template<typename T>
inline Vector1<T>& operator *= (Vector1<T>& v, T a)
{
    v.x *= a;
    return v;
}

template<typename T>
inline const Vector1<T> operator / (const Vector1<T>& v1, const Vector1<T>& v2)
{
    return Vector1<T>(v1.x / v2.x);
}

template<typename T>
inline Vector1<T>& operator /= (Vector1<T>& v1, const Vector1<T>& v2)
{
    v1.x /= v2.x;
    return v1;
}

template<typename T>
inline const Vector1<T> operator / (const Vector1<T>& v, T a)
{
    return Vector1<T>(v.x / a);
}

template<typename T>
inline Vector1<T>& operator /= (Vector1<T>& v, T a)
{
    v.x /= a;
    return v;
}

template<typename T>
inline bool operator == (const Vector1<T>& v1, const Vector1<T>& v2)
{
    return (v1.x == v2.x);
}

template<typename T>
inline bool operator != (const Vector1<T>& v1, const Vector1<T>& v2)
{
    return (v1.x != v2.x);
}

template<typename T>
inline bool operator < (const Vector1<T>& v1, const Vector1<T>& v2)
{
    return (v1.x < v2.x);
}

template<typename T>
inline bool operator <= (const Vector1<T>& v1, const Vector1<T>& v2)
{
    return (v1.x <= v2.x);
}

template<typename T>
inline bool operator > (const Vector1<T>& v1, const Vector1<T>& v2)
{
    return (v1.x > v2.x);
}

template<typename T>
inline bool operator >= (const Vector1<T>& v1, const Vector1<T>& v2)
{
    return (v1.x >= v2.x);
}

template<typename T>
inline std::ostream& operator<<(std::ostream& out, const Vector1<T>& v)
{
    return out << "(" << v.x << ")";
}

template<typename T>
inline T dot(const Vector1<T>& v1, const Vector1<T>& v2)
{
    return v1.x * v2.x;
}


// Vector of 2 components of arithmetic type T (T = int, double, etc.)
// with access by .x, .y, and [], provides basic arithmetic operations
template <typename T>
struct Vector2
{
    T x, y;

    Vector2() :
        x(0), y(0) {}

    Vector2(T _x, T _y) :
        x(_x), y(_y) {}

    template<typename U>
    Vector2(const Vector2<U>& other) :
        x(other.x), y(other.y) {}

    inline T operator[](int idx) const
    {
        return *((T*)this + idx);
    }

    inline T& operator[](int idx)
    {
        return *((T*)this + idx);
    }

    inline T volume() const
    {
        return x * y;
    }

    inline T norm() const
    {
        return sqrt(x * x + y * y);
    }

    inline T norm2() const
    {
        return x * x + y * y;
    }

    inline std::string toString()
    {
        std::stringstream s;
        s << *this;
        return s.str();
    }
};

template<typename T>
inline const Vector2<T> operator + (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return Vector2<T>(v1.x + v2.x, v1.y + v2.y);
}

template<typename T>
inline Vector2<T>& operator += (Vector2<T>& v1, const Vector2<T>& v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    return v1;
}

template<typename T>
inline const Vector2<T> operator - (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return Vector2<T>(v1.x - v2.x, v1.y - v2.y);
}

template<typename T>
inline Vector2<T>& operator -= (Vector2<T>& v1, const Vector2<T>& v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    return v1;
}

template<typename T>
inline const Vector2<T> operator * (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return Vector2<T>(v1.x * v2.x, v1.y * v2.y);
}

template<typename T>
inline Vector2<T>& operator *= (Vector2<T>& v1, const Vector2<T>& v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    return v1;
}

template<typename T>
inline const Vector2<T> operator * (const Vector2<T>& v, T a)
{
    return Vector2<T>(v.x * a, v.y * a);
}

template<typename T>
inline const Vector2<T> operator * (T a, const Vector2<T>& v)
{
    return Vector2<T>(v.x * a, v.y * a);
}

template<typename T>
inline Vector2<T>& operator *= (Vector2<T>& v, T a)
{
    v.x *= a;
    v.y *= a;
    return v;
}

template<typename T>
inline const Vector2<T> operator / (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return Vector2<T>(v1.x / v2.x, v1.y / v2.y);
}

template<typename T>
inline Vector2<T>& operator /= (Vector2<T>& v1, const Vector2<T>& v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    return v1;
}

template<typename T>
inline const Vector2<T> operator / (const Vector2<T>& v, T a)
{
    return Vector2<T>(v.x / a, v.y / a);
}

template<typename T>
inline Vector2<T>& operator /= (Vector2<T>& v, T a)
{
    v.x /= a;
    v.y /= a;
    return v;
}

template<typename T>
inline bool operator == (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y);
}

template<typename T>
inline bool operator != (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x != v2.x) || (v1.y != v2.y);
}

template<typename T>
inline bool operator < (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x < v2.x) && (v1.y < v2.y);
}

template<typename T>
inline bool operator <= (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x <= v2.x) && (v1.y <= v2.y);
}

template<typename T>
inline bool operator > (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x > v2.x) && (v1.y > v2.y);
}

template<typename T>
inline bool operator >= (const Vector2<T>& v1, const Vector2<T>& v2)
{
    return (v1.x >= v2.x) && (v1.y >= v2.y);
}

template<typename T>
inline std::ostream& operator<<(std::ostream& out, const Vector2<T>& v)
{
    return out << "(" << v.x << ", " << v.y << ")";
}

template<typename T>
inline T dot(const Vector2<T>& v1, const Vector2<T>& v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}


// Vector of 3 components of arithmetic type T (T = int, double, etc.)
// with access by .x, .y, .z and [], provides basic arithmetic operations
template <typename T>
struct Vector3
{
    T x, y, z;

    Vector3() :
        x(0), y(0), z(0) {}

    Vector3(T _x, T _y, T _z) :
        x(_x), y(_y), z(_z) {}

    template<typename U>
    Vector3(const Vector3<U>& other) :
        x(other.x), y(other.y), z(other.z) {}

    Vector3(Vector2<T> v) :
        x(v.x), y(v.y), z(0) {}

    __forceinline T operator[](int idx) const
    {
        return *((T*)this + idx);
    }

    __forceinline T& operator[](int idx)
    {
        return *((T*)this + idx);
    }

    __forceinline T volume() const
    {
        return x * y * z;
    }

    __forceinline T norm() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    __forceinline T norm2() const
    {
        return x * x + y * y + z * z;
    }

    __forceinline void normalize()
    {
        T no = norm();
        if (no != 0) {
            x /= no;
            y /= no;
            z /= no;
        }
    }

    inline std::string toString()
    {
        std::stringstream s;
        s << *this;
        return s.str();
    }
};

template<typename T>
__forceinline const Vector3<T> operator + (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return Vector3<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

template<typename T>
__forceinline Vector3<T>& operator += (Vector3<T>& v1, const Vector3<T>& v2)
{
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;
    return v1;
}

template<typename T>
__forceinline const Vector3<T> operator - (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return Vector3<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

template<typename T>
__forceinline const Vector3<T> operator - (const Vector3<T>& v1)
{
    return Vector3<T>(-v1.x, -v1.y, -v1.z);
}

template<typename T>
__forceinline Vector3<T>& operator -= (Vector3<T>& v1, const Vector3<T>& v2)
{
    v1.x -= v2.x;
    v1.y -= v2.y;
    v1.z -= v2.z;
    return v1;
}

template<typename T>
__forceinline const Vector3<T> operator * (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return Vector3<T>(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

template<typename T>
__forceinline Vector3<T>& operator *= (Vector3<T>& v1, const Vector3<T>& v2)
{
    v1.x *= v2.x;
    v1.y *= v2.y;
    v1.z *= v2.z;
    return v1;
}

template<typename T>
__forceinline const Vector3<T> operator * (const Vector3<T>& v, T a)
{
    return Vector3<T>(v.x * a, v.y * a, v.z * a);
}

template<typename T>
__forceinline const Vector3<T> operator * (T a, const Vector3<T>& v2)
{
    return Vector3<T>(a * v2.x, a * v2.y, a * v2.z);
}

template<typename T>
__forceinline Vector3<T>& operator *= (Vector3<T>& v, T a)
{
    v.x *= a;
    v.y *= a;
    v.z *= a;
    return v;
}

template<typename T>
__forceinline const Vector3<T> operator / (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return Vector3<T>(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
}

template<typename T>
__forceinline Vector3<T>& operator /= (Vector3<T>& v1, const Vector3<T>& v2)
{
    v1.x /= v2.x;
    v1.y /= v2.y;
    v1.z /= v2.z;
    return v1;
}

template<typename T>
__forceinline const Vector3<T> operator / (const Vector3<T>& v, T a)
{
    return Vector3<T>(v.x / a, v.y / a, v.z / a);
}

template<typename T>
__forceinline Vector3<T>& operator /= (Vector3<T>& v, T a)
{
    v.x /= a;
    v.y /= a;
    v.z /= a;
    return v;
}

template<typename T>
__forceinline bool operator == (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z);
}

template<typename T>
__forceinline bool operator != (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x != v2.x) || (v1.y != v2.y) || (v1.z != v2.z);
}

template<typename T>
__forceinline bool operator < (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x < v2.x) && (v1.y < v2.y) && (v1.z < v2.z);
}

template<typename T>
__forceinline bool operator <= (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x <= v2.x) && (v1.y <= v2.y) && (v1.z <= v2.z);
}

template<typename T>
__forceinline bool operator > (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x > v2.x) && (v1.y > v2.y) && (v1.z > v2.z);
}

template<typename T>
__forceinline bool operator >= (const Vector3<T>& v1, const Vector3<T>& v2)
{
    return (v1.x >= v2.x) && (v1.y >= v2.y) && (v1.z >= v2.z);
}

template<typename T>
inline std::ostream& operator<<(std::ostream& out, const Vector3<T>& v)
{
    return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

template<typename T>
__forceinline T dot(const Vector3<T>& v1, const Vector3<T>& v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template<typename T>
__forceinline const Vector3<T> cross(const Vector3<T>& v1, const Vector3<T>& v2)
{
    return Vector3<T>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x);
}


template<Dimension dimension, typename T>
struct VectorTypeHelper {
};

template<typename T>
struct VectorTypeHelper<One, T> {
    typedef Vector1<T> Type;
};

template<typename T>
struct VectorTypeHelper<Two, T> {
    typedef Vector2<T> Type;
};

template<typename T>
struct VectorTypeHelper<Three, T> {
    typedef Vector3<T> Type;
};


template<int dimension, typename T>
struct VectorIntTypeHelper {
};

template<typename T>
struct VectorIntTypeHelper<1, T> {
    typedef Vector1<T> Type;
};

template<typename T>
struct VectorIntTypeHelper<2, T> {
    typedef Vector2<T> Type;
};

template<typename T>
struct VectorIntTypeHelper<3, T> {
    typedef Vector3<T> Type;
};


template<typename Vector>
struct VectorDimensionHelper {
};

template<typename T>
struct VectorDimensionHelper<Vector1<T> > {
    static const int dimension = 1;
};

template<typename T>
struct VectorDimensionHelper<Vector2<T> > {
    static const int dimension = 2;

};

template<typename T>
struct VectorDimensionHelper<Vector3<T> > {
    static const int dimension = 3;
};

typedef Vector3<int> Int3;
typedef Vector3<FP> FP3;
typedef Vector3<complexFP> ComplexFP3;

inline Int3 operator%(const Int3& v1, const Int3& v2) {
    return Int3((v1.x + v2.x) % v2.x, (v1.y + v2.y) % v2.y, (v1.z + v2.z) % v2.z);
}

template<typename Real>
inline const Vector1<int> floor(const Vector1<Real>& v)
{
    return Vector1<int>(static_cast<int>(std::floor(v.x)));
}

template<typename Real>
inline const Vector2<int> floor(const Vector2<Real>& v)
{
    return Vector2<int>(static_cast<int>(std::floor(v.x)),
        static_cast<int>(std::floor(v.y)));
}

template<typename Real>
inline const Vector3<int> floor(const Vector3<Real>& v)
{
    return Vector3<int>(static_cast<int>(std::floor(v.x)),
        static_cast<int>(std::floor(v.y)),
        static_cast<int>(std::floor(v.z)));
}

template<typename Real>
inline const Vector1<int> truncate(const Vector1<Real>& v)
{
    return Vector1<int>(static_cast<int>(v.x));
}

template<typename Real>
inline const Vector2<int> truncate(const Vector2<Real>& v)
{
    return Vector2<int>(static_cast<int>(v.x),
        static_cast<int>(v.y));
}

template<typename Real>
inline const Vector3<int> truncate(const Vector3<Real>& v)
{
    return Vector3<int>(static_cast<int>(v.x),
        static_cast<int>(v.y),
        static_cast<int>(v.z));
}

template<typename Real>
inline Real inverse(const Real a)
{
    return static_cast<Real>(1) / a;
}

template<typename Real>
inline const Vector2<Real> inverse(const Vector2<Real>& v)
{
    return Vector2<Real>(static_cast<Real>(1) / v.x, static_cast<Real>(1) / v.y);
}

template<typename Real>
inline const Vector3<Real> inverse(const Vector3<Real>& v)
{
    return Vector3<Real>(static_cast<Real>(1) / v.x, static_cast<Real>(1) / v.y, static_cast<Real>(1) / v.z);
}

template<Dimension dimension, typename T>
struct OnesHelper {
    static typename VectorTypeHelper<dimension, T>::Type get();
};

template<typename T>
struct OnesHelper<One, T> {
    static typename VectorTypeHelper<One, T>::Type get() { return typename VectorTypeHelper<One, T>::Type(static_cast<T>(1)); }
};

template<typename T>
struct OnesHelper<Two, T> {
    static typename VectorTypeHelper<Two, T>::Type get() { return typename VectorTypeHelper<Two, T>::Type(static_cast<T>(1), static_cast<T>(1)); }
};

template<typename T>
struct OnesHelper<Three, T> {
    static typename VectorTypeHelper<Three, T>::Type get() { return typename VectorTypeHelper<Three, T>::Type(static_cast<T>(1), static_cast<T>(1), static_cast<T>(1)); }
};

template<Dimension dimension, typename T>
inline typename VectorTypeHelper<dimension, T>::Type ones()
{
    return OnesHelper<dimension, T>::get();
}


template<typename T>
struct ScalarType {
};

template<typename T>
struct ScalarType<Vector1<T> > {
    typedef T Type;
};

template<typename T>
struct ScalarType<Vector2<T> > {
    typedef T Type;
};

template<typename T>
struct ScalarType<Vector3<T> > {
    typedef T Type;
};

inline Int3 remainder(const Int3& v1, const Int3& v2)
{
    return Int3(v1.x % v2.x, v1.y % v2.y, v1.z % v2.z);
}

inline FP dist(const FP3& v1, const FP3& v2)
{
    return (v1 - v2).norm();
}

inline const FP3 operator * (const FP3& v1, const Int3& v2)
{
    return FP3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

inline const FP3 operator * (const Int3& v1, const FP3& v2)
{
    return FP3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

inline FP sqr(const FP3& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline FP SP(const FP3& v1, const FP3& v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

__forceinline FP dot(const FP3& v1, const FP3& v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

__forceinline FP3 cross(const FP3& v1, const FP3& v2)
{
    return FP3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x);
}

inline const FP3 VP(const FP3& v1, const FP3& v2)
{
    return FP3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x);
}

typedef FP Vector3<FP>::* MemberOfFP3;
typedef complexFP Vector3<complexFP>::* MemberOfComplexFP3;


} // namespace pfc
