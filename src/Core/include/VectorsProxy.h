#pragma once

#include "Dimension.h"
#include "FP.h"
#include "Vectors.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <functional>

using namespace std;

namespace pfc {

    // Vector of 1 component of arithmetic type T (T = int, double, etc.)
    // with access by .x, and [], provides basic arithmetic operations
    // This is essentially a wrapper around T with interface compatible with Vector2 and Vector3
    template <typename T>
    struct Vector1Proxy
    {
        reference_wrapper<T> x;

        Vector1Proxy(T& _x) : x(_x)
        {}

        Vector1Proxy(Vector1<T>& other) : x(other.x)
        {}

        Vector1Proxy(const Vector1Proxy<T>& other) : x(other.x)
        {}

        inline T operator[](int) const
        {
            return x.get();
        }

        inline T& operator[](int)
        {
            return x.get();
        }

        inline T volume() const
        {
            return x.get();
        }

        inline T norm() const
        {
            T x_ = x.get();
            return x_ * ((x_ > 0) ? static_cast<T>(1) : static_cast<T>(-1));
        }

        inline T norm2() const
        {
            T x_ = x.get();
            return x_ * x_;
        }

        inline Vector1Proxy<T>& operator= (const Vector1<T>& v)
        {
            x.get() = v.x;
            return *this;
        }

        inline Vector1<T> toVector()
        {
            return Vector1<T>(x.get());
        }
    };

    template<typename T>
    inline const Vector1<T> operator + (const Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        return Vector1<T>(v1.x.get() + v2.x.get());
    }

    template<typename T>
    inline Vector1Proxy<T>& operator += (Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        v1.x.get() += v2.x.get();
        return v1;
    }

    template<typename T>
    inline const Vector1<T> operator - (const Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        return Vector1<T>(v1.x.get() - v2.x.get());
    }

    template<typename T>
    inline Vector1Proxy<T>& operator -= (Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        v1.x.get() -= v2.x.get();
        return v1;
    }

    template<typename T>
    inline const Vector1<T> operator * (const Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        return Vector1<T>(v1.x.get() * v2.x.get());
    }

    template<typename T>
    inline const Vector1<T> operator * (const Vector1Proxy<T>& v, T a)
    {
        return Vector1<T>(v.x.get() * a);
    }

    template<typename T>
    inline const Vector1<T> operator * (T a, const Vector1Proxy<T>& v)
    {
        return Vector1<T>(v.x.get() * a);
    }

    template<typename T>
    inline Vector1Proxy<T>& operator *= (Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        v1.x.get() *= v2.x.get();
        return v1;
    }

    template<typename T>
    inline Vector1Proxy<T>& operator *= (Vector1Proxy<T>& v, T a)
    {
        v.x.get() *= a;
        return v;
    }

    template<typename T>
    inline const Vector1<T> operator / (const Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        return Vector1<T>(v1.x.get() / v2.x.get());
    }

    template<typename T>
    inline Vector1Proxy<T>& operator /= (Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        v1.x.get() /= v2.x.get();
        return v1;
    }

    template<typename T>
    inline const Vector1<T> operator / (const Vector1Proxy<T>& v, T a)
    {
        return Vector1<T>(v.x.get() / a);
    }

    template<typename T>
    inline Vector1Proxy<T>& operator /= (Vector1Proxy<T>& v, T a)
    {
        v.x.get() /= a;
        return v;
    }

    template<typename T>
    inline bool operator == (const Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        return (v1.x.get() == v2.x.get());
    }

    template<typename T>
    inline bool operator != (const Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        return (v1.x.get() != v2.x.get());
    }

    template<typename T>
    inline bool operator < (const Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        return (v1.x.get() < v2.x.get());
    }

    template<typename T>
    inline bool operator <= (const Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        return (v1.x.get() <= v2.x.get());
    }

    template<typename T>
    inline bool operator > (const Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        return (v1.x.get() > v2.x.get());
    }

    template<typename T>
    inline bool operator >= (const Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        return (v1.x.get() >= v2.x.get());
    }

    template<typename T>
    inline std::ostream& operator<<(std::ostream& out, const Vector1Proxy<T>& v)
    {
        return out << "(" << v.x.get() << ")";
    }

    template<typename T>
    inline T dot(const Vector1Proxy<T>& v1, const Vector1Proxy<T>& v2)
    {
        return v1.x.get() * v2.x.get();
    }


    // Vector of 2 components of arithmetic type T (T = int, double, etc.)
    // with access by .x, .y, and [], provides basic arithmetic operations
    template <typename T>
    struct Vector2Proxy
    {
        reference_wrapper<T> x, y;

        Vector2Proxy(T& _x, T& _y) : x(_x), y(_y)
        {}

        Vector2Proxy(Vector2<T>& other) : x(other.x), y(other.y)
        {}

        Vector2Proxy(const Vector2Proxy<T>& other) : x(other.x), y(other.y)
        {}

        inline T operator[](int idx) const
        {
            return (*((reference_wrapper<T> *)this + idx)).get();
        }

        inline T& operator[](int idx)
        {
            return (*((reference_wrapper<T> *)this + idx)).get();
        }

        inline T volume() const
        {
            return x.get() * y.get();
        }

        inline T norm() const
        {
            T x_ = x.get(), y_ = y.get();
            return sqrt(x_ * x_ + y_ * y_);
        }

        inline T norm2() const
        {
            T x_ = x.get(), y_ = y.get();
            return x_ * x_ + y_ * y_;
        }

        inline Vector2Proxy<T>& operator= (const Vector2<T>& v)
        {
            x.get() = v.x;
            y.get() = v.y;
            return *this;
        }

        inline Vector2<T> toVector()
        {
            return Vector2<T>(x.get(), y.get());
        }
    };

    template<typename T>
    inline const Vector2<T> operator + (const Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        return Vector2<T>(v1.x.get() + v2.x.get(), v1.y.get() + v2.y.get());
    }

    template<typename T>
    inline Vector2Proxy<T>& operator += (Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        v1.x.get() += v2.x.get();
        v1.y.get() += v2.y.get();
        return v1;
    }

    template<typename T>
    inline const Vector2<T> operator - (const Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        return Vector2<T>(v1.x.get() - v2.x.get(), v1.y.get() - v2.y.get());
    }

    template<typename T>
    inline Vector2Proxy<T>& operator -= (Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        v1.x.get() -= v2.x.get();
        v1.y.get() -= v2.y.get();
        return v1;
    }

    template<typename T>
    inline const Vector2<T> operator * (const Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        return Vector2<T>(v1.x.get() * v2.x.get(), v1.y.get() * v2.y.get());
    }

    template<typename T>
    inline Vector2Proxy<T>& operator *= (Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        v1.x.get() *= v2.x.get();
        v1.y.get() *= v2.y.get();
        return v1;
    }

    template<typename T>
    inline const Vector2<T> operator * (const Vector2Proxy<T>& v, T a)
    {
        return Vector2<T>(v.x.get() * a, v.y.get() * a);
    }

    template<typename T>
    inline const Vector2<T> operator * (T a, const Vector2Proxy<T>& v)
    {
        return Vector2<T>(v.x.get() * a, v.y.get() * a);
    }

    template<typename T>
    inline Vector2Proxy<T>& operator *= (Vector2Proxy<T>& v, T a)
    {
        v.x.get() *= a;
        v.y.get() *= a;
        return v;
    }

    template<typename T>
    inline const Vector2<T> operator / (const Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        return Vector2<T>(v1.x.get() / v2.x.get(), v1.y.get() / v2.y.get());
    }

    template<typename T>
    inline Vector2Proxy<T>& operator /= (Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        v1.x.get() /= v2.x.get();
        v1.y.get() /= v2.y.get();
        return v1;
    }

    template<typename T>
    inline const Vector2<T> operator / (const Vector2Proxy<T>& v, T a)
    {
        return Vector2<T>(v.x.get() / a, v.y.get() / a);
    }

    template<typename T>
    inline Vector2Proxy<T>& operator /= (Vector2Proxy<T>& v, T a)
    {
        v.x.get() /= a;
        v.y.get() /= a;
        return v;
    }

    template<typename T>
    inline bool operator == (const Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        return (v1.x.get() == v2.x.get()) && (v1.y.get() == v2.y.get());
    }

    template<typename T>
    inline bool operator != (const Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        return (v1.x.get() != v2.x.get()) || (v1.y.get() != v2.y.get());
    }

    template<typename T>
    inline bool operator < (const Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        return (v1.x.get() < v2.x.get()) && (v1.y.get() < v2.y.get());
    }

    template<typename T>
    inline bool operator <= (const Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        return (v1.x.get() <= v2.x.get()) && (v1.y.get() <= v2.y.get());
    }

    template<typename T>
    inline bool operator > (const Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        return (v1.x.get() > v2.x.get()) && (v1.y.get() > v2.y.get());
    }

    template<typename T>
    inline bool operator >= (const Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        return (v1.x.get() >= v2.x.get()) && (v1.y.get() >= v2.y.get());
    }

    template<typename T>
    inline std::ostream& operator<<(std::ostream& out, const Vector2Proxy<T>& v)
    {
        return out << "(" << v.x.get() << ", " << v.y.get() << ")";
    }

    template<typename T>
    inline T dot(const Vector2Proxy<T>& v1, const Vector2Proxy<T>& v2)
    {
        return v1.x.get() * v2.x.get() + v1.y.get() * v2.y.get();
    }


    // Vector of 3 components of arithmetic type T (T = int, double, etc.)
    // with access by .x, .y, .z and [], provides basic arithmetic operations
    template <typename T>
    struct Vector3Proxy
    {
        reference_wrapper<T> x, y, z;

        Vector3Proxy(T& _x, T& _y, T& _z) : x(_x), y(_y), z(_z)
        {}

        Vector3Proxy(Vector3<T>& other) : 
            x(other.x), y(other.y), z(other.z) {}

        Vector3Proxy(const Vector3Proxy<T>& other) : 
            x(other.x), y(other.y), z(other.z) {}


        inline T operator[](int idx) const
        {
            return (*((reference_wrapper<T> *)this + idx)).get();
        }

        inline T& operator[](int idx)
        {
            return (*((reference_wrapper<T> *)this + idx)).get();
        }

        inline T volume() const
        {
            return x.get() * y.get() * z.get();
        }

        inline T norm() const
        {
            T x_ = x.get(), y_ = y.get(), z_ = z.get();
            return sqrt(x_ * x_ + y_ * y_ + z_ * z_);
        }

        inline T norm2() const
        {
            T x_ = x.get(), y_ = y.get(), z_ = z.get();
            return x_ * x_ + y_ * y_ + z_ * z_;
        }

        inline Vector3Proxy<T>& operator= (const Vector3<T>& v)
        {
            x.get() = v.x;
            y.get() = v.y;
            z.get() = v.z;
            return *this;
        }

        inline Vector3<T> toVector()
        {
            return Vector3<T>(x.get(), y.get(), z.get());
        }
    };

    template<typename T>
    inline const Vector3<T> operator + (const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return Vector3<T>(v1.x.get() + v2.x.get(), v1.y.get() + v2.y.get(), v1.z.get() + v2.z.get());
    }

    template<typename T>
    inline Vector3Proxy<T>& operator += (Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        v1.x.get() += v2.x.get();
        v1.y.get() += v2.y.get();
        v1.z.get() += v2.z.get();
        return v1;
    }
    
    template<typename T>
    inline Vector3Proxy<T>& operator += (Vector3Proxy<T>& v1, const Vector3<T>& v2)
    {
        v1.x.get() += v2.x;
        v1.y.get() += v2.y;
        v1.z.get() += v2.z;
        return v1;
    }

    template<typename T>
    inline const Vector3<T> operator - (const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return Vector3<T>(v1.x.get() - v2.x.get(), v1.y.get() - v2.y.get(), v1.z.get() - v2.z.get());
    }

    template<typename T>
    inline Vector3Proxy<T>& operator -= (Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        v1.x.get() -= v2.x.get();
        v1.y.get() -= v2.y.get();
        v1.z.get() -= v2.z.get();
        return v1;
    }

    template<typename T>
    inline const Vector3<T> operator * (const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return Vector3<T>(v1.x.get() * v2.x.get(), v1.y.get() * v2.y.get(), v1.z.get() * v2.z.get());
    }

    template<typename T>
    inline Vector3Proxy<T>& operator *= (Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        v1.x.get() *= v2.x.get();
        v1.y.get() *= v2.y.get();
        v1.z.get() *= v2.z.get();
        return v1;
    }

    template<typename T>
    inline const Vector3<T> operator * (const Vector3Proxy<T>& v, T a)
    {
        return Vector3<T>(v.x.get() * a, v.y.get() * a, v.z.get() * a);
    }

    template<typename T>
    inline const Vector3<T> operator * (T a, const Vector3Proxy<T>& v2)
    {
        return Vector3<T>(a * v2.x.get(), a * v2.y.get(), a * v2.z.get());
    }

    template<typename T>
    inline Vector3Proxy<T>& operator *= (Vector3Proxy<T>& v, T a)
    {
        v.x.get() *= a;
        v.y.get() *= a;
        v.z.get() *= a;
        return v;
    }

    template<typename T>
    inline const Vector3<T> operator / (const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return Vector3<T>(v1.x.get() / v2.x.get(), v1.y.get() / v2.y.get(), v1.z.get() / v2.z.get());
    }

    template<typename T>
    inline Vector3Proxy<T>& operator /= (Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        v1.x.get() /= v2.x.get();
        v1.y.get() /= v2.y.get();
        v1.z.get() /= v2.z.get();
        return v1;
    }

    template<typename T>
    inline const Vector3<T> operator / (const Vector3Proxy<T>& v, T a)
    {
        return Vector3<T>(v.x.get() / a, v.y.get() / a, v.z.get() / a);
    }

    template<typename T>
    inline Vector3Proxy<T>& operator /= (Vector3Proxy<T>& v, T a)
    {
        v.x.get() /= a;
        v.y.get() /= a;
        v.z.get() /= a;
        return v;
    }

    template<typename T>
    inline bool operator == (const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return (v1.x.get() == v2.x.get()) && (v1.y.get() == v2.y.get()) && (v1.z.get() == v2.z.get());
    }

    template<typename T>
    inline bool operator != (const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return (v1.x.get() != v2.x.get()) || (v1.y.get() != v2.y.get()) || (v1.z.get() != v2.z.get());
    }

    template<typename T>
    inline bool operator < (const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return (v1.x.get() < v2.x.get()) && (v1.y.get() < v2.y.get()) && (v1.z.get() < v2.z.get());
    }

    template<typename T>
    inline bool operator <= (const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return (v1.x.get() <= v2.x.get()) && (v1.y.get() <= v2.y.get()) && (v1.z.get() <= v2.z.get());
    }

    template<typename T>
    inline bool operator > (const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return (v1.x.get() > v2.x.get()) && (v1.y.get() > v2.y.get()) && (v1.z.get() > v2.z.get());
    }

    template<typename T>
    inline bool operator >= (const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return (v1.x.get() >= v2.x.get()) && (v1.y.get() >= v2.y.get()) && (v1.z.get() >= v2.z.get());
    }

    template<typename T>
    inline std::ostream& operator<<(std::ostream& out, const Vector3Proxy<T>& v)
    {
        return out << "(" << v.x.get() << ", " << v.y.get() << ", " << v.z.get() << ")";
    }

    template<typename T>
    inline T dot(const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return v1.x.get() * v2.x.get() + v1.y.get() * v2.y.get() + v1.z.get() * v2.z.get();
    }

    template<typename T>
    inline const Vector3<T> cross(const Vector3Proxy<T>& v1, const Vector3Proxy<T>& v2)
    {
        return Vector3<T>(v1.y.get() * v2.z.get() - v1.z.get() * v2.y.get(), 
            v1.z.get() * v2.x.get() - v1.x.get() * v2.z.get(),
            v1.x.get() * v2.y.get() - v1.y.get() * v2.x.get());
    }


    template<Dimension dimension, typename T>
    struct ProxyVectorTypeHelper {
    };

    template<typename T>
    struct ProxyVectorTypeHelper<One, T> {
        typedef Vector1Proxy<T> Type;
    };

    template<typename T>
    struct ProxyVectorTypeHelper<Two, T> {
        typedef Vector2Proxy<T> Type;
    };

    template<typename T>
    struct ProxyVectorTypeHelper<Three, T> {
        typedef Vector3Proxy<T> Type;
    };


    template<int dimension, typename T>
    struct ProxyVectorIntTypeHelper {
    };

    template<typename T>
    struct ProxyVectorIntTypeHelper<1, T> {
        typedef Vector1Proxy<T> Type;
    };

    template<typename T>
    struct ProxyVectorIntTypeHelper<2, T> {
        typedef Vector2Proxy<T> Type;
    };

    template<typename T>
    struct ProxyVectorIntTypeHelper<3, T> {
        typedef Vector3Proxy<T> Type;
    };


    template<typename T>
    struct VectorDimensionHelper<Vector1Proxy<T> > {
        static const int dimension = 1;
    };

    template<typename T>
    struct VectorDimensionHelper<Vector2Proxy<T> > {
        static const int dimension = 2;
    };

    template<typename T>
    struct VectorDimensionHelper<Vector3Proxy<T> > {
        static const int dimension = 3;
    };

    typedef Vector3Proxy<int> Int3Proxy;
    typedef Vector3Proxy<FP> FP3Proxy;

    template<typename Real>
    inline const Vector1<int> floor(const Vector1Proxy<Real>& v)
    {
        return Vector1<int>(static_cast<int>(std::floor(v.x.get())));
    }

    template<typename Real>
    inline const Vector2<int> floor(const Vector2Proxy<Real>& v)
    {
        return Vector2<int>(static_cast<int>(std::floor(v.x.get())),
            static_cast<int>(std::floor(v.y.get())));
    }

    template<typename Real>
    inline const Vector3<int> floor(const Vector3Proxy<Real>& v)
    {
        return Vector3<int>(static_cast<int>(std::floor(v.x.get())),
            static_cast<int>(std::floor(v.y.get())),
            static_cast<int>(std::floor(v.z.get())));
    }

    template<typename Real>
    inline const Vector1<int> truncate(const Vector1Proxy<Real>& v)
    {
        return Vector1<int>(static_cast<int>(v.x.get()));
    }

    template<typename Real>
    inline const Vector2<int> truncate(const Vector2Proxy<Real>& v)
    {
        return Vector2<int>(static_cast<int>(v.x.get()),
            static_cast<int>(v.y.get()));
    }

    template<typename Real>
    inline const Vector3<int> truncate(const Vector3Proxy<Real>& v)
    {
        return Vector3<int>(static_cast<int>(v.x.get()),
            static_cast<int>(v.y.get()),
            static_cast<int>(v.z.get()));
    }

    template<typename Real>
    inline const Vector1<Real> inverse(const Vector1Proxy<Real>& v)
    {
        return Vector1<Real>(static_cast<Real>(1) / v.x.get());
    }

    template<typename Real>
    inline const Vector2<Real> inverse(const Vector2Proxy<Real>& v)
    {
        return Vector2<Real>(static_cast<Real>(1) / v.x.get(), 
            static_cast<Real>(1) / v.y.get());
    }

    template<typename Real>
    inline const Vector3<Real> inverse(const Vector3Proxy<Real>& v)
    {
        return Vector3<Real>(static_cast<Real>(1) / v.x.get(), 
            static_cast<Real>(1) / v.y.get(), 
            static_cast<Real>(1) / v.z.get());
    }

    template<typename T>
    struct ScalarType<Vector1Proxy<T>> {
        typedef T Type;
    };

    template<typename T>
    struct ScalarType<Vector2Proxy<T>> {
        typedef T Type;
    };

    template<typename T>
    struct ScalarType<Vector3Proxy<T>> {
        typedef T Type;
    };

    inline Int3 remainder(const Int3Proxy& v1, const Int3Proxy& v2)
    {
        return Int3(v1.x.get() % v2.x.get(), v1.y.get() % v2.y.get(), v1.z.get() % v2.z.get());
    }

    inline FP dist(const FP3Proxy& v1, const FP3Proxy& v2)
    {
        return (v1 - v2).norm();
    }

} // namespace pfc
