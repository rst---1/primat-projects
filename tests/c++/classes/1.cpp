/*
 * =====================================================================================
 *
 *       Filename:  1.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06.06.2013 10:30:13
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
// #include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "../../../../projects/common/common.h"
#include <math.h>
#include <fstream>
#include <iostream>


// #include <vector>

// class A
// {
//     struct B {int b;};
//     public:
//         std::vector<B> b;
// };
// 
// enum class foo : size_t
// {
//     a,
//     b,
//     c,
//     d,
//     e,
//     j
// };

// struct 
// {
//     int a;
//     double b;
// } mmm()
// {
//     struct 
// {
//     int a;
//     double b;
// } M;
//     return M;
//     };

// auto bar (int a) -> decltype (struct {int a;})
// {
//     int b;
//     return a + 1;
// };

// struct Super
// {
// };
// 
// struct Children : Super
// {
//     int a;
// };
// 
// template <typename T>
// int bar ()
// {
//     T a;
//     return 1; 
// };
// 
// struct foo
// {
//     private:
//     struct ret
//     {
//         int a;
//     };
//     public:
//     foo () {};
//     ret operator () () {ret r; r.a = 10; return r;};
// };
// 
// template <typename T>
// struct A
// {
//     T b;
//     T a;
//     int val;
// };
// 
// template <typename T>
// struct B
// {
//     A<int> b;
// };
// 
// #define FOR(iter, begin, end) for(size_t iter = begin; iter < end; ++iter)
// 
// // template<> bar<int>();
// //
// class AAA
// {
//     public:
//         AAA () {a = 10;};
//         void foo () {auto bar = [this] () {++a;}; bar();};
//         int a;
// };
// 
// int bazz (int in)
// {
//     return in;
// };
// 
// // template <typename T>
// // struct Trans
// // {
// //     T operator () (T in) {return in};
// // };
// // 
// // struct Param
// // {
// // 
// // };
// 
// template <typename T>
// T trans (T in)
// {
//     return in;
// };
// 
// template <typename T>
// struct Fint
// {
//     T i;
//     Fint(T in){i = in;};
//     // operator T () {i = in;};
// };
// 
// 
// void operator++ (Fint<int> f) {Fint<int> a = f;};
// 
// struct Param {};
// template <typename T>
// constexpr T operator== (const Param p, const T a) { return a;};
// #define param(a) Param()
// // int operator = (T in) { i = in; return in;};
// // auto trans [] (
// // 
// // void operator "" _p(const char* str, size_t size)
// // {
// // };
// 
// // #define param(a,b) b
// 
// struct bbb 
// {
//     bbb(int i) {a = i;};
//     int a;
// };
// 
// void fooaaa(bbb a)
// {
// };
// 
// class I {
//     public:
//         virtual void f() = 0;
//         virtual void g() = 0;
// };
// 
// class A : public I {
//     public:
//         void f() { std::cout << "A: вызываем метод f()" << std::endl; }
//         void g() { std::cout << "A: вызываем метод g()" << std::endl; }
// };
// 
// class B : public I {
//     public:
//         void f() { std::cout << "B: вызываем метод f()" << std::endl; }
//         void g() { std::cout << "B: вызываем метод g()" << std::endl; }
// };
// 
// class Delegat 
// {
//     public:
//     void set_d(I i) { m = i;};
//     void f() {m.f();};
//     void g() {m.g();};
//     private:
//     I m;
// };

// constexpr int FOR(long long int



// constexpr double power(double a, long int i)
// {
//     return (i > 0 ? a * power(a, i - 1) : 1);
// };
// 
// constexpr double to_zero(double a)
// {
//     return (a > 0 ? to_zero(a - 1.0) : 0);
// };
// 
// struct foo {int x; int y;};
// struct foo bar = {x : 1, y : 15};
// struct foo bar = {y : 15, x : 1};
// 
// 
// struct foo {int x; int y; int operator()(){return x + y;};};
// foo{x : 10, y : 15}();
// foo{y : 15, x : 10}();
// 

namespace prmt
{
    template <int spacedim>
        class Point
        {
            Point () {};
        };

    template <>
        class Point<1>
        {
            public:
                Point () : _x(0.0) {};
                Point (dbl x) { _x = x; };

                dbl& x() { return _x; };

            private:
                dbl _x;
        };

    template <>
        class Point<2>
        {
            public:
                Point () : _x(0.0), _y(0.0) {};
                Point (dbl x, dbl y) { _x = x; _y = y; };

                dbl& x() { return _x; };
                dbl& y() { return _y; };

                dbl x() const { return _x; };
                dbl y() const { return _y; };

            private:
                dbl _x;
                dbl _y;
        };

    template <>
        class Point<3>
        {
            public:
                Point () : _x(0.0), _y(0.0), _z(0.0) {};
                Point (dbl x, dbl y, dbl z) { _x = x; _y = y; _z = z;};

                dbl& x() { return _x; };
                dbl& y() { return _y; };
                dbl& z() { return _z; };

            private:
                dbl _x;
                dbl _y;
                dbl _z;
        };

    class Parameter
    {
        public:
            Parameter (dbl val) { set_val (val); };
            void set (dbl val) { set_val (val); };
            operator dbl () { return value; };
        private:
            void set_val (dbl val) 
            { 
                assert((val > 0.0 - 1e-10) and (val < 1.0 + 1e-10));
                // if ((val > 0.0 - 1e-10) and (val < 1.0 + 1e-10))
                value = val;
                // else
            };
            dbl value = 0.0;
    };

    template <int spacedim>
        class Curve
        {
            Curve () {};
        };

    template <>
        class Curve<1>
        {
            public:
                // Curve () : fp(0.0), lp(0.0) {};
                Curve (Point<1> first_point, Point<1> last_point) 
                {
                    fp = first_point;
                    lp = last_point;
                };

                Point<1>& first_point() { return fp; };
                Point<1>& last_point() { return lp; };
                Point<1> operator() (Parameter parameter) 
                { 
                    return (lp.x() - fp.x()) * parameter + fp.x();
                };
            private:
                Point<1> fp;
                Point<1> lp;
        };

    template <>
        class Curve<2>
        {
            public:
                // Curve (Point<2> first_point, Point<2> last_point) 
                // {
                //     fp = first_point;
                //     lp = last_point;
                // };

                // virtual Point<2>& first_point() { return fp; };
                // virtual Point<2>& last_point() { return lp; };
                virtual Point<2> operator() (Parameter parameter) = 0;
            protected:
                // Curve () : fp(0.0, 0.0), lp(0.0, 0.0) {};
                // Point<2> fp;
                // Point<2> lp;
        };

    class PiecewiseLinearCurve : Curve<2>
    {
        public:

            PiecewiseLinearCurve (vec<Point<2>>& points) : 
                point(points), parametric_point(points.size()) 
        {
            parametric_point[0] = 0.0;
            FOR(i, 1, point.size())
            {
                parametric_point[i] = parametric_point[i-1] +
                    sqrt(
                            pow(point[i-1].x() - point[i].x(), 2.0) +
                            pow(point[i-1].y() - point[i].y(), 2.0)
                        );
            };
            FOR(i, 1, point.size())
                parametric_point[i] /= parametric_point.back();
        }; 

            // virtual Point<2>& first_point() { return point.front(); };
            // virtual Point<2>& last_point() { return point.back(); };
            virtual Point<2> operator () (Parameter parameter) override
            {
                FOR(i, 1, point.size())
                {
                    if (
                            (parameter > parametric_point[i-1] - 1e-10) and
                            (parameter < parametric_point[i]   + 1e-10)
                       )
                    {
                        dbl t = 
                            (parameter - parametric_point[i-1]) /
                            (parametric_point[i] - parametric_point[i-1]);
                        printf("t=%f ", t);
                        return Point<2>(
                                point[i-1].x() * (1.0 - t) + point[i].x() * t,
                                point[i-1].y() * (1.0 - t) + point[i].y() * t);
                    };
                };
                return Point<2>(0.0,0.0);

            };

            // protected:
            const vec<Point<2>> point;

        protected:
            vec<dbl> parametric_point;
    };

    template<int spacedim>
        class LoopCondition 
        {
            public:
                LoopCondition (Point<spacedim> w, vec<Point<spacedim>> b) :
                    substitutiv(w),     
                    substitutable(b) {};

                Point<spacedim> substitutiv;
                vec<Point<spacedim>> substitutable;
        };
};

using namespace prmt;

struct Test
{
    int a[4];
    union 
    {
    int b;
    int c;
    };
};

class A
{
    public:
    int s;
    friend class B;
};

class B
{
    public:
        int b;
};


int main ()
{
    Test t = Test{1,2,3,4,{b : 5}};
    printf("%d %d\n", t.a[2], t.b);
    vec<int>({1,2,3});
    int bbb[4] = {0, 1, 2, bbb[0]};
    vec<prmt::Point<2>> vc;
    Point<3> p;
    goto label;
    LABEL(label);
    p.z() = 10.0;
    printf("%f\n", p.z());
    Parameter p1(0.2);
    Curve<1> cu(Point<1>(-1.0), Point<1>(-5.0));
    printf("%f %f %f\n", 
            cu.first_point().x(),
            cu.last_point().x(),
            cu(Parameter(0.5)).x());
    std::vector<Point<2>> v = {
        Point<2>(0.0, 0.0),
        Point<2>(1.0, 1.0),
        Point<2>(2.0, 0.0)};
    PiecewiseLinearCurve cu2(v);
    {
        std::ofstream of1("curve.gpd");
        of1;
        double a = 0.0;
        while (a <= 1.)
        {
            printf("a=%f ", a);
            Point<2> p = cu2(Parameter(a));
            of1 << p.x() << " " << p.y() << std::endl;
            a += 0.1;
        };
    };

    // printf("%f %f %f\n", 
    //         cu2.first_point().x(),
    //         cu2.last_point().x(),
    //         0.0);
    // int i = 2;
    // auto a = [&i](){printf("%d\n", i);};
    // a();
    // i = 3;
    // a();
    // const long int degree = 20;
    // printf("%lf\n", power(2.0, degree));
    // printf("%lf\n", to_zero(power(2.0, degree)));
    // return d(2.1, 1000000000);
    // printf("%ld\n", i);
    // Delegat delegat;
    // A a;
    // delegate.set_d(A);
    // delegate.f();

    // Fint<int> fint = 10; 
    // fint == 20;
    // mock m;
    // m << 20;
    // int b = bazz(param(a, 10));
    // int b = bazz(param("a") = 10);
    // int b = bazz("a"_p (10));
    // int b = bazz(trans(10));
    // int b = bazz(param(blabla) == 10);
    // int c = bazz("blabla" == 10);
    // int b = bazz(/*a*/ 10));
    // // int b = ({int a = 10;a;});
    // printf("%d\n", b);
    // fooaaa(bbb(10));
    // bbb toto = 10;

    // auto a = foo()();
    // a.a;
    // decltype(a) b;
    // A<B<A<int>>> c;
    // A<A<A<int>>> d;
    // d.a.a.val = 9;
    // d.b.val = 10;
    // d.a.b.val = 11;

    // auto test = [] () {
    //     struct AA
    //     {
    //         AA () {};
    //         int a;
    //     };
    //     AA a;
    //     a.a = 10;
    //     return AA();
    // };

    // auto aaa = test ();
    // decltype(aaa) bbb;
    // aaa.a;
    // test().a;
    // FOR(i, 0, 10) int uuu = 0;
    // int abc = 10;
    // int dfj = 100;
    // auto fooo = [&abc] () {
    //     return struct {auto foooo = [&] () {++abc;};
    //     ++abc;
    //     return foooo;
    // };
    // printf("%d\n", abc);
    // fooo();
    // fooo()();
    // printf("%d\n", abc);
    // 
    // struct { int a; } bub[2];
    // c.a.b.a;
    // A a;
    // printf("%ld\n", sizeof(foo));
    // foo f;
    // Super s;
    // Children ch;
    // ch.a = 10;
    // s = ch;
    // s.a;
    // bar<struct{int a}>();
    // for (auto i : f)
    //     printf("%d\n", i);
    return 0;
};
