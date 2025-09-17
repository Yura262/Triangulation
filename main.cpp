#include <vector>
#include <iostream>


class Point{
    int x;
    int y;

public:

    Point(int x, int y): x(x), y(y) {}
    void print() const {
        std::cout << "(" << x << ", " << y << ")\n";
    }

};

class Triangle{
    Point a;
    Point b;
    Point c;
public:
    Triangle(Point a, Point b, Point c): a(a), b(b), c(c) {}
    void print() const {
        std::cout << "Triangle vertices:\n";
        a.print();
        b.print();
        c.print();
    }
};


