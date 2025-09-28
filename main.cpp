#include <vector>
#include <iostream>

const double EPSILON = 1e-9;

class Point
{
public:
    double x;
    double y;

    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}
    void print() const
    {
        std::cout << "(" << x << ", " << y << ")\n";
    }
    bool operator==(const Point &other) const
    {
        return x == other.x && y == other.y;
    }
    bool operator!=(const Point &other) const
    {
        return x != other.x || y != other.y;
    }
    Point operator-(const Point &other) const
    {
        return Point(x - other.x, y - other.y);
    }
};
class Point3
{
public:
    double x;
    double y;
    double z;

    Point3() : x(0), y(0), z(0) {}
    Point3(double x, double y, double z) : x(x), y(y), z(z) {}
    void print() const
    {
        std::cout << "(" << x << ", " << y << ", " << z << ")\n";
    }
    bool operator==(const Point3 &other) const
    {
        return x == other.x && y == other.y && z == other.z;
    }
    bool operator!=(const Point3 &other) const
    {
        return x != other.x || y != other.y || z != other.z;
    }
    Point3 operator-(const Point3 &other) const
    {
        return Point3(x - other.x, y - other.y, z - other.z);
    }
    Point3 CrossProduct(const Point3 &v) const
    {
        return {
        y * v.z - z * v.y,
        z * v.x - x * v.z,
        x * v.y - y * v.x
    };
    }
};
class Edge
{
public:
    Point start;
    Point end;
    Edge(Point start, Point end) : start(start), end(end) {}
    
    bool operator==(const Edge &other) const
    {
        return (start == other.start && end == other.end) ||
               (start == other.end && end == other.start);
    }
};

class Triangle
{
public:
    Point a;
    Point b;
    Point c;
    Triangle() : a(), b(), c() {}
    Triangle(Point a, Point b, Point c) : a(a), b(b), c(c) {}
    void print() const
    {
        std::cout << "Triangle vertices:\n";
        a.print();
        b.print();
        c.print();
    }
    Point centroid() const
    {
        double cx = (a.x + b.x + c.x) / 3;
        double cy = (a.y + b.y + c.y) / 3;
        return Point(cx, cy);
    }
    bool operator==(const Triangle &other) const
    {
        return (a == other.a || a == other.b || a == other.c) &&
               (b == other.a || b == other.b || b == other.c) &&
               (c == other.a || c == other.b || c == other.c);
    }
    double Radius() const
    {
        double ab = std::sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2));
        double bc = std::sqrt(std::pow(c.x - b.x, 2) + std::pow(c.y - b.y, 2));
        double ca = std::sqrt(std::pow(a.x - c.x, 2) + std::pow(a.y - c.y, 2));
        double s = (ab + bc + ca) / 2;
        double area = std::sqrt(s * (s - ab) * (s - bc) * (s - ca));
        return (ab * bc * ca) / (4 * area);
    }

    bool IsPointInTriangleCircle(const Point &p) const
    {
        Point3 pa(a.x, a.y, a.x * a.x + a.y * a.y);
        Point3 pb(b.x, b.y, b.x * b.x + b.y * b.y);
        Point3 pc(c.x, c.y, c.x * c.x + c.y * c.y);

        double parab = p.x * p.x + p.y * p.y;
        Point3 N = (pb - pa).CrossProduct(pc - pa); 
        N.print();
        double plane = (pa.z-(N.x * (p.x - pa.x) + N.y * (p.y - pa.y)) / N.z);
        return parab < plane;
    }
    double SignedArea() const
    {
        return a.x * (b.y - c.y) - a.y * (b.x - c.x) + 1 * (b.x * c.y - c.x * b.y);
    }
    bool IsPointOnEdge(const Point &p) const
    {
        double signedAreaAB = Triangle(a, b, p).SignedArea();
        double signedAreaBC = Triangle(b, c, p).SignedArea();
        double signedAreaCA = Triangle(c, a, p).SignedArea();

        return (std::abs(signedAreaAB) < EPSILON &&
                std::abs(signedAreaBC) < EPSILON &&
                std::abs(signedAreaCA) < EPSILON);
    }
    Edge GetEdgeByPoint(const Point &p) const
    {
        double signedAreaAB = Triangle(a, b, p).SignedArea();
        double signedAreaBC = Triangle(b, c, p).SignedArea();
        double signedAreaCA = Triangle(c, a, p).SignedArea();

        if (std::abs(signedAreaAB) < EPSILON)
            return Edge(a, b);
        else if (std::abs(signedAreaBC) < EPSILON)
            return Edge(b, c);
        else if (std::abs(signedAreaCA) < EPSILON)
            return Edge(c, a);
    }

    Point oppositeVertex(const Edge &edge){
        if (a != edge.start && a != edge.end)
            return a;
        else if (b != edge.start && b != edge.end)
            return b;
        else
            return c;
    }
    bool containsVertex(const Point &p) const{
        return p == a || p == b || p == c;
    }
    bool containsEdge(const Edge &edge) const
    {
        return (containsVertex(edge.start) && containsVertex(edge.end));
    }


};
std::vector<Triangle> FlipEdge(const Edge &edge, const Triangle &tri1, const Triangle &tri2)
{
    // Find the opposite vertices in each triangle
    Point opp1, opp2;
    if (tri1.a != edge.start && tri1.a != edge.end)
        opp1 = tri1.a;
    else if (tri1.b != edge.start && tri1.b != edge.end)
        opp1 = tri1.b;
    else
        opp1 = tri1.c;

    if (tri2.a != edge.start && tri2.a != edge.end)
        opp2 = tri2.a;
    else if (tri2.b != edge.start && tri2.b != edge.end)
        opp2 = tri2.b;
    else
        opp2 = tri2.c;

    // Create the new triangles formed by flipping the edge
    Triangle newTri1(opp1, opp2, edge.start);
    Triangle newTri2(opp1, opp2, edge.end);

    return {newTri1, newTri2};
}

class Triangulation
{
    std::vector<Triangle> triangles;

public:
    void addTriangle(const Triangle &triangle)
    {
        triangles.push_back(triangle);
    };

    std::vector<Triangle> getTriangles()
    {
        return triangles;
    };
    void print() const
    {
        for (const auto &triangle : triangles)
        {
            triangle.print();
        }
    }
    void DeleteTriangle(const Triangle &t){
        triangles.erase(std::remove(triangles.begin(), triangles.end(), t), triangles.end());
    }
    void InsertPoint(const Point &point)
    {
        // 1. Find Triangle containing the point
        // 2. if point lies on edge remove the edge and make edges to 4 points
        // 3. else Subdivide the triangle into 3 new triangles - make edges to 3 points
        // 4. legalize all required edges recursively

        // 1:
        Triangle t = FindTriangleByPoint(point, triangles[0]);//triangle we inserting into

        // if point on edge
        // find ajacent triangle
        // remove edge
        // make edges to 4 points
        // legalize them
        if (t.IsPointOnEdge(point))
        {
            Edge CommonEdgeOfTwoTriangles=t.GetEdgeByPoint(point);

            Triangle oppositeT = FindAjasentTriangleByEdge(CommonEdgeOfTwoTriangles, t);

            // Find the opposite vertex in the triangle
            Point oppositeVertex = oppositeT.oppositeVertex(CommonEdgeOfTwoTriangles);

            DeleteTriangle(t);
            DeleteTriangle(oppositeT);

            addTriangle(Triangle (t.a, t.b, point));
            addTriangle(Triangle(oppositeVertex,t.b,  point));
            addTriangle(Triangle(oppositeVertex, t.c, point));
            addTriangle(Triangle(t.c, t.a, point));

            LegalizeEdge(point, Edge(t.a, t.b));
            LegalizeEdge(point, Edge(t.b, oppositeVertex));
            LegalizeEdge(point, Edge(oppositeVertex,t.c));
            LegalizeEdge(point, Edge(t.c, t.a));
            return;
        }

        // else:

        // Remove the old triangle and add the new ones
        // triangles.erase(std::remove(triangles.begin(), triangles.end(), triangleInsertingInto), triangles.end());
        DeleteTriangle(t);

        addTriangle(Triangle(t.a, t.b, point));
        addTriangle(Triangle(t.b, t.c, point));
        addTriangle(Triangle(t.c, t.a, point));

        LegalizeEdge(point, Edge(t.a, t.b));
        LegalizeEdge(point, Edge(t.b, t.c));
        LegalizeEdge(point, Edge(t.c, t.a));
    }

    std::vector<Triangle> FindAjasentTrianglesByEdge(const Edge&edge)const
    {
        std::vector<Triangle> result;
        for (const auto &triangle : triangles)
        {
            if (triangle.containsEdge(edge))
            {
                result.push_back(triangle);
            }
        }
        return result;
    }
    bool OnExteriourFace(const Edge &edge){
        int count = 0;
        for (const auto &triangle : triangles)
        {
            if ((triangle.a == edge.start || triangle.a == edge.end) && (triangle.b == edge.start || triangle.b == edge.end) ||
                (triangle.a == edge.start || triangle.a == edge.end) && (triangle.c == edge.start || triangle.c == edge.end) ||
                (triangle.b == edge.start || triangle.b == edge.end) && (triangle.c == edge.start || triangle.c == edge.end))
            {
                count++;
            }
        }
        return count < 2; // If the edge is part of less than 2 triangles, it's on the exterior
    }
    
    Triangle FindAjasentTriangleByEdge(const Edge &edge, const Triangle &currentTriangle)
    {
        for (auto &triangle : triangles)
        {
            if (triangle == currentTriangle)
                continue;
            if (triangle.containsEdge(edge))
                return triangle;
        }
        for (auto &triangle : triangles)
            triangle.print();
        throw std::runtime_error("No adjacent triangle found for the given edge.");
    }

    Triangle FindTriangleByPoint(const Point &p, const Triangle &startingTriangle)
    {
        Triangle prevt = startingTriangle;
        Triangle t = startingTriangle;
        while (true)
        {
            // triangle t: a b c
            double signedAreaAB = Triangle(t.a, t.b, p).SignedArea();
            double signedAreaBC = Triangle(t.b, t.c, p).SignedArea();
            double signedAreaCA = Triangle(t.c, t.a, p).SignedArea();

            // Point is inside the triangle or on the edge
            if (signedAreaAB > 0 && signedAreaBC > 0 && signedAreaCA > 0 &&
                std::abs(signedAreaAB) < EPSILON &&
                std::abs(signedAreaBC) < EPSILON &&
                std::abs(signedAreaCA) < EPSILON)
            {

                return t;
            }

            prevt = t;
            // check if clockwise (outside of the triangle)
            if (signedAreaAB < 0)
                t = FindAjasentTriangleByEdge(Edge(t.a, t.b), t);
            else if (signedAreaBC < 0)
                t = FindAjasentTriangleByEdge(Edge(t.b, t.c), t);
            else if (signedAreaCA < 0)
                t = FindAjasentTriangleByEdge(Edge(t.c, t.a), t);

            if (t == prevt)
                return t;
        }
    }
    void LegalizeEdge(const Point &point, const Edge &edge)
{
    // Find the triangle that contains the edge
    std::vector<Triangle> adjTriangles = FindAjasentTrianglesByEdge(edge);
    if (adjTriangles.size() != 2)
    {
        // Edge is on the exterior, no need to legalize
        return;
    }

    Triangle oppositeT;
    Triangle currentT;
    if (adjTriangles[0].containsVertex(point))
    {
        oppositeT = adjTriangles[0];
        currentT = adjTriangles[1];
    }
    else
    {
        oppositeT = adjTriangles[1];
        currentT = adjTriangles[0];
    }

    // Find the opposite vertex in the triangle
    Point oppositeVertex = oppositeT.oppositeVertex(edge);

    // Check if the circumcircle of the triangle contains the point
    Triangle testTriangle(edge.start, edge.end, point);
    if (currentT.IsPointInTriangleCircle(oppositeVertex))
    {
        // Flip the edge
        DeleteTriangle(oppositeT);
        DeleteTriangle(currentT);

        std::vector<Triangle> outTriangles = FlipEdge(edge, oppositeT, currentT);
        for (const auto &tri : outTriangles)
        {
            addTriangle(tri);
        }

        // Recursively legalize the new edges
        LegalizeEdge(point, Edge(edge.start, oppositeVertex));
        LegalizeEdge(point, Edge(edge.end, oppositeVertex));
    }
}
void RemoveSuperTriangleTriangles(const Triangle &superTriangle)
{
    auto it = triangles.begin();
    while (it != triangles.end())
    {
        if (it->containsVertex(superTriangle.a) ||
            it->containsVertex(superTriangle.b) ||
            it->containsVertex(superTriangle.c))
        {
            it = triangles.erase(it);
        }
        else
        {
            ++it;
        }
    }
    
}

};

Triangle createSuperTriangle(std::vector<Point> points)
{
    double minX = points[0].x;
    double minY = points[0].y;
    double maxX = points[0].x;
    double maxY = points[0].y;

    for (const auto &point : points)
    {
        if (point.x < minX)
            minX = point.x;
        if (point.y < minY)
            minY = point.y;
        if (point.x > maxX)
            maxX = point.x;
        if (point.y > maxY)
            maxY = point.y;
    }

    double dx = maxX - minX;
    double dy = maxY - minY;
    double deltaMax = std::max(dx, dy);
    double midX = (minX + maxX) / 2;
    double midY = (minY + maxY) / 2;

    Point p1(midX - 20 * deltaMax, midY - deltaMax);
    Point p2(midX, midY + 20 * deltaMax);
    Point p3(midX + 20 * deltaMax, midY - deltaMax);

    return Triangle(p1, p2, p3);
}


// Example usage
int main() {
    try {
        
        
        // Create super triangle
        Triangle tr(Point(1,1),Point(3,3),Point(2,4));

        // std::cout<<"Tring Result:"<<tr.IsPointInTriangleCircle(Point(8,2))<<std::endl;
        // std::cout<<"Tring Result:"<<tr.IsPointInTriangleCircle(Point(1,2))<<std::endl;






        // // Initialize triangulation with super triangle
        // Triangulation triangulation;
        // triangulation.addTriangle(superTriangle);
        
        // // Insert all points
        // for (const auto& point : points) {
        //     triangulation.InsertPoint(point);
        // }
        
        // // Remove triangles containing super triangle vertices
        // triangulation.RemoveTrianglesWithSuperTriangleVertices(superTriangle);
        
        // // Print result
        // triangulation.print();
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    
    return 0;
}