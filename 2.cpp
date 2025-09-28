#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

const double EPSILON = 1e-9;

class Point
{
public:
    int x;
    int y;

    Point() : x(0), y(0) {}
    Point(int x, int y) : x(x), y(y) {}
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
    double x;  // Changed to double for better precision
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
    Point3 CrossProduct(const Point3 &p2) const  // Added const
    {
        return Point3(
            y * p2.z - z * p2.y,
            z * p2.x - x * p2.z,
            x * p2.y - y * p2.x);
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
        int cx = (a.x + b.x + c.x) / 3;
        int cy = (a.y + b.y + c.y) / 3;
        return Point(cx, cy);
    }
    bool operator==(const Triangle &other) const
    {
        // Check if triangles have the same vertices (in any order)
        int count = 0;
        if (a == other.a || a == other.b || a == other.c) count++;
        if (b == other.a || b == other.b || b == other.c) count++;
        if (c == other.a || c == other.b || c == other.c) count++;
        return count == 3;
    }
    double Radius() const
    {
        double ab = std::sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2));
        double bc = std::sqrt(std::pow(c.x - b.x, 2) + std::pow(c.y - b.y, 2));
        double ca = std::sqrt(std::pow(a.x - c.x, 2) + std::pow(a.y - c.y, 2));
        double s = (ab + bc + ca) / 2;
        double area = std::sqrt(s * (s - ab) * (s - bc) * (s - ca));
        if (area < EPSILON) return 0;  // Degenerate triangle
        return (ab * bc * ca) / (4 * area);
    }
    
    bool IsPointInCircumcircle(const Point &p) const
    {
        // Using determinant method for in-circle test
        double ax = a.x - p.x;
        double ay = a.y - p.y;
        double bx = b.x - p.x;
        double by = b.y - p.y;
        double cx = c.x - p.x;
        double cy = c.y - p.y;
        
        double det = (ax * ax + ay * ay) * (bx * cy - cx * by) -
                     (bx * bx + by * by) * (ax * cy - cx * ay) +
                     (cx * cx + cy * cy) * (ax * by - bx * ay);
        
        return det > EPSILON;
    }
    
    double SignedArea() const
    {
        return 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
    }
    
    bool IsPointInside(const Point &p) const
    {
        double area = std::abs(SignedArea());
        double area1 = std::abs(Triangle(p, b, c).SignedArea());
        double area2 = std::abs(Triangle(a, p, c).SignedArea());
        double area3 = std::abs(Triangle(a, b, p).SignedArea());
        
        return std::abs(area - (area1 + area2 + area3)) < EPSILON;
    }
    
    bool IsPointOnEdge(const Point &p) const
    {
        // Check if point is on edge AB
        double crossAB = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
        if (std::abs(crossAB) < EPSILON) {
            double dotAB = (p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y);
            double lenAB = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);
            if (dotAB >= -EPSILON && dotAB <= lenAB + EPSILON) return true;
        }
        
        // Check if point is on edge BC
        double crossBC = (c.x - b.x) * (p.y - b.y) - (c.y - b.y) * (p.x - b.x);
        if (std::abs(crossBC) < EPSILON) {
            double dotBC = (p.x - b.x) * (c.x - b.x) + (p.y - b.y) * (c.y - b.y);
            double lenBC = (c.x - b.x) * (c.x - b.x) + (c.y - b.y) * (c.y - b.y);
            if (dotBC >= -EPSILON && dotBC <= lenBC + EPSILON) return true;
        }
        
        // Check if point is on edge CA
        double crossCA = (a.x - c.x) * (p.y - c.y) - (a.y - c.y) * (p.x - c.x);
        if (std::abs(crossCA) < EPSILON) {
            double dotCA = (p.x - c.x) * (a.x - c.x) + (p.y - c.y) * (a.y - c.y);
            double lenCA = (a.x - c.x) * (a.x - c.x) + (a.y - c.y) * (a.y - c.y);
            if (dotCA >= -EPSILON && dotCA <= lenCA + EPSILON) return true;
        }
        
        return false;
    }
    
    Edge GetEdgeContainingPoint(const Point &p) const
    {
        // Check if point is on edge AB
        double crossAB = (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
        if (std::abs(crossAB) < EPSILON) {
            double dotAB = (p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y);
            double lenAB = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);
            if (dotAB >= -EPSILON && dotAB <= lenAB + EPSILON) 
                return Edge(a, b);
        }
        
        // Check if point is on edge BC
        double crossBC = (c.x - b.x) * (p.y - b.y) - (c.y - b.y) * (p.x - b.x);
        if (std::abs(crossBC) < EPSILON) {
            double dotBC = (p.x - b.x) * (c.x - b.x) + (p.y - b.y) * (c.y - b.y);
            double lenBC = (c.x - b.x) * (c.x - b.x) + (c.y - b.y) * (c.y - b.y);
            if (dotBC >= -EPSILON && dotBC <= lenBC + EPSILON) 
                return Edge(b, c);
        }
        
        // Must be on edge CA
        return Edge(c, a);
    }

    Point oppositeVertex(const Edge &edge) const
    {
        if (a != edge.start && a != edge.end)
            return a;
        else if (b != edge.start && b != edge.end)
            return b;
        else
            return c;
    }
    
    bool containsVertex(const Point &p) const
    {
        return p == a || p == b || p == c;
    }
    
    bool containsEdge(const Edge &edge) const
    {
        return (containsVertex(edge.start) && containsVertex(edge.end));
    }
};

class Triangulation
{
    std::vector<Triangle> triangles;

public:
    void addTriangle(const Triangle &triangle)
    {
        triangles.push_back(triangle);
    }

    std::vector<Triangle> getTriangles() const
    {
        return triangles;
    }
    
    void print() const
    {
        std::cout << "Triangulation with " << triangles.size() << " triangles:\n";
        for (size_t i = 0; i < triangles.size(); i++)
        {
            std::cout << "Triangle " << i << ":\n";
            triangles[i].print();
        }
    }
    
    void DeleteTriangle(const Triangle &t)
    {
        triangles.erase(std::remove(triangles.begin(), triangles.end(), t), triangles.end());
    }
    
    Triangle* FindTriangleContainingPoint(const Point &p)
    {
        for (auto &triangle : triangles)
        {
            if (triangle.IsPointInside(p))
                return &triangle;
        }
        return nullptr;
    }
    
    std::vector<Triangle> FindTrianglesWithEdge(const Edge &edge) const
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
    
    Triangle* FindAdjacentTriangle(const Edge &edge, const Triangle &currentTriangle)
    {
        for (auto &triangle : triangles)
        {
            if (triangle == currentTriangle)
                continue;
            if (triangle.containsEdge(edge))
                return &triangle;
        }
        return nullptr;
    }
    
    void InsertPoint(const Point &point)
    {
        // Find triangle containing the point
        Triangle* containingTriangle = FindTriangleContainingPoint(point);
        if (!containingTriangle)
        {
            std::cerr << "Point is outside triangulation\n";
            return;
        }
        
        Triangle t = *containingTriangle;  // Make a copy
        
        // Check if point is on an edge
        if (t.IsPointOnEdge(point))
        {
            Edge sharedEdge = t.GetEdgeContainingPoint(point);
            Triangle* adjacentTriangle = FindAdjacentTriangle(sharedEdge, t);
            
            if (adjacentTriangle)
            {
                Triangle adjT = *adjacentTriangle;  // Make a copy
                Point oppositeInT = t.oppositeVertex(sharedEdge);
                Point oppositeInAdjT = adjT.oppositeVertex(sharedEdge);
                
                // Remove old triangles
                DeleteTriangle(t);
                DeleteTriangle(adjT);
                
                // Create 4 new triangles
                addTriangle(Triangle(sharedEdge.start, point, oppositeInT));
                addTriangle(Triangle(point, sharedEdge.end, oppositeInT));
                addTriangle(Triangle(sharedEdge.start, point, oppositeInAdjT));
                addTriangle(Triangle(point, sharedEdge.end, oppositeInAdjT));
                
                // Legalize edges
                LegalizeEdge(point, Edge(sharedEdge.start, oppositeInT));
                LegalizeEdge(point, Edge(sharedEdge.end, oppositeInT));
                LegalizeEdge(point, Edge(sharedEdge.start, oppositeInAdjT));
                LegalizeEdge(point, Edge(sharedEdge.end, oppositeInAdjT));
            }
            return;
        }
        
        // Point is inside triangle - split into 3
        DeleteTriangle(t);
        
        addTriangle(Triangle(t.a, t.b, point));
        addTriangle(Triangle(t.b, t.c, point));
        addTriangle(Triangle(t.c, t.a, point));
        
        // Legalize edges
        LegalizeEdge(point, Edge(t.a, t.b));
        LegalizeEdge(point, Edge(t.b, t.c));
        LegalizeEdge(point, Edge(t.c, t.a));
    }
    
    void LegalizeEdge(const Point &point, const Edge &edge)
    {
        std::vector<Triangle> adjacentTriangles = FindTrianglesWithEdge(edge);
        
        if (adjacentTriangles.size() != 2)
            return;  // Edge is on boundary
        
        Triangle* tri1 = nullptr;
        Triangle* tri2 = nullptr;
        
        // Find which triangle contains the new point
        for (auto &tri : adjacentTriangles)
        {
            if (tri.containsVertex(point))
                tri1 = &tri;
            else
                tri2 = &tri;
        }
        
        if (!tri1 || !tri2)
            return;
        
        Point oppositeVertex = tri2->oppositeVertex(edge);
        
        // Check Delaunay condition
        if (tri1->IsPointInCircumcircle(oppositeVertex))
        {
            // Need to flip edge
            Triangle newTri1 = *tri1;  // Copy
            Triangle newTri2 = *tri2;  // Copy
            
            DeleteTriangle(newTri1);
            DeleteTriangle(newTri2);
            
            // Create new triangles with flipped edge
            addTriangle(Triangle(point, edge.start, oppositeVertex));
            addTriangle(Triangle(point, oppositeVertex, edge.end));
            
            // Recursively legalize new edges
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

Triangle createSuperTriangle(const std::vector<Point> &points)
{
    if (points.empty())
        throw std::runtime_error("Cannot create super triangle from empty point set");
    
    int minX = points[0].x;
    int minY = points[0].y;
    int maxX = points[0].x;
    int maxY = points[0].y;

    for (const auto &point : points)
    {
        if (point.x < minX) minX = point.x;
        if (point.y < minY) minY = point.y;
        if (point.x > maxX) maxX = point.x;
        if (point.y > maxY) maxY = point.y;
    }

    int dx = maxX - minX;
    int dy = maxY - minY;
    int deltaMax = std::max(dx, dy);
    int midX = (minX + maxX) / 2;
    int midY = (minY + maxY) / 2;

    Point p1(midX - 20 * deltaMax, midY - deltaMax);
    Point p2(midX, midY + 20 * deltaMax);
    Point p3(midX + 20 * deltaMax, midY - deltaMax);

    return Triangle(p1, p2, p3);
}

// Example usage
int main() {
    try {
        std::vector<Point> points = {
            Point(0, 0), 
            Point(4, 0), 
            Point(2, 3),
            Point(1, 1),
            Point(3, 1)
        };
        
        // Create super triangle
        Triangle superTriangle = createSuperTriangle(points);
        
        // Initialize triangulation with super triangle
        Triangulation triangulation;
        triangulation.addTriangle(superTriangle);
        
        // Insert all points
        for (const auto& point : points) {
            triangulation.InsertPoint(point);
        }
        
        // Remove triangles containing super triangle vertices
        triangulation.RemoveTrianglesWithSuperTriangleVertices(superTriangle);
        
        // Print result
        triangulation.print();
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    
    return 0;
}