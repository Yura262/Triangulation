#include <vector>
#include <iostream>


class Point{
public:
    int x;
    int y;


    Point(int x, int y): x(x), y(y) {}
    void print() const {
        std::cout << "(" << x << ", " << y << ")\n";
    }
    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }
};


class Edge{
    public:
    Point start;
    Point end;
    Edge(Point start, Point end): start(start), end(end) {}};

class Triangle{
public:
    Point a;
    Point b;
    Point c;
    Triangle(Point a, Point b, Point c): a(a), b(b), c(c) {}
    void print() const {
        std::cout << "Triangle vertices:\n";
        a.print();
        b.print();
        c.print();
    }
    Point centroid() const {
        int cx = (a.x + b.x + c.x) / 3;
        int cy = (a.y + b.y + c.y) / 3;
        return Point(cx, cy);
    }
    double Radius() const {
        double ab = std::sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2));
        double bc = std::sqrt(std::pow(c.x - b.x, 2) + std::pow(c.y - b.y, 2));
        double ca = std::sqrt(std::pow(a.x - c.x, 2) + std::pow(a.y - c.y, 2));
        double s = (ab + bc + ca) / 2;
        double area = std::sqrt(s * (s - ab) * (s - bc) * (s - ca));
        return (ab * bc * ca) / (4 * area);
    }

    bool IsPointInTriangle(const Point& p) const {
        // Barycentric Technique
        int alpha = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y)) /
                    ((b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y));
        int beta = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y)) /
                   ((c.y - a.y) * (b.x - c.x) + (a.x - c.x) * (b.y - c.y));
        int gamma = 1.0f - alpha - beta;

        return alpha >= 0 && beta >= 0 && gamma >= 0;
    }


};


class Triangulation{
    std::vector<Triangle> triangles;
public:
    void addTriangle(const Triangle& triangle) {
        triangles.push_back(triangle);
    };

    std::vector<Triangle> getTriangles() {
        return triangles;
    };
    void print() const {
        for (const auto& triangle : triangles) {
            triangle.print();
        }
    }

    void InsertPoint(const Point& point) {
        //1. Find Triangle containing the point
        //2. if point lies on edge remove the edge and make edges to 4 points //TODO------------------
        //3. else Subdivide the triangle into 3 new triangles - make edges to 3 points
        //4. legalize all required edges recursively

        std::vector<Edge> edgesToLegalize; 
        //1:
        for (const auto& triangle : triangles) {
            if (triangle.IsPointInTriangle(point)) {
                Triangle t1(triangle.a, triangle.b, point);
                Triangle t2(triangle.b, triangle.c, point);
                Triangle t3(triangle.c, triangle.a, point);
                
                // Remove the old triangle and add the new ones
                triangles.erase(std::remove(triangles.begin(), triangles.end(), triangle), triangles.end());
                triangles.push_back(t1);
                triangles.push_back(t2);
                triangles.push_back(t3);

                edgesToLegalize.push_back(Edge(triangle.a, triangle.b));
                edgesToLegalize.push_back(Edge(triangle.b, triangle.c));
                edgesToLegalize.push_back(Edge(triangle.c, triangle.a));


                break;
            }
        }
        //TODO------------------
        for (const auto& edge : edgesToLegalize) {
            if edge
            Triangle adjacentTriangle = FindAjasentTriangleByEdge(edge.start, edge.end);
            // Legalize edge if necessary
        }

    }
    
    Triangle FindAjasentTriangleByEdge(const Point& edgeStart, const Point& edgeEnd) {
            for (const auto& triangle : triangles) {
                if ((triangle.a == edgeStart || triangle.a == edgeEnd) && (triangle.b == edgeStart || triangle.b == edgeEnd) || 
                (triangle.a == edgeStart || triangle.a == edgeEnd) && (triangle.c == edgeStart || triangle.c == edgeEnd) || 
                (triangle.b == edgeStart || triangle.b == edgeEnd) && (triangle.c == edgeStart || triangle.c == edgeEnd)) 
                {
                return triangle;
                }
            }   
    }
    
};

Triangle createSuperTriangle(std::vector<Point> points) {
    int minX = points[0].x;
    int minY = points[0].y;
    int maxX = points[0].x;
    int maxY = points[0].y;

    for (const auto& point : points) {
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



int main() {
    std::vector<Point> points = {Point(0, 0), Point(1, 2), Point(2, 1), Point(3, 3)};
    Triangle superTriangle = createSuperTriangle(points);
    superTriangle.print();
    std::cout << "Centroid: ";
    superTriangle.centroid().print();
    std::cout << "Circumradius: " << superTriangle.Radius() << "\n";

    Triangulation triangulation;
    for (const auto& point : points) {
        triangulation.InsertPoint(point);
    }    
    triangulation.print();

    return 0;
}