#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <stack>
#include <map>
#include <stdexcept>

const double EPSILON = 1e-9;
constexpr double PI = 3.14159265358979323846;

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
    Point operator+(const Point& other) const {
        return Point{x + other.x, y + other.y};
    }
    Point operator/(const double other) const {
        return Point{x / other, y / other};
    }
    double dist2(const Point &p2) const
    {
        double dx = x - p2.x;
        double dy = y - p2.y;
        return dx * dx + dy * dy;
    }
    friend std::ostream &operator<<(std::ostream &os, const Point &obj)
    {
        os << "(" << obj.x << ", " << obj.y << ")";
        return os;
    }
    bool operator<(const Point &other) const
    {
        if (std::abs(x - other.x) > EPSILON)
            return x < other.x;
        return y < other.y;
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
            x * v.y - y * v.x};
    }
};
class Edge
{
public:
    Point start;
    Point end;
    Edge() : start(), end() {}
    Edge(Point s, Point e) : start(s), end(e)
    {
        // normalize ordering so Edge(a,b) and Edge(b,a) map to the same key
        if (end < start) std::swap(start, end);
    }

    bool operator==(const Edge &other) const
    {
        return (start == other.start && end == other.end);
    }
    void print() const
    {
        start.print();
        end.print();
    }
    bool operator<(const Edge &other) const
    {
        if (start < other.start)
            return true;
        if (other.start < start)
            return false;
        return end < other.end;
    }
};
static double orient(const Point &a, const Point &b, const Point &c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}
// robust in-circle test using determinant and orientation
static long double inCircleDeterminant(const Point &a, const Point &b, const Point &c, const Point &p)
{
    auto sq = [](const Point &q) -> long double
    { return (long double)q.x * q.x + (long double)q.y * q.y; };

    long double ax = a.x, ay = a.y, bx = b.x, by = b.y, cx = c.x, cy = c.y, px = p.x, py = p.y;
    long double a2 = ax * ax + ay * ay;
    long double b2 = bx * bx + by * by;
    long double c2 = cx * cx + cy * cy;
    long double p2 = px * px + py * py;

    long double m00 = ax - px;
    long double m01 = ay - py;
    long double m02 = a2 - p2;

    long double m10 = bx - px;
    long double m11 = by - py;
    long double m12 = b2 - p2;

    long double m20 = cx - px;
    long double m21 = cy - py;
    long double m22 = c2 - p2;

    long double det = m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) + m02 * (m10 * m21 - m11 * m20);
    return det;
}

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
    friend std::ostream &operator<<(std::ostream &os, const Triangle &obj)
    {
        os << "Triangle: " << obj.a << " " << obj.b << " " << obj.c << ", \n";
        return os;
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
        double area = std::sqrt(std::max(0.0, s * (s - ab) * (s - bc) * (s - ca)));
        if (area < EPSILON) return 0.0;
        return (ab * bc * ca) / (4 * area);
    }
    Point circumcircleCenter() const
    {
        Point A = a;
        Point B = b;
        Point C = c;
        double d = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
        if (std::abs(d) < 1e-12)
        {
            return Point(0, 0); // collinear case
        }
        double ux = ((A.x * A.x + A.y * A.y) * (B.y - C.y) +
                     (B.x * B.x + B.y * B.y) * (C.y - A.y) +
                     (C.x * C.x + C.y * C.y) * (A.y - B.y)) /
                    d;

        double uy = ((A.x * A.x + A.y * A.y) * (C.x - B.x) +
                     (B.x * B.x + B.y * B.y) * (A.x - C.x) +
                     (C.x * C.x + C.y * C.y) * (B.x - A.x)) /
                    d;

        Point center{ux, uy};

        return center;
    }
    double radiusOfCircle() const
    {
        // returns squared radius
        return circumcircleCenter().dist2(a);
    }

    double SignedArea() const
    {
        return (c.x - a.x) * (b.y - a.y) - (b.x - a.x) * (c.y - a.y);
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

        throw std::runtime_error("Point is not on any edge of the triangle");
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
        // check AB
        double areaABP = orient(a, b, p);
        if (std::abs(areaABP) < EPSILON)
        {
            // check bounding box
            if ((p.x >= std::min(a.x, b.x) - EPSILON && p.x <= std::max(a.x, b.x) + EPSILON) &&
                (p.y >= std::min(a.y, b.y) - EPSILON && p.y <= std::max(a.y, b.y) + EPSILON))
                return true;
        }
        // BC
        double areaBCP = orient(b, c, p);
        if (std::abs(areaBCP) < EPSILON)
        {
            if ((p.x >= std::min(b.x, c.x) - EPSILON && p.x <= std::max(b.x, c.x) + EPSILON) &&
                (p.y >= std::min(b.y, c.y) - EPSILON && p.y <= std::max(b.y, c.y) + EPSILON))
                return true;
        }
        // CA
        double areaCAP = orient(c, a, p);
        if (std::abs(areaCAP) < EPSILON)
        {
            if ((p.x >= std::min(c.x, a.x) - EPSILON && p.x <= std::max(c.x, a.x) + EPSILON) &&
                (p.y >= std::min(c.y, a.y) - EPSILON && p.y <= std::max(c.y, a.y) + EPSILON))
                return true;
        }
        return false;
    }
    bool IsPointInTriangleCircle(const Point &p) const
    {
        double orientation = orient(a, b, c);
        long double det = inCircleDeterminant(a, b, c, p);

        const long double EPS_IN = 1e-12L;
        if (std::fabsl(orientation) < EPS_IN)
        {
            // Degenerate triangle (almost collinear) -> treat as not containing point
            return false;
        }
        if (orientation > 0)
        {
            return det > EPS_IN;
        }
        else
        {
            return det < -EPS_IN;
        }
    }
    std::vector<Edge> edges() const
    {
        std::vector<Edge> edges;
        edges.push_back(Edge(a, b));
        edges.push_back(Edge(b, c));
        edges.push_back(Edge(c, a));
        return edges;
    }
    double area() const
    {
        return std::abs((a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) / 2.0);
    }

    double angleFromSides(double a_len, double b_len, double c_len) const
    {
        // angle opposite to side 'a_len'
        double cosAngle = (b_len * b_len + c_len * c_len - a_len * a_len) / (2 * b_len * c_len);
        cosAngle = std::clamp(cosAngle, -1.0, 1.0);
        return std::acos(cosAngle) * 180.0 / PI;
    }

    // Function to calculate the minimum angle of a triangle
    double minAngleDeg() const
    {
        // use real side lengths (not squared)
        double aa = std::sqrt(b.dist2(c)); // side opposite A
        double bb = std::sqrt(a.dist2(c)); // side opposite B
        double cc = std::sqrt(a.dist2(b)); // side opposite C

        double angleA = angleFromSides(aa, bb, cc);
        double angleB = angleFromSides(bb, aa, cc);
        double angleC = angleFromSides(cc, aa, bb);

        return std::min({angleA, angleB, angleC});
    }

    Edge longestEdge() const
    {
        double edgeAB = a.dist2(b); // squared distances OK for comparison
        double edgeBC = b.dist2(c);
        double edgeCA = c.dist2(a);
        Point edgeStart;
        Point edgeEnd;
        if (edgeAB >= edgeBC && edgeAB >= edgeCA) {
            edgeStart = a;
            edgeEnd = b;
        } else if (edgeBC >= edgeAB && edgeBC >= edgeCA) {
            edgeStart = b;
            edgeEnd = c;
        } else {
            edgeStart = c;
            edgeEnd = a;
        }
        return Edge(edgeStart,edgeEnd);
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
    double min_angle_deg = 50;
    double max_area = 5;

public:
    void setMinAngle(double angle) { min_angle_deg = angle; }
    void setMaxArea(double area) { max_area = area; }
    double getMinAngle() const { return min_angle_deg; }
    double getMaxArea() const { return max_area; }
    void addTriangle(const Triangle &triangle)
    {
        triangles.push_back(triangle);
    };
    std::vector<Point> getPoints() const
    {
        std::vector<Point> points;
        for (const auto &t : triangles)
        {
            points.push_back(t.a);
            points.push_back(t.b);
            points.push_back(t.c);
        }
        return points;
    }
    std::vector<Triangle> getTriangles() const
    {
        return triangles;
    };
    void print() const
    {
        std::cout << "Printing triangulation\n";
        for (const auto &triangle : triangles)
        {
            triangle.print();
        }
        std::cout << std::endl;
    }
    void DeleteTriangle(const Triangle &t)
    {
        triangles.erase(std::remove(triangles.begin(), triangles.end(), t), triangles.end());
    }
    bool test() const
    {
        std::vector<Point> points = getPoints();
        const double eps = 1e-9; // tolerance
        for (const auto &t : triangles)
        {
            Point center = t.circumcircleCenter();
            double r2 = t.radiusOfCircle();
            if (r2 < 0)
                continue; // skip collinear (practically won't happen)

            for (auto &pi : points)
            {
                if (pi == t.a || pi == t.b || pi == t.c)
                    continue;
                double d2 = center.dist2(pi);
                if (std::abs(d2 - r2) < eps)
                {
                    std::cout << "Triangle " << t
                              << " has extra point " << pi
                              << " on its circumcircle.\n";
                    return false;
                }
            }
        }
        return true;
    }
    void InsertPoint(const Point &point)
    {
        std::cout << "Inserting point ";
        point.print();

        if (triangles.empty()) throw std::runtime_error("No triangles in triangulation");

        Triangle t = FindTriangleByPoint(point, triangles[0]); // triangle we are inserting into

        if (t.IsPointOnEdge(point))
        {
            std::cout << "point is on the edge" << std::endl;
            Edge CommonEdgeOfTwoTriangles = t.GetEdgeByPoint(point);

            Triangle oppositeT = FindAjasentTriangleByEdge(CommonEdgeOfTwoTriangles, t);

            Point oppositeVertex = oppositeT.oppositeVertex(CommonEdgeOfTwoTriangles);

            DeleteTriangle(t);
            DeleteTriangle(oppositeT);

            addTriangle(Triangle(t.a, t.b, point));
            addTriangle(Triangle(oppositeVertex, t.b, point));
            addTriangle(Triangle(oppositeVertex, t.c, point));
            addTriangle(Triangle(t.c, t.a, point));

            LegalizeEdge(point, Edge(t.a, t.b));
            LegalizeEdge(point, Edge(t.b, oppositeVertex));
            LegalizeEdge(point, Edge(oppositeVertex, t.c));
            LegalizeEdge(point, Edge(t.c, t.a));
            return;
        }
        std::cout << "point is not on edge" << std::endl;

        DeleteTriangle(t);

        addTriangle(Triangle(t.a, t.b, point));
        addTriangle(Triangle(t.b, t.c, point));
        addTriangle(Triangle(t.c, t.a, point));

        LegalizeEdge(point, Edge(t.a, t.b));
        LegalizeEdge(point, Edge(t.b, t.c));
        LegalizeEdge(point, Edge(t.c, t.a));
        test();
    }

    std::vector<Point> getBoundaryVertices(const Point &start, const Point &end) const
    {
        // Count edges
        std::map<Edge, int> edgeCount;
        for (const auto &tri : triangles)
        {
            for (const auto &e : tri.edges())
            {
                edgeCount[e]++;
            }
        }

        std::map<Point, std::vector<Point>> adj;
        for (const auto &[edge, count] : edgeCount)
        {
            if (count == 1) // boundary edge
            {
                adj[edge.start].push_back(edge.end);
                adj[edge.end].push_back(edge.start);
            }
        }

        std::vector<Point> path;
        Point cur = start;
        Point prev;
        while (!(cur == end))
        {
            path.push_back(cur);
            bool moved = false;
            auto it = adj.find(cur);
            if (it == adj.end()) break;
            for (const auto &nxt : it->second)
            {
                if (!(nxt == prev))
                {
                    prev = cur;
                    cur = nxt;
                    moved = true;
                    break;
                }
            }
            if (!moved) break; // cannot proceed
        }
        path.push_back(end);
        return path;
    }

    std::vector<Triangle> enforceQuality()
    {
        int i=0;
        if (max_area <= 0 && min_angle_deg <= 0)
            return triangles;

        bool improved = true;
        int iteration = 0, max_iterations = 10000;

        while (improved && iteration++ < max_iterations)
        {
            improved = false;

            // iterate by index, copy triangle by value to avoid invalidation issues when InsertPoint modifies triangles vector
            for (size_t idx = 0; idx < triangles.size(); ++idx)
            {
                Triangle tri = triangles[idx]; // copy

                double area = tri.area();
                double minAngle = tri.minAngleDeg();

                if ((max_area > 0 && area > max_area) ||
                    (min_angle_deg > 0 && minAngle < min_angle_deg))
                {
                    Point newPoint;

                    if (max_area > 0 && area > max_area)
                    {
                        // Centroid refinement for area
                        newPoint = (tri.a + tri.b + tri.c) / 3.0;
                    }
                    else
                    {
                        // Min angle violation → split longest edge
                        Edge e = tri.longestEdge();
                        newPoint = (e.start + e.end) / 2.0;
                    }

                    InsertPoint(newPoint); // retriangulate with new point
                    i++;
                    improved = true;
                    break; // restart scanning, since triangulation changed
                }
            }
        }
        std::cout<<"insertedpoints "<<i<<"\n";

        if (iteration >= max_iterations)
        {
            std::cerr << "Warning: enforceQuality reached iteration cap\n";
        }

        return triangles;
    }

    std::vector<Triangle> FindAjasentTrianglesByEdge(const Edge &edge) const
    {
        std::vector<Triangle> result;
        for (const auto &triangle : triangles)
        {
            if (triangle.containsEdge(edge))
                result.push_back(triangle);
        }
        return result;
    }
    bool OnExteriourFace(const Edge &edge) const
    {
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

    Triangle FindAjasentTriangleByEdge(const Edge &edge, const Triangle &currentTriangle) const
    {
        for (const auto &triangle : triangles)
        {
            if (triangle == currentTriangle)
                continue;
            if (triangle.containsEdge(edge))
                return triangle;
        }
        // print triangles for debugging -> throw
        throw std::runtime_error("No adjacent triangle found for the given edge.");
    }

    Triangle FindTriangleByPoint(const Point &p, const Triangle & /*startingTriangle*/) const
    {
        for (const auto &triangle : triangles)
        {
            if (triangle.IsPointInside(p) || triangle.IsPointOnEdge(p))
                return triangle;
        }
        throw std::runtime_error("No triangle contains the given point");
    }

    void LegalizeEdge(const Point &insertedPoint, const Edge &startEdge)
    {
        // iterative, stack-based legalization
        std::stack<Edge> st;
        st.push(startEdge);

        while (!st.empty())
        {
            Edge edge = st.top();
            st.pop();

            auto adj = FindAjasentTrianglesByEdge(edge);
            if (adj.size() != 2)
            {
                // boundary edge; nothing to do
                continue;
            }

            // identify which triangle contains the inserted point
            Triangle triWithPoint = adj[0].containsVertex(insertedPoint) ? adj[0] : adj[1];
            Triangle otherTri = adj[0].containsVertex(insertedPoint) ? adj[1] : adj[0];

            // opposite vertex in the neighboring triangle
            Point opp = otherTri.oppositeVertex(edge);

            // test: is opp inside circumcircle of triWithPoint?
            if (triWithPoint.IsPointInTriangleCircle(opp))
            {
                // flip
                Point opp1 = triWithPoint.oppositeVertex(edge); // opposite in triWithPoint
                Point opp2 = otherTri.oppositeVertex(edge);     // opposite in otherTri

                // remove old triangles
                DeleteTriangle(triWithPoint);
                DeleteTriangle(otherTri);

                // create flipped triangles (edge becomes opp1-opp2)
                Triangle newT1(opp1, opp2, edge.start);
                Triangle newT2(opp1, opp2, edge.end);
                addTriangle(newT1);
                addTriangle(newT2);

                // After flip, edges that might become illegal are:
                st.push(Edge(opp1, edge.start));
                st.push(Edge(opp1, edge.end));
                st.push(Edge(opp2, edge.start));
                st.push(Edge(opp2, edge.end));
            }
        }
    }

    void RemoveSuperTriangleTriangles(const Triangle &superTriangle)
    {
        // remove any triangle which contains any vertex of superTriangle
        // careful: use an index-based erase to avoid invalidating iteration
        std::vector<Triangle> kept;
        for (const auto &it : triangles)
        {
            if (it.containsVertex(superTriangle.a) ||
                it.containsVertex(superTriangle.b) ||
                it.containsVertex(superTriangle.c))
            {
                // skip (remove)
            }
            else
            {
                kept.push_back(it);
            }
        }
        triangles.swap(kept);
    }
};

Triangle createSuperTriangle(const std::vector<Point> &points)
{
    if (points.empty()) throw std::runtime_error("createSuperTriangle: empty points");

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
    if (deltaMax < EPSILON) deltaMax = 1.0; // avoid degenerate
    double midX = (minX + maxX) / 2;
    double midY = (minY + maxY) / 2;

    Point p1(midX - 20 * deltaMax, midY - deltaMax);
    Point p2(midX, midY + 20 * deltaMax);
    Point p3(midX + 20 * deltaMax, midY - deltaMax);

    return Triangle(p1, p2, p3);
}

// Example usage
int main()
{
    try
    {
        Point p1 = Point(1, 1);
        Point p2 = Point(3, 3);
        Point p3 = Point(2, 4);

        Triangle superT = createSuperTriangle({p1, p2, p3});
        superT.print();

        Triangulation triangulation;
        triangulation.addTriangle(superT);
        triangulation.InsertPoint(p1);
        triangulation.print();
        triangulation.InsertPoint(p2);
        triangulation.print();
        triangulation.InsertPoint(p3);
        triangulation.print();

        // set a small max area to force splitting if you want to test
        triangulation.setMaxArea(0.5);
        triangulation.enforceQuality();
        triangulation.print();
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
