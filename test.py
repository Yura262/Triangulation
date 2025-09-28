import numpy as np

def lift(p):
    x,y = p
    return np.array([x, y, x*x + y*y], dtype=float)

def plane_z_at(x, y, A, B, C):
    # A,B,C are 3D lifted points (numpy arrays)
    N = np.cross(B - A, C - A)   # normal vector
    print(N)
    Nx, Ny, Nz = N
    Ax, Ay, Az = A
    if abs(Nz) < 1e-12:
        raise ValueError("Nz nearly zero: circle may be undefined (colinear points).")
    return Az - (Nx*(x - Ax) + Ny*(y - Ay)) / Nz

# Example:
A = lift((1,1))
B = lift((3,3))
C = lift((2,4))

for p in [(8,2), (1,2)]:
    x,y = p
    parab = x*x + y*y
    zplane = plane_z_at(x, y, A, B, C)
    
    inside = parab < zplane
    print(p, "parab=", parab, "plane=", zplane, "inside?", inside)
