import numpy as np
import glob
import os
import re

# ===== INPUT =====
# Try to read the latest export file produced by visualize.py
# Format expected (examples):

def _read_latest_export():
    files = glob.glob('triangulation_export_*.txt')
    if not files:
        return None, None, None
    latest = max(files, key=os.path.getmtime)
    pts = []
    tris = []
    # Regexes accept floats and optional spaces
    pre = re.compile(r'^point\s*\d*\s*\(\s*([\-+0-9.eE]+)\s*,\s*([\-+0-9.eE]+)\s*\)')
    tre = re.compile(r'^triangle\s*\d*\s*\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)')
    try:
        with open(latest, 'r') as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith('#'):
                    continue
                m = pre.match(s)
                if m:
                    x = float(m.group(1))
                    y = float(m.group(2))
                    pts.append((x, y))
                    continue
                m2 = tre.match(s)
                if m2:
                    p1 = int(m2.group(1))
                    p2 = int(m2.group(2))
                    p3 = int(m2.group(3))
                    tris.append((p1, p2, p3))
                    continue
    except Exception as e:
        print(f"Failed to read export file '{latest}': {e}")
        return None, None, None

    return pts if pts else None, tris if tris else None, latest


# Attempt to load latest export;
pts, tris, fname = _read_latest_export()
if pts is not None and tris is not None:
    points = pts
    triangles = tris
    print(f"Loaded points/triangles from: {fname}")
else:
    input(f"File import failed ({fname}).")
    exit()

# ===== CONSTANTS =====
a1 = 5     # coefficient for x-derivative
a2 = 8     # coefficient for y-derivative
f_const = 2
N = len(points)

# Initialize global matrices
K = np.zeros((N, N))
M = np.zeros((N, N))
F = np.zeros(N)

# ===== LOOP OVER ELEMENTS =====
for e, (g1, g2, g3) in enumerate(triangles):
    x1, y1 = points[g1]
    x2, y2 = points[g2]
    x3, y3 = points[g3]

    # --- 1. Compute area using determinant ---
    det = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
    area = 0.5 * abs(det)
    if area <= 0:
        raise ValueError(f"Element {e}: Invalid orientation or zero area (A={area})")

    # --- 2. Compute b_i, c_i ---
    b1 = y2 - y3
    b2 = y3 - y1
    b3 = y1 - y2

    c1 = x3 - x2
    c2 = x1 - x3
    c3 = x2 - x1

    b = [b1, b2, b3]
    c = [c1, c2, c3]

    # --- 3. Local matrices ---
    Kloc = np.zeros((3, 3))
    Mloc = np.zeros((3, 3))
    for a in range(3):
        for b_ in range(3):
            Kloc[a, b_] = (1 / (4 * area)) * (a1 * b[a] * b[b_] + a2 * c[a] * c[b_])
            Mloc[a, b_] = (area / 12.0) * (2 if a == b_ else 1)
    # print(f"Element {e}: Kloc =\n{Kloc}\nMloc =\n{Mloc}\n\n")
    with open("local_stiffness_matrices_Kloc.txt", "a") as f:
        f.write(f"Element {e}:\n{Kloc}\n\n")
    # for row in Kloc:    
    #     print(sum(row))
    # --- 4. Local right-hand side vector ---
    Floc = np.array([2 * area / 3.0] * 3)

    # --- 5. Assembly into global matrices ---
    global_indices = [g1, g2, g3]
    for a in range(3):
        I = global_indices[a]
        F[I] += Floc[a]
        for b_ in range(3):
            J = global_indices[b_]
            K[I, J] += Kloc[a, b_]
            M[I, J] += Mloc[a, b_]

# ===== AFTER ASSEMBLY =====
A = K + M  # final system matrix

# ===== PRINT RESULTS =====
np.set_printoptions(precision=4, suppress=True)

print("=== Global Stiffness Matrix K ===")
print(K)
print("\n=== Global Mass Matrix M ===")
print(M)
print("\n=== Global Right-hand Side F ===")
print(F)
print("\n=== Final System Matrix A = K + M ===")
print(A)
input()

np.savetxt("global_stiffness_matrix_K.txt", K, fmt="%.6f")
np.savetxt("final_system_matrix_A.txt", A, fmt="%.6f")