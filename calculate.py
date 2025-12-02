import numpy as np
import glob
import os
import re
INCLUDENEYMAN=True
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


a1 = 5
a2 = 8
f_const = 2
N = len(points)


K = np.zeros((N, N))
M = np.zeros((N, N))
F = np.zeros(N)


for e, (g1, g2, g3) in enumerate(triangles):
    x1, y1 = points[g1]
    x2, y2 = points[g2]
    x3, y3 = points[g3]

    # area
    det = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
    area = 0.5 * abs(det)
    if area <= 0:
        raise ValueError(f"Element {e}: Invalid orientation or zero area (A={area})")

    # b_i, c_i
    b1 = y2 - y3
    b2 = y3 - y1
    b3 = y1 - y2

    c1 = x3 - x2
    c2 = x1 - x3
    c3 = x2 - x1

    b = [b1, b2, b3]
    c = [c1, c2, c3]

    # Local matrices
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

    Floc = np.array([f_const * area / 3.0] * 3)

    #global matrix
    global_indices = [g1, g2, g3]
    for a in range(3):
        I = global_indices[a]
        F[I] += Floc[a]
        for b_ in range(3):
            J = global_indices[b_]
            K[I, J] += Kloc[a, b_]
            M[I, J] += Mloc[a, b_]


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
np.savetxt("global_rhs_F.txt", F, fmt="%.6f")



import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d import Axes3D


class FEMBoundaryApp:
    def __init__(self, root):
        self.root = root
        self.root.title("FEM Boundary Conditions Visualizer")
        self.root.geometry("1100x700")

        # 1. Ініціалізація даних (Сітка та Матриці)
        self.init_fem_data()

        # 2. Налаштування графічного інтерфейсу
        self.setup_ui()

        # 3. Первинна візуалізація
        self.plot_mesh()

    def init_fem_data(self):
        """
        Генерує просту тестову сітку (прямокутник) та порожні глобальні матриці.
        """
        # Створення точок (сітка 4x4 для прикладу, від -2 до 2)
        x = np.linspace(-2, 2, 5)
        y = np.linspace(0, 3, 4)
        xx, yy = np.meshgrid(x, y)
        # self.points = pts#np.vstack([xx.ravel(), yy.ravel()]).T  # shape (N, 2)

        # Проста тріангуляція Делоне для створення елементів
        # from scipy.spatial import Delaunay
        # tri = Delaunay(self.points)
        # self.triangles = tris#tri.simplices
        self.points = np.array(pts, dtype=float)
        self.triangles = np.array(tris, dtype=int)
        
        
        self.num_nodes = len(self.points)
        # Ініціалізація глобальної матриці жорсткості (A або K) та вектора сили (F)
        # Заповнюємо випадковими/діагональними даними для наочності, ніби система зібрана
        self.global_matrix = A#np.eye(self.num_nodes) * 4.0 + np.random.rand(self.num_nodes, self.num_nodes) * 0.5
        # Робимо симетричною (як це часто буває в FEM)
        #self.global_matrix = (self.global_matrix + self.global_matrix.T) / 2
        
        self.global_vector = F# np.zeros(self.num_nodes)

        # Зберігаємо копію "чистих" матриць, щоб можна було скидати зміни
        self.original_matrix = self.global_matrix.copy()
        self.original_vector = self.global_vector.copy()

        # Список для зберігання виділених вузлів
        self.highlighted_nodes = []
        if INCLUDENEYMAN:
            self.dirichlet_nodes = [] # Список точок з умовою Діріхле
            self.neumann_nodes = []

    def setup_ui(self):
        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.pack(side=tk.LEFT, fill=tk.Y)

        ttk.Label(control_frame, text="Керування Границею", font=("Arial", 12, "bold")).pack(pady=10)

        # --- Координати ---
        input_frame = ttk.LabelFrame(control_frame, text="Геометрія границі")
        input_frame.pack(fill=tk.X, pady=5)

        ttk.Label(input_frame, text="Point 1 (x,y):").grid(row=0, column=0)
        self.entry_x1 = ttk.Entry(input_frame, width=5); self.entry_x1.insert(0, "0.0"); self.entry_x1.grid(row=0, column=1)
        self.entry_y1 = ttk.Entry(input_frame, width=5); self.entry_y1.insert(0, "0.0"); self.entry_y1.grid(row=0, column=2)

        ttk.Label(input_frame, text="Point 2 (x,y):").grid(row=1, column=0)
        self.entry_x2 = ttk.Entry(input_frame, width=5); self.entry_x2.insert(0, "1.0"); self.entry_x2.grid(row=1, column=1)
        self.entry_y2 = ttk.Entry(input_frame, width=5); self.entry_y2.insert(0, "0.0"); self.entry_y2.grid(row=1, column=2)
        
        btn_frame = ttk.Frame(control_frame)
        btn_frame.pack(pady=10)
        
        
        if INCLUDENEYMAN:
        # --- Вибір типу умови ---
            type_frame = ttk.LabelFrame(control_frame, text="Тип умови")
            type_frame.pack(fill=tk.X, pady=10)
        
            self.bc_type = tk.StringVar(value="Dirichlet")
            # Combobox для вибору
            combo = ttk.Combobox(type_frame, textvariable=self.bc_type, state="readonly", 
                                values=("Dirichlet (Фіксоване значення)", "Neumann (Потік/Сила)"))
            combo.pack(fill=tk.X, padx=5, pady=5)
            combo.bind("<<ComboboxSelected>>", self.update_label_text)

            self.val_label = ttk.Label(type_frame, text="Значення u =")
            self.val_label.pack(side=tk.LEFT, padx=5)
            
            self.entry_val = ttk.Entry(type_frame, width=10)
            self.entry_val.pack(side=tk.LEFT, padx=5)
            self.entry_val.insert(0, "10.0")
            
            ttk.Button(btn_frame, text="Застосувати", command=self.apply_boundary_condition).pack(side=tk.LEFT, padx=2)
            ttk.Button(btn_frame, text="Скинути Все", command=self.reset_data).pack(side=tk.LEFT, padx=2)
        else:
            # Значення ГУ
            val_frame = ttk.LabelFrame(control_frame, text="Умова (Діріхле)")
            val_frame.pack(fill=tk.X, pady=10)
            ttk.Label(val_frame, text="u =").pack(side=tk.LEFT)
            self.entry_val = ttk.Entry(val_frame, width=10)
            self.entry_val.pack(side=tk.LEFT, padx=5)
            self.entry_val.insert(0, "10.0")

            # Кнопки
            btn_frame = ttk.Frame(control_frame)
            btn_frame.pack(pady=10)
            ttk.Button(btn_frame, text="Застосувати", command=self.apply_boundary_condition).pack(side=tk.LEFT, padx=2)
            ttk.Button(btn_frame, text="Скинути", command=self.reset_data).pack(side=tk.LEFT, padx=2)
        # --- Кнопки ---
        
        

        ttk.Separator(control_frame).pack(fill='x', pady=10)
        ttk.Button(control_frame, text="ОБЧИСЛИТИ (2D + 3D)", command=self.solve_and_visualize).pack(fill=tk.X, pady=5)

        self.log_text = tk.Text(control_frame, height=15, width=35, font=("Consolas", 8))
        self.log_text.pack(fill=tk.BOTH, expand=True, pady=10)

        # Права панель
        plot_frame = ttk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def update_label_text(self, event):
        """Змінює підпис поля вводу залежно від вибору"""
        if "Dirichlet" in self.bc_type.get():
            self.val_label.config(text="Значення u =")
        else:
            self.val_label.config(text="Потік/Сила F +=")

    def is_point_on_segment(self, point, p1, p2, tolerance=1e-5):
        """
        Перевіряє, чи лежить точка 'point' на відрізку між p1 та p2.
        Використовує перевірку відстаней: dist(A,C) + dist(C,B) == dist(A,B)
        """
        px, py = point
        x1, y1 = p1
        x2, y2 = p2

        # Відстані
        d_p1_p = np.sqrt((px - x1)**2 + (py - y1)**2)
        d_p_p2 = np.sqrt((px - x2)**2 + (py - y2)**2)
        d_p1_p2 = np.sqrt((x1 - x2)**2 + (y1 - y2)**2)

        return np.isclose(d_p1_p + d_p_p2, d_p1_p2, atol=tolerance)

    def apply_boundary_condition(self):
        try:
            x1, y1 = float(self.entry_x1.get()), float(self.entry_y1.get())
            x2, y2 = float(self.entry_x2.get()), float(self.entry_y2.get())
            val = float(self.entry_val.get())
        except ValueError:
            return

        indices = []
        if INCLUDENEYMAN:
            bc_mode = self.bc_type.get()
            
            # Знаходимо вузли
            for idx, pt in enumerate(self.points):
                if self.is_point_on_segment(pt, (x1, y1), (x2, y2)):
                    indices.append(idx)
                    # Додаємо у відповідний список для візуалізації
                    if "Dirichlet" in bc_mode:
                        self.dirichlet_nodes.append(pt)
                        # Видаляємо з Неймана, якщо він там був, щоб уникнути конфлікту кольорів
                        self.neumann_nodes = [p for p in self.neumann_nodes if not np.array_equal(p, pt)]
                    else:
                        self.neumann_nodes.append(pt)
                        self.dirichlet_nodes = [p for p in self.dirichlet_nodes if not np.array_equal(p, pt)]

            if not indices:
                self.log("Вузлів на лінії не знайдено.")
                return

            self.log(f"\n--- Apply {bc_mode.split()[0]} ---")
            
            for i in indices:
                if "Dirichlet" in bc_mode:
                    # 1. Dirichlet: Модифікація матриці (Penalty method або Substitution)
                    # Тут використовуємо Substitution: рядок = 0, діагональ = 1, RHS = val
                    self.global_matrix[i, :] = 0.0
                    self.global_matrix[i, i] = 1.0
                    self.global_vector[i] = val
                    self.log(f"Node {i}: Fixed u = {val}")
                else:
                    # 2. Neumann: Модифікація ТІЛЬКИ вектора сили
                    # Ми додаємо значення до існуючої сили (накладання потоків)
                    self.global_vector[i] += val
                    self.log(f"Node {i}: Added load/flux += {val}")

            self.plot_mesh((x1, y1), (x2, y2))
        else:
            u_val = float(self.entry_val.get())
            self.highlighted_nodes = []
            for idx, pt in enumerate(self.points):
                if self.is_point_on_segment(pt, (x1, y1), (x2, y2)):
                    indices.append(idx)
                    self.highlighted_nodes.append(pt)

            if indices:
                self.log(f"ГУ застосовано до {len(indices)} вузлів.")
                for i in indices:
                    self.global_matrix[i, :] = 0.0
                    self.global_matrix[i, i] = 1.0
                    self.global_vector[i] = u_val
                self.plot_mesh((x1, y1), (x2, y2))
            else:
                self.log("Вузлів на лінії не знайдено.")
        

    def plot_mesh(self, segment_start=None, segment_end=None):
        self.ax.clear()

        # Малювання тріангуляції
        self.ax.triplot(self.points[:, 0], self.points[:, 1], self.triangles, 'g-', lw=0.5, alpha=0.5)
        
        # Малювання всіх вузлів синім
        self.ax.plot(self.points[:, 0], self.points[:, 1], 'bo', markersize=4, label='Вузли')

        # Якщо є виділені вузли - малюємо червоним
        if self.highlighted_nodes:
            hx = [p[0] for p in self.highlighted_nodes]
            hy = [p[1] for p in self.highlighted_nodes]
            self.ax.plot(hx, hy, 'ro', markersize=8, label='Граничні вузли')

        # Малювання лінії границі, яку ввів користувач
        if segment_start and segment_end:
            self.ax.plot([segment_start[0], segment_end[0]], 
                         [segment_start[1], segment_end[1]], 
                         'r--', lw=2, label='Лінія границі')

        # Підписи індексів вузлів (щоб користувач бачив i)
        for i, (px, py) in enumerate(self.points):
            self.ax.text(px, py + 0.1, str(i), fontsize=9, color='black', ha='center')

        self.ax.set_title("FEM Сітка та Граничні умови")
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.legend()
        self.ax.grid(True)
        self.ax.axis('equal')
        self.canvas.draw()

    def reset_data(self):
        """Скидає матриці до початкового стану"""
        self.global_matrix = self.original_matrix.copy()
        self.global_vector = self.original_vector.copy()
        self.highlighted_nodes = []
        self.plot_mesh()
        self.log("\n--- МАТРИЦІ СКИНУТО ---")

    def log(self, message):
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
    def solve_and_visualize(self):
        self.log("\n--- Start Solver ---")
        try:
            u_solution = np.linalg.solve(self.global_matrix, self.global_vector)
            print("Розв'язок u:", u_solution)
            # Створюємо вікно результатів
            result_window = tk.Toplevel(self.root)
            result_window.title("Результати: 2D Heatmap та 3D Поверхня")
            result_window.geometry("1200x600") # Робимо ширшим для двох графіків

            # Створюємо Figure з двома підграфіками (1 рядок, 2 колонки)
            # subplot_kw={'projection': '3d'} потрібен для другого графіка
            fig = plt.figure(figsize=(12, 6))

            # --- 1. Лівий графік: 2D Heatmap ---
            ax2d = fig.add_subplot(1, 2, 1)
            contour = ax2d.tricontourf(
                self.points[:, 0], self.points[:, 1], self.triangles, 
                u_solution, levels=500, cmap="viridis"
            )
            ax2d.triplot(self.points[:, 0], self.points[:, 1], self.triangles, 'w-', lw=0.5, alpha=0.3)
            plt.colorbar(contour, ax=ax2d, shrink=0.7, label="u(x,y)")
            ax2d.set_title("2D Результат (Heatmap)")
            ax2d.set_xlabel("X"); ax2d.set_ylabel("Y")
            ax2d.axis('equal')

            # --- 2. Правий графік: 3D Surface ---
            ax3d = fig.add_subplot(1, 2, 2, projection='3d')
            
            # Використовуємо plot_trisurf
            surf = ax3d.plot_trisurf(
                self.points[:, 0], 
                self.points[:, 1], 
                u_solution, 
                triangles=self.triangles, 
                cmap='viridis', 
                edgecolor='none', 
                linewidth=0.2,
                antialiased=True
            )
            
            # Додаємо colorbar і для 3D (опціонально)
            plt.colorbar(surf, ax=ax3d, shrink=0.7, label="u(x,y)")
            
            ax3d.set_title("3D Результат (Деформація/Поле)")
            ax3d.set_xlabel("X")
            ax3d.set_ylabel("Y")
            ax3d.set_zlabel("Z (Value)")

            # Оптимізація вигляду
            plt.tight_layout()

            # Вбудовування в вікно
            canvas_res = FigureCanvasTkAgg(fig, master=result_window)
            canvas_res.draw()
            canvas_res.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        except np.linalg.LinAlgError as e:
            messagebox.showerror("Error", f"Solver failed: {e}")
            
        self.log("--- Solver Finished ---")
if __name__ == "__main__":
    root = tk.Tk()
    app = FEMBoundaryApp(root)
    root.mainloop()