import numpy as np
import glob
import os
import re
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg



class FEMApp:
    def __init__(self, root):
        self.root = root
        self.root.title("FEM")
        self.root.geometry("1250x850")


        self.points = None
        self.triangles = None
        self.filename = "Не вибрано"
        

        self.K = None
        self.M = None
        self.F = None
        self.global_matrix = None 
        self.global_vector = None 


        self.dirichlet_nodes = []
        # self.neumann_nodes = []

        self.setup_ui()
        self.load_data_auto()

    def setup_ui(self):
        # --- Ліва панель керування ---
        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.pack(side=tk.LEFT, fill=tk.Y)

        # === Блок Файлу ===
        file_frame = ttk.LabelFrame(control_frame, text="Файл даних")
        file_frame.pack(fill=tk.X, pady=5)
        
        self.lbl_filename = ttk.Label(file_frame, text=f"Файл: {self.filename}", font=("Arial", 8), wraplength=250)
        self.lbl_filename.pack(pady=2, fill=tk.X)
        
        btn_file_box = ttk.Frame(file_frame)
        btn_file_box.pack(fill=tk.X, pady=2)
        
        ttk.Button(btn_file_box, text="Обрати файл...", command=self.browse_file).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=2)
        ttk.Button(btn_file_box, text="Авто (Останній)", command=self.load_data_auto).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=2)

        # === Параметри ===
        phys_frame = ttk.LabelFrame(control_frame, text="Параметри")
        phys_frame.pack(fill=tk.X, pady=5)

        ttk.Label(phys_frame, text="a1:").grid(row=0, column=0, padx=5, pady=2)
        self.entry_a1 = ttk.Entry(phys_frame, width=8); self.entry_a1.insert(0, "5.0"); self.entry_a1.grid(row=0, column=1)

        ttk.Label(phys_frame, text="a2:").grid(row=1, column=0, padx=5, pady=2)
        self.entry_a2 = ttk.Entry(phys_frame, width=8); self.entry_a2.insert(0, "8.0"); self.entry_a2.grid(row=1, column=1)

        ttk.Label(phys_frame, text="f_const:").grid(row=2, column=0, padx=5, pady=2)
        self.entry_f = ttk.Entry(phys_frame, width=8); self.entry_f.insert(0, "2.0"); self.entry_f.grid(row=2, column=1)

        ttk.Button(phys_frame, text="Застосувати", command=self.rebuild_system).grid(row=3, column=0, columnspan=2, pady=5, sticky="ew")

        # === Геометрія ===
        geo_frame = ttk.LabelFrame(control_frame, text="Вибір границі")
        geo_frame.pack(fill=tk.X, pady=5)

        # Поля вводу
        coord_grid = ttk.Frame(geo_frame)
        coord_grid.pack(fill=tk.X, padx=2)
        
        ttk.Label(coord_grid, text="P1 (x,y):").grid(row=0, column=0)
        self.entry_x1 = ttk.Entry(coord_grid, width=5); self.entry_x1.insert(0, "0.0"); self.entry_x1.grid(row=0, column=1)
        self.entry_y1 = ttk.Entry(coord_grid, width=5); self.entry_y1.insert(0, "0.0"); self.entry_y1.grid(row=0, column=2)

        ttk.Label(coord_grid, text="P2 (x,y):").grid(row=1, column=0)
        self.entry_x2 = ttk.Entry(coord_grid, width=5); self.entry_x2.insert(0, "1.0"); self.entry_x2.grid(row=1, column=1)
        self.entry_y2 = ttk.Entry(coord_grid, width=5); self.entry_y2.insert(0, "0.0"); self.entry_y2.grid(row=1, column=2)

        # кнопки вибору заготовок
        presets_frame = ttk.Frame(geo_frame)
        presets_frame.pack(fill=tk.X, pady=5)
        
        # P1=(0,0), P2=(1,0)
        ttk.Button(presets_frame, text="Низ (0,0->1,0)", 
                   command=lambda: self.set_geo_coords(0,0,1,0)).grid(row=0, column=0, sticky="ew", padx=1, pady=1)
        
        # P1=(0,0), P2=(0,2)
        ttk.Button(presets_frame, text="Ліво (0,0->0,2)", 
                   command=lambda: self.set_geo_coords(0,0,0,2)).grid(row=0, column=1, sticky="ew", padx=1, pady=1)
        
        # P1=(0,2), P2=(1,2)
        ttk.Button(presets_frame, text="Верх (0,2->1,2)", 
                   command=lambda: self.set_geo_coords(0,2,1,2)).grid(row=1, column=0, sticky="ew", padx=1, pady=1)
        
        # P1=(1,0), P2=(1,2)
        ttk.Button(presets_frame, text="Право (1,0->1,2)", 
                   command=lambda: self.set_geo_coords(1,0,1,2)).grid(row=1, column=1, sticky="ew", padx=1, pady=1)

        # === Граничні умови ===
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
        
        # bc_frame = ttk.LabelFrame(control_frame, text="Граничні умови")
        # bc_frame.pack(fill=tk.X, pady=5)

        # self.bc_type = tk.StringVar(value="Dirichlet")
        # combo = ttk.Combobox(bc_frame, textvariable=self.bc_type, state="readonly", 
        #                      values=("Dirichlet (Фіксоване u)", "Neumann (Потік/Сила)"))
        # combo.pack(fill=tk.X, padx=5, pady=5)
        # combo.bind("<<ComboboxSelected>>", self.update_bc_label)

        # val_frame = ttk.Frame(bc_frame)
        # val_frame.pack(fill=tk.X, pady=2)
        # self.val_label = ttk.Label(val_frame, text="Значення u =")
        # self.val_label.pack(side=tk.LEFT, padx=5)
        # self.entry_val = ttk.Entry(val_frame, width=10)
        # self.entry_val.insert(0, "0.0")
        # self.entry_val.pack(side=tk.LEFT, padx=5)

        # btn_frame = ttk.Frame(bc_frame)
        # btn_frame.pack(pady=5)
        # ttk.Button(btn_frame, text="Застосувати умову", command=self.apply_boundary_condition).pack(side=tk.LEFT, padx=2)
        # ttk.Button(btn_frame, text="Скинути умови", command=self.reset_bcs_only).pack(side=tk.LEFT, padx=2)

        # === Обчислення ===
        ttk.Separator(control_frame).pack(fill='x', pady=10)
        ttk.Button(control_frame, text="ОБЧИСЛИТИ ", command=self.solve_and_visualize).pack(fill=tk.X, pady=5)

        # Лог
        self.log_text = tk.Text(control_frame, height=12, width=40, font=("Consolas", 8))
        self.log_text.pack(fill=tk.BOTH, expand=True, pady=5)

        # --- Права панель (Plot) ---
        plot_frame = ttk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def set_geo_coords(self, x1, y1, x2, y2):
        """Допоміжна функція для кнопок швидкого вводу"""
        self.entry_x1.delete(0, tk.END); self.entry_x1.insert(0, str(x1))
        self.entry_y1.delete(0, tk.END); self.entry_y1.insert(0, str(y1))
        self.entry_x2.delete(0, tk.END); self.entry_x2.insert(0, str(x2))
        self.entry_y2.delete(0, tk.END); self.entry_y2.insert(0, str(y2))
        #self.log(f"Встановлено координати: ({x1},{y1}) -> ({x2},{y2})")
        
    def reset_data(self):
        """Скидає матриці до початкового стану"""
        self.global_matrix = self.original_matrix.copy()
        self.global_vector = self.original_vector.copy()
        self.highlighted_nodes = []
        self.plot_mesh()
        self.log("\nматриці скинуто.")
        
    def log(self, msg):
        self.log_text.insert(tk.END, msg + "\n")
        self.log_text.see(tk.END)

    def update_bc_label(self, event):
        if "Dirichlet" in self.bc_type.get():
            self.val_label.config(text="Значення u =")
        else:
            self.val_label.config(text="Потік/Сила F +=")

    # --- ФУНКЦІЇ ЗАВАНТАЖЕННЯ ДАНИХ ---
    def load_data_auto(self):
        """Шукає найновіший файл автоматично"""
        files = glob.glob('triangulation_export_*.txt')
        if not files:
            self.log("Авто: файли triangulation_export_*.txt не знайдено.")
            return
        latest = max(files, key=os.path.getmtime)
        self.parse_file(latest)

    def browse_file(self):
        """Відкриває діалог вибору файлу"""
        fname = filedialog.askopenfilename(filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if fname:
            self.parse_file(fname)

    def parse_file(self, filepath):
        """Парсить вибраний файл"""
        pts = []
        tris = []
        pre = re.compile(r'^point\s*\d*\s*\(\s*([\-+0-9.eE]+)\s*,\s*([\-+0-9.eE]+)\s*\)')
        tre = re.compile(r'^triangle\s*\d*\s*\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)')
        
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    s = line.strip()
                    if not s or s.startswith('#'):
                        continue
                    m = pre.match(s)
                    if m:
                        pts.append((float(m.group(1)), float(m.group(2))))
                        continue
                    m2 = tre.match(s)
                    if m2:
                        tris.append((int(m2.group(1)), int(m2.group(2)), int(m2.group(3))))
                        continue
            
            if not pts or not tris:
                raise ValueError("Файл порожній або має неправильний формат.")

            self.points = np.array(pts, dtype=float)
            self.triangles = np.array(tris, dtype=int)
            self.filename = os.path.basename(filepath)
            self.lbl_filename.config(text=f"Файл: {self.filename}")
            
            self.log(f"Завантажено: {self.filename} ({len(pts)} точок, {len(tris)} трикутників)")
            
            
            
            
            self.rebuild_system()

        except Exception as e:
            messagebox.showerror("Помилка читання", f"Не вдалося прочитати файл:\n{e}")




    def rebuild_system(self):
        self.build_matrices()
        self.plot_mesh()

    def reset_bcs_only(self):
        if self.K is None: return
        self.global_matrix = self.K + self.M
        self.global_vector = self.F.copy()
        self.dirichlet_nodes = []
        self.neumann_nodes = []
        self.plot_mesh()
        self.log("Граничні умови скинуто.")

    def build_matrices(self):
        if self.points is None: return
        
        try:
            a1 = float(self.entry_a1.get())
            a2 = float(self.entry_a2.get())
            f_const = float(self.entry_f.get())
        except ValueError:
            messagebox.showerror("Помилка", "Перевірте a1, a2, f_const")
            return

        N = len(self.points)
        K = np.zeros((N, N))
        M = np.zeros((N, N))
        F = np.zeros(N)

        for e, (g1, g2, g3) in enumerate(self.triangles):
            x1, y1 = self.points[g1]
            x2, y2 = self.points[g2]
            x3, y3 = self.points[g3]

            det = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
            area = 0.5 * abs(det)
            if area <= 1e-12: continue

            b = [y2 - y3, y3 - y1, y1 - y2]
            c = [x3 - x2, x1 - x3, x2 - x1]

            Kloc = np.zeros((3, 3))
            Mloc = np.zeros((3, 3))
            
            for a in range(3):
                for b_ in range(3):
                    Kloc[a, b_] = (1 / (4 * area)) * (a1 * b[a] * b[b_] + a2 * c[a] * c[b_])
                    Mloc[a, b_] = (area / 12.0) * (2 if a == b_ else 1)

            Floc = np.array([f_const * area / 3.0] * 3)

            global_indices = [g1, g2, g3]
            for a in range(3):
                I = global_indices[a]
                F[I] += Floc[a]
                for b_ in range(3):
                    J = global_indices[b_]
                    K[I, J] += Kloc[a, b_]
                    M[I, J] += Mloc[a, b_]

        self.K = K
        self.M = M
        self.F = F
        self.global_matrix = self.K + self.M
        self.global_vector = self.F.copy()
        self.dirichlet_nodes = []
        self.neumann_nodes = []

    def is_point_on_segment(self, point, p1, p2, tolerance=1e-5):
        px, py = point
        x1, y1 = p1
        x2, y2 = p2
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
        
        
        u_val = float(self.entry_val.get())
        self.dirichlet_nodes = []
        for idx, pt in enumerate(self.points):
            if self.is_point_on_segment(pt, (x1, y1), (x2, y2)):
                indices.append(idx)
                self.dirichlet_nodes.append(pt)

        if indices:
            self.log(f"ГУ застосовано до {len(indices)} вузлів.")
            for i in indices:
                self.global_matrix[i, :] = 0.0
                self.global_matrix[i, i] = 1.0
                self.global_vector[i] = u_val
            self.plot_mesh((x1, y1), (x2, y2))
        else:
            self.log("Вузлів на лінії не знайдено.")
            
            
        # bc_mode = self.bc_type.get()
        # for idx, pt in enumerate(self.points):
        #     if self.is_point_on_segment(pt, (x1, y1), (x2, y2)):
        #         indices.append(idx)
        #         if "Dirichlet" in bc_mode:
        #             self.dirichlet_nodes.append(pt)
        #             self.neumann_nodes = [p for p in self.neumann_nodes if not np.array_equal(p, pt)]
        #         else:
        #             self.neumann_nodes.append(pt)
        #             self.dirichlet_nodes = [p for p in self.dirichlet_nodes if not np.array_equal(p, pt)]

        # if not indices:
        #     self.log("Вузлів на лінії не знайдено.")
        #     return

        # self.log(f"Застосовано {bc_mode.split()[0]} до {len(indices)} вузлів.")
        
        # for i in indices:
        #     if "Dirichlet" in bc_mode:
        #         self.global_matrix[i, :] = 0.0
        #         self.global_matrix[i, i] = 1.0
        #         self.global_vector[i] = val
        #     else:
        #         self.global_vector[i] += val

        # self.plot_mesh((x1, y1), (x2, y2))

    def plot_mesh(self, segment_start=None, segment_end=None):
        self.ax.clear()
        
        self.ax.triplot(self.points[:, 0], self.points[:, 1], self.triangles, 'g-', lw=0.5, alpha=0.5)
        self.ax.plot(self.points[:, 0], self.points[:, 1], 'bo', markersize=3, label='Вузли', alpha=0.6)

        if self.dirichlet_nodes:
            d_arr = np.array(self.dirichlet_nodes)
            self.ax.plot(d_arr[:, 0], d_arr[:, 1], 'ro', markersize=6, label='Dirichlet')
        
        if self.neumann_nodes:
            n_arr = np.array(self.neumann_nodes)
            self.ax.plot(n_arr[:, 0], n_arr[:, 1], 'mo', markersize=6, label='Neumann')

        if segment_start and segment_end:
            self.ax.plot([segment_start[0], segment_end[0]], 
                         [segment_start[1], segment_end[1]], 
                         'k--', lw=2, label='Selection')

        for i, (px, py) in enumerate(self.points):
            self.ax.text(px, py, str(i), fontsize=8, color='black')

        self.ax.set_title("")
        self.ax.legend(loc='upper right')
        self.ax.axis('equal')
        self.canvas.draw()

    def solve_and_visualize(self):
        self.log("\nSolving System")
        try:
            u_solution = np.linalg.solve(self.global_matrix, self.global_vector)
            print("Розв'язок u:", u_solution)
            res_win = tk.Toplevel(self.root)
            res_win.title("Результати")
            res_win.geometry("1200x600")

            fig = plt.figure(figsize=(12, 6))

            # 2d heatmap
            ax2d = fig.add_subplot(1, 2, 1)
            contour = ax2d.tricontourf(self.points[:, 0], self.points[:, 1], self.triangles, u_solution, levels=100, cmap="viridis")
            ax2d.triplot(self.points[:, 0], self.points[:, 1], self.triangles, 'w-', lw=0.2, alpha=0.5)
            plt.colorbar(contour, ax=ax2d, label="u(x,y)")
            ax2d.set_title("2D Heatmap")
            ax2d.axis('equal')

            # 3D surface
            ax3d = fig.add_subplot(1, 2, 2, projection='3d')
            surf = ax3d.plot_trisurf(self.points[:, 0], self.points[:, 1], u_solution, 
                                     triangles=self.triangles, cmap='viridis', 
                                     edgecolor='none', linewidth=0, antialiased=True)
            plt.colorbar(surf, ax=ax3d, shrink=0.7)
            ax3d.set_title("3D Surface")
            
            
            x_range = np.ptp(self.points[:, 0]) # max - min
            y_range = np.ptp(self.points[:, 1])
            z_range = max(x_range, y_range) 
            
            ax3d.set_box_aspect((x_range, y_range, z_range))

            ax3d.set_xlabel("X")
            ax3d.set_ylabel("Y")
            ax3d.set_zlabel("u")

            canvas_res = FigureCanvasTkAgg(fig, master=res_win)
            canvas_res.draw()
            canvas_res.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            

        except np.linalg.LinAlgError as e:
            messagebox.showerror("Помилка", f"помилка")

if __name__ == "__main__":
    root = tk.Tk()
    app = FEMApp(root)
    root.mainloop()