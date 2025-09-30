#!/usr/bin/env python3
"""
Interactive Delaunay Triangulation Visualizer
Extended with:
 - enforce quality: angle & area inputs + button
 - boundary extraction: start/end inputs (x,y) + highlight button
 - manual point insertion: x & y inputs + button
 - a checkbox to toggle boundary highlight
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.widgets import Button, CheckButtons, TextBox
from matplotlib.patches import Rectangle

# Add your build directory to path
module_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'build/Debug'))
sys.path.insert(0, module_dir)

try:
    import triangulation
except ImportError:
    print("Could not import triangulation module. Make sure it's compiled.")
    print(f"Looking in: {module_dir}")
    sys.exit(1)


class DelaunayGUI:
    def __init__(self):
        self.enable_click_insert = False   # new flag

        # Initialize figure and axis
        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        self.fig.canvas.manager.set_window_title('Interactive Delaunay Triangulation')

        # Adjust layout to make room for controls
        plt.subplots_adjust(left=0.1, bottom=0.25, right=0.85, top=0.95)

        # Data storage
        self.points = []
        self.triangulation_obj = None
        self.triangles = []

        # Display options
        self.show_vertices = True
        self.show_vertex_numbers = True
        self.show_triangles = True
        self.show_triangle_numbers = True
        self.show_edges = True
        self.show_circumcircles = False
        self.show_boundary_highlight = False
        self._skip_recompute = False
        # Boundary path (list of (x,y)) when user requests highlight
        self.boundary_path = None
        self.boundary_artists = []

        # Plot elements (store references for updating)
        self.vertex_scatter = None
        self.vertex_labels = []
        self.triangle_patches = []
        self.triangle_labels = []
        self.edge_lines = []
        self.circumcircles = []

        # Setup the plot and UI
        self.setup_plot()
        self.setup_controls()

        # Connect mouse click event
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

        # Initial draw
        self.update_plot()

    def setup_plot(self):
        """Setup the initial plot appearance"""
        self.ax.set_xlabel('X', fontsize=12)
        self.ax.set_ylabel('Y', fontsize=12)
        self.ax.set_title('Click to add points', fontsize=14, fontweight='bold')
        self.ax.grid(True, alpha=0.3)
        self.ax.set_aspect('equal')

        # Set initial limits
        self.ax.set_xlim(-10, 10)
        self.ax.set_ylim(-10, 10)
    def setup_controls(self):
        """Setup UI control buttons, textboxes and checkboxes"""
          # Checkbox for display options
    
        # --- background panels (frames) ---
        # Insert point frame
        self.fig.patches.append(Rectangle((0.05, 0.13), 0.88, 0.06,
                                        transform=self.fig.transFigure,
                                        facecolor="lightgrey", alpha=0.3, zorder=-1))

        # Enforce quality frame
        self.fig.patches.append(Rectangle((0.65, 0.07), 0.28, 0.06,
                                        transform=self.fig.transFigure,
                                        facecolor="lightgrey", alpha=0.3, zorder=-1))

        # Boundary frame
        self.fig.patches.append(Rectangle((0.05, 0.07), 0.60, 0.06,
                                        transform=self.fig.transFigure,
                                        facecolor="lightgrey", alpha=0.3, zorder=-1))

        # --- buttons row ---
        ax_clear = plt.axes([0.10, 0.14, 0.10, 0.04])
        self.btn_clear = Button(ax_clear, 'Clear All', color='lightcoral', hovercolor='red')
        self.btn_clear.on_clicked(self.clear_all)

        ax_reset = plt.axes([0.21, 0.14, 0.10, 0.04])
        self.btn_reset = Button(ax_reset, 'Reset View', color='lightblue', hovercolor='blue')
        self.btn_reset.on_clicked(self.reset_view)

        ax_random = plt.axes([0.32, 0.14, 0.10, 0.04])
        self.btn_random = Button(ax_random, 'Add Random', color='lightgreen', hovercolor='green')
        self.btn_random.on_clicked(self.add_random_points)

        ax_update = plt.axes([0.43, 0.14, 0.10, 0.04])
        self.btn_update = Button(ax_update, 'Update Plot', color='lightblue', hovercolor='blue')
        self.btn_update.on_clicked(self.update_plot)

        # --- Insert point fields ---
        ax_ins_x = plt.axes([0.67, 0.14, 0.05, 0.04])
        self.tb_ins_x = TextBox(ax_ins_x, 'X ', initial='')
        ax_ins_y = plt.axes([0.73, 0.14, 0.05, 0.04])
        self.tb_ins_y = TextBox(ax_ins_y, 'Y ', initial='')
        ax_ins_btn = plt.axes([0.80, 0.14, 0.07, 0.04])
        self.btn_insert = Button(ax_ins_btn, 'Insert', color='lightgreen')
        self.btn_insert.on_clicked(self.insert_point_from_inputs)

        # --- Enforce quality (angle & area) ---
        ax_angle = plt.axes([0.69,0.08, 0.05, 0.04])
        self.tb_angle = TextBox(ax_angle, 'Angle°', initial='')
        ax_area = plt.axes([0.79, 0.08, 0.05, 0.04])
        self.tb_area = TextBox(ax_area, 'Max Area', initial='')
        ax_enforce = plt.axes([0.85, 0.08, 0.07, 0.04])
        self.btn_enforce = Button(ax_enforce, 'Enforce', color='lightcoral')
        self.btn_enforce.on_clicked(self.enforce_quality_from_inputs)

        # # --- Boundary input (start & end as "x,y") ---
        # ax_bstart = plt.axes([0.11, 0.08, 0.15, 0.04])
        # self.tb_bstart = TextBox(ax_bstart, 'Start (x,y)', initial='')
        # ax_bend = plt.axes([0.35, 0.08, 0.15, 0.04])
        # self.tb_bend = TextBox(ax_bend, 'End (x,y)', initial='')
        # ax_bbtn = plt.axes([0.53, 0.08, 0.11, 0.04])
        # self.btn_boundary = Button(ax_bbtn, 'Highlight', color='orange')
        # self.btn_boundary.on_clicked(self.highlight_boundary_from_inputs)

        # Checkbox for display options
        ax_check = plt.axes([0.87, 0.30, 0.12, 0.35])  # a bit taller
        labels = ['Vertices', 'Vertex Numbers', 'Triangles',
                'Triangle Numbers', 'Edges', 'Circumcircles',
                'Boundary Highlight', 'Click Insert']   # added here
        visibility = [self.show_vertices, self.show_vertex_numbers, self.show_triangles,
                    self.show_triangle_numbers, self.show_edges,
                    self.show_circumcircles, self.show_boundary_highlight,
                    self.enable_click_insert]            # added here
        self.check = CheckButtons(ax_check, labels, visibility)
        self.check.on_clicked(self.update_visibility)

        # Info text
        self.info_text = self.ax.text(0.02, 0.98, '', transform=self.ax.transAxes,
                                    fontsize=10, verticalalignment='top',
                                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    def on_click(self, event):
        """Handle mouse click events"""
        # Only process clicks within the main axes
        # if not self.enable_click_insert:   # 🔒 only allow if enabled
        #     return
        if event.inaxes != self.ax:
            return

        # Only process left clicks
        if event.button != 1:
            return

        # Add new point
        x, y = event.xdata, event.ydata
        self.add_point(x, y)

    def add_point(self, x, y):
        """Add a point and update triangulation"""
        # Check if point already exists
        for px, py in self.points:
            if px == x and py == y:
                return

        self.points.append((x, y))
        self.compute_triangulation()
        self.update_plot()

    def insert_point_from_inputs(self, event):
        """Insert point from X/Y textboxes"""
        try:
            sx = self.tb_ins_x.text.strip()
            sy = self.tb_ins_y.text.strip()
            if sx == '' or sy == '':
                return
            x = float(sx)
            y = float(sy)
            self.add_point(x, y)
            # clear boxes for convenience
            self.tb_ins_x.set_val('')
            self.tb_ins_y.set_val('')
        except Exception as e:
            print("Insert point parse error:", e)

    def enforce_quality_from_inputs(self, event):
        """Read min angle & max area; call C++ enforce routine and update points/triangles"""
        try:
            a_str = self.tb_angle.text.strip()
            ar_str = self.tb_area.text.strip()
            angle = float(a_str) if a_str != '' else 0.0
            area = float(ar_str) if ar_str != '' else 0.0

            # make sure triangulation exists
            self.compute_triangulation()
            if self.triangulation_obj is None:
                print("No triangulation available to enforce quality on")
                return

            # Set parameters and call enforceQuality on the existing triangulation object.
            try:
                self.triangulation_obj.setMinAngle(angle)
                self.triangulation_obj.setMaxArea(area)

                # The C++ method may modify the triangulation in-place and possibly return a list.
                # We call it and then re-read triangles/vertices from the triangulation object.
                _ret = self.triangulation_obj.enforceQuality()

                # Read triangles from the (now modified) triangulation object
                try:
                    cpp_triangles = self.triangulation_obj.getTriangles()
                except Exception as ee:
                    print("getTriangles failed after enforceQuality:", ee)
                    cpp_triangles = []

                # Build a new unique points list from vertices present in returned triangles.
                new_points = []
                index_map = {}  # (x,y) -> index in new_points
                new_triangles = []

                for tri in cpp_triangles:
                    # ensure a, b, c exist on tri (depending on your binding names)
                    verts = (tri.a, tri.b, tri.c)
                    tri_idxs = []
                    for v in verts:
                        key = (float(v.x), float(v.y))
                        if key not in index_map:
                            index_map[key] = len(new_points)
                            new_points.append(key)
                        tri_idxs.append(index_map[key])
                    if len(tri_idxs) == 3:
                        new_triangles.append(tri_idxs)

                # Replace current points/triangles with the enforced version
                if new_points:
                    self.points = new_points
                    self.triangles = new_triangles
                else:
                    # If no triangles returned, fall back to recomputing from points
                    print("Warning: enforceQuality produced no triangles; keeping current points/triangles")

                # Prevent immediate recompute in update_plot (compute_triangulation rebuilds from self.points
                # which may otherwise re-create the triangulation differently). We'll skip one recompute and
                # let plot use the triangles we just extracted.
                self._skip_recompute = True

            except Exception as e:
                print("enforceQuality call failed:", e)

            # update display (triangulation changed)
            self.update_plot()
        except Exception as e:
            print("Enforce quality parse error:", e)
            
    def highlight_boundary_from_inputs(self, event):
        """Parse start/end, compute boundary via C++ and draw it"""
        try:
            s = self.tb_bstart.text.strip()
            e = self.tb_bend.text.strip()
            if s == '' or e == '':
                print("Please enter both start and end points as x,y")
                return
            sx, sy = [float(v.strip()) for v in s.split(',')]
            ex, ey = [float(v.strip()) for v in e.split(',')]

            self.compute_triangulation()
            if self.triangulation_obj is None:
                print("No triangulation available for boundary extraction")
                return

            pstart = triangulation.Point(sx, sy)
            pend = triangulation.Point(ex, ey)

            try:
                py_pts = self.triangulation_obj.getBoundaryVertices(pstart, pend)
            except Exception as ee:
                print("getBoundaryVertices failed:", ee)
                return

            # convert to python coords
            path = [(pt.x, pt.y) for pt in py_pts]
            if not path:
                print("Empty boundary returned")
                return

            self.boundary_path = path
            # enable highlight checkbox state and drawing
            self.show_boundary_highlight = True
            self.update_plot()
        except Exception as e:
            print("Boundary parse/draw error:", e)

    def compute_triangulation(self):
        """Compute the Delaunay triangulation"""
        if len(self.points) < 3:
            self.triangles = []
            self.triangulation_obj = None
            return

        try:
            # Create points for C++
            cpp_points = [triangulation.Point(x, y) for x, y in self.points]

            # Create super triangle
            super_triangle = triangulation.createSuperTriangle(cpp_points)

            # Create triangulation and add super triangle
            self.triangulation_obj = triangulation.Triangulation()
            self.triangulation_obj.addTriangle(super_triangle)

            # Insert all points
            for point in cpp_points:
                self.triangulation_obj.InsertPoint(point)
            self.triangulation_obj.RemoveSuperTriangleTriangles(super_triangle)
            
            # Get triangles
            cpp_triangles = self.triangulation_obj.getTriangles()

            # Convert triangles to Python format as indices
            self.triangles = []
            for tri in cpp_triangles:
                indices = []
                for vertex in [tri.a, tri.b, tri.c]:
                    for i, (px, py) in enumerate(self.points):
                        if px == vertex.x and py == vertex.y:
                            indices.append(i)
                            break
                if len(indices) == 3:
                    self.triangles.append(indices)

        except Exception as e:
            print(f"Triangulation error: {e}")
            self.triangles = []
            self.triangulation_obj = None

    def update_plot(self, event=None):
        """Update the plot with current points and triangulation"""
        # If a previous operation (enforceQuality) already produced an authoritative
        # triangulation/points, skip one recomputation so we don't overwrite the C++ result.
        if self._skip_recompute:
            # consume the skip flag and do not call compute_triangulation() this time
            self._skip_recompute = False
        else:
            self.compute_triangulation()
        # Clear old plot elements
        self.clear_plot_elements()

        if not self.points:
            self.info_text.set_text('Click to add points\nPoints: 0\nTriangles: 0')
            self.fig.canvas.draw_idle()
            return

        # Convert points to numpy array
        points_array = np.array(self.points)

        # Draw triangles
        if self.show_triangles and self.triangles:
            for i, tri_indices in enumerate(self.triangles):
                if len(tri_indices) == 3:
                    triangle_points = points_array[tri_indices]

                    # Draw triangle fill
                    poly = Polygon(triangle_points, alpha=0.2, facecolor='lightblue',
                                   edgecolor='none')
                    self.ax.add_patch(poly)
                    self.triangle_patches.append(poly)

                    # Draw triangle number
                    if self.show_triangle_numbers:
                        centroid = triangle_points.mean(axis=0)
                        label = self.ax.text(centroid[0], centroid[1], str(i),
                                             ha='center', va='center',
                                             fontsize=9, color='blue',
                                             bbox=dict(boxstyle='round,pad=0.3',
                                                       facecolor='white', alpha=0.7))
                        self.triangle_labels.append(label)

        # Draw edges
        if self.show_edges and self.triangles:
            drawn_edges = set()
            for tri_indices in self.triangles:
                if len(tri_indices) == 3:
                    for j in range(3):
                        edge = tuple(sorted([tri_indices[j], tri_indices[(j + 1) % 3]]))
                        if edge not in drawn_edges:
                            drawn_edges.add(edge)
                            p1 = points_array[edge[0]]
                            p2 = points_array[edge[1]]
                            line, = self.ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                                                 'b-', linewidth=1, alpha=0.6)
                            self.edge_lines.append(line)

        # Draw circumcircles
        if self.show_circumcircles and self.triangles:
            for tri_indices in self.triangles:
                if len(tri_indices) == 3:
                    tri_points = points_array[tri_indices]
                    center, radius = self.compute_circumcircle(tri_points)
                    if center is not None and radius is not None:
                        circle = plt.Circle(center, radius, fill=False,
                                            edgecolor='green', alpha=0.3, linestyle='--')
                        self.ax.add_patch(circle)
                        self.circumcircles.append(circle)

        # Draw vertices
        if self.show_vertices:
            self.vertex_scatter = self.ax.scatter(points_array[:, 0], points_array[:, 1],
                                                 c='red', s=50, zorder=5, picker=True)

        # Draw vertex numbers
        if self.show_vertex_numbers:
            for i, (x, y) in enumerate(self.points):
                label = self.ax.text(x + 0.2, y + 0.2, str(i),
                                    fontsize=10, color='darkred', fontweight='bold')
                self.vertex_labels.append(label)

        # Draw boundary highlight if requested
        if self.show_boundary_highlight and self.boundary_path:
            xs = [p[0] for p in self.boundary_path]
            ys = [p[1] for p in self.boundary_path]
            # draw polyline
            line, = self.ax.plot(xs, ys, linestyle='-', linewidth=2, marker='o', markersize=6, color='orange', zorder=6)
            self.boundary_artists.append(line)
            # also highlight vertices with scatter
            scatter = self.ax.scatter(xs, ys, c='orange', s=80, edgecolors='k', zorder=7)
            self.boundary_artists.append(scatter)

        # Update info text
        self.info_text.set_text(f'Click to add points\nPoints: {len(self.points)}\n'
                                f'Triangles: {len(self.triangles)}')

        # Update axis limits if needed
        self.auto_scale()

        # Redraw
        self.fig.canvas.draw_idle()

    def clear_plot_elements(self):
        """Clear all plot elements"""
        # Remove vertex scatter
        if self.vertex_scatter:
            self.vertex_scatter.remove()
            self.vertex_scatter = None

        # Remove vertex labels
        for label in self.vertex_labels:
            label.remove()
        self.vertex_labels = []

        # Remove triangle patches
        for patch in self.triangle_patches:
            patch.remove()
        self.triangle_patches = []

        # Remove triangle labels
        for label in self.triangle_labels:
            label.remove()
        self.triangle_labels = []

        # Remove edge lines
        for line in self.edge_lines:
            line.remove()
        self.edge_lines = []

        # Remove circumcircles
        for circle in self.circumcircles:
            circle.remove()
        self.circumcircles = []

        # Remove boundary artists
        for a in self.boundary_artists:
            try:
                a.remove()
            except Exception:
                pass
        self.boundary_artists = []

    def clear_all(self, event):
        """Clear all points and triangulation"""
        self.points = []
        self.triangles = []
        self.triangulation_obj = None
        self.boundary_path = None
        self.update_plot()

    def reset_view(self, event):
        """Reset the view to default"""
        self.ax.set_xlim(-10, 10)
        self.ax.set_ylim(-10, 10)
        self.fig.canvas.draw_idle()

    def auto_scale(self):
        """Auto-scale axes to fit all points"""
        if not self.points:
            return

        points_array = np.array(self.points)
        margin = 2
        x_min, x_max = points_array[:, 0].min() - margin, points_array[:, 0].max() + margin
        y_min, y_max = points_array[:, 1].min() - margin, points_array[:, 1].max() + margin

        # Keep aspect ratio equal
        x_range = x_max - x_min
        y_range = y_max - y_min
        if x_range > y_range:
            y_center = (y_min + y_max) / 2
            y_min = y_center - x_range / 2
            y_max = y_center + x_range / 2
        else:
            x_center = (x_min + x_max) / 2
            x_min = x_center - y_range / 2
            x_max = x_center + y_range / 2

        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)

    def add_random_points(self, event):
        """Add random points"""
        n_points = 5
        x_min, x_max = self.ax.get_xlim()
        y_min, y_max = self.ax.get_ylim()

        for _ in range(n_points):
            x = int(np.random.uniform(x_min + 1, x_max - 1))
            y = int(np.random.uniform(y_min + 1, y_max - 1))
            self.add_point(x, y)

    def add_grid_points(self, event):
        """Add points in a grid pattern"""
        x_min, x_max = self.ax.get_xlim()
        y_min, y_max = self.ax.get_ylim()

        # Create a 3x3 grid in the center of the view
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2
        spacing = min(x_max - x_min, y_max - y_min) / 6

        for i in range(-1, 2):
            for j in range(-1, 2):
                x = int(x_center + i * spacing)
                y = int(y_center + j * spacing)
                self.add_point(x, y)

    def update_visibility(self, label):
        """Update visibility of plot elements based on checkboxes"""
        if label == 'Vertices':
            self.show_vertices = not self.show_vertices
        elif label == 'Vertex Numbers':
            self.show_vertex_numbers = not self.show_vertex_numbers
        elif label == 'Triangles':
            self.show_triangles = not self.show_triangles
        elif label == 'Triangle Numbers':
            self.show_triangle_numbers = not self.show_triangle_numbers
        elif label == 'Edges':
            self.show_edges = not self.show_edges
        elif label == 'Circumcircles':
            self.show_circumcircles = not self.show_circumcircles
        elif label == 'Boundary Highlight':
            self.show_boundary_highlight = not self.show_boundary_highlight

        self.update_plot()

    def compute_circumcircle(self, tri_points):
        """Compute circumcircle center and radius for a triangle"""
        try:
            p1, p2, p3 = tri_points

            ax, ay = p1
            bx, by = p2
            cx, cy = p3

            d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
            if abs(d) < 1e-10:
                return None, None

            ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) +
                  (cx * cx + cy * cy) * (ay - by)) / d
            uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) +
                  (cx * cx + cy * cy) * (bx - ax)) / d

            radius = np.sqrt((ax - ux) ** 2 + (ay - uy) ** 2)

            return (ux, uy), radius
        except:
            return None, None


def main():
    """Main function to run the GUI"""
    print("Delaunay Triangulation Interactive Visualizer")
    print("=" * 50)
    print("Instructions:")
    print("- Click on the plot to add points")
    print("- Use textboxes to insert points, enforce quality, or highlight a boundary")
    print("- Use the checkboxes to toggle display options (Boundary Highlight toggles the boundary view)")
    print("- Use buttons to clear, reset view, or add random/grid points")
    print("=" * 50)

    gui = DelaunayGUI()
    plt.show()


if __name__ == "__main__":
    main()
