#!/usr/bin/env python3
"""
interactive_projection.py

An interactive matplotlib-based application for visualizing polygon projections.
This tool allows users to:
- Create and manipulate convex polygons on a canvas by dragging vertices
- Adjust the projection direction
- See real-time updates of exact expected distance calculations
- Toggle between interior and boundary projections
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.lines import Line2D
from matplotlib.widgets import Button, RadioButtons, Slider
import matplotlib.gridspec as gridspec
from polygon_projection import PolygonProjection, PolygonRegion, generate_regular_polygon


class InteractivePolygonProjection:
    """
    Interactive visualization tool for polygon projections.
    """
    
    def __init__(self):
        # Initialize figure and layout
        self.fig = plt.figure(figsize=(12, 8))
        self.fig.canvas.manager.set_window_title('Interactive Polygon Projection')
        gs = gridspec.GridSpec(3, 3)
        
        # Main canvas
        self.ax_canvas = self.fig.add_subplot(gs[:2, :2])
        self.ax_canvas.set_xlim(-3, 3)
        self.ax_canvas.set_ylim(-3, 3)
        self.ax_canvas.set_aspect('equal')
        self.ax_canvas.grid(True)
        self.ax_canvas.set_title('Polygon Editor - Drag vertices or direction')
        
        # Results panel
        self.ax_results = self.fig.add_subplot(gs[:2, 2])
        self.ax_results.axis('off')
        self.results_text = self.ax_results.text(
            0.05, 0.95, 'Results will appear here',
            verticalalignment='top', 
            transform=self.ax_results.transAxes
        )
        
        # Controls panel
        self.ax_controls = self.fig.add_subplot(gs[2, :])
        self.ax_controls.axis('off')
        
        # Polygon initialization
        self.vertices = generate_regular_polygon(n_sides=6, radius=2.0)
        self.polygon_patch = MplPolygon(
            self.vertices, alpha=0.2, fc='b', ec='blue', zorder=1
        )
        self.ax_canvas.add_patch(self.polygon_patch)
        
        # Direction initialization
        self.direction_start = np.array([0.0, 0.0])
        self.direction_end = np.array([1.0, 0.5])
        self.direction = self.direction_end - self.direction_start
        
        # Direction line
        self.direction_line = Line2D(
            [self.direction_start[0], self.direction_end[0]], 
            [self.direction_start[1], self.direction_end[1]], 
            color='red', linewidth=2, zorder=2
        )
        self.ax_canvas.add_line(self.direction_line)
        
        # Vertex points for dragging
        self.vertex_points = []
        self.update_vertex_points()
        
        # Direction endpoints for dragging
        self.direction_start_point = Line2D(
            [self.direction_start[0]], [self.direction_start[1]], 
            marker='o', color='red', markersize=8, zorder=3
        )
        self.direction_end_point = Line2D(
            [self.direction_end[0]], [self.direction_end[1]], 
            marker='o', color='red', markersize=8, zorder=3
        )
        self.ax_canvas.add_line(self.direction_start_point)
        self.ax_canvas.add_line(self.direction_end_point)
        
        # Projection axis (parallel to the direction vector but passing through the segment)
        # Calculate initial position of the projection axis
        norm = np.linalg.norm(self.direction)
        if norm > 0:
            normalized = self.direction / norm
            extension = 3.0
            midpoint = (self.direction_start + self.direction_end) / 2
            
            # Initialize projection axis
            self.projection_axis = Line2D(
                [midpoint[0] - extension * normalized[0], midpoint[0] + extension * normalized[0]], 
                [midpoint[1] - extension * normalized[1], midpoint[1] + extension * normalized[1]],
                color='black', alpha=0.7, linestyle='-', zorder=1
            )
        else:
            self.projection_axis = Line2D([0, 0], [0, 0], color='black', alpha=0.7, zorder=1)
            
        self.ax_canvas.add_line(self.projection_axis)
        
        # Initialize RadioButtons for region selection
        region_ax = plt.axes([0.05, 0.2, 0.15, 0.1])
        self.region_selector = RadioButtons(
            region_ax, ('Interior', 'Boundary'),
            activecolor='blue'
        )
        self.region_selector.on_clicked(self.update_region)
        self.region = PolygonRegion.INTERIOR
        
        # Monte Carlo samples slider (logarithmic scale)
        samples_ax = plt.axes([0.55, 0.25, 0.35, 0.03])
        # Use log scale from 10^2 to 10^6
        self.samples_slider = Slider(
            samples_ax, 'MC Samples', 2, 6, 
            valinit=5
        )
        
        # Override the format string dynamically
        def format_func(x, pos):
            return f"{int(10**x):,}"
        
        self.samples_slider.valtext.set_text(format_func(5, None))
        
        # Function to convert from log slider to actual value
        def logarithmic_update(val):
            # Convert from slider value to number of samples (10^val)
            samples = int(10**val)
            # Update display text with formatted value
            self.samples_slider.valtext.set_text(format_func(val, None))
            
        # Connect log update function
        self.samples_slider.on_changed(logarithmic_update)
        self.n_samples = 100000
        
        # Monte Carlo calculation button
        mc_button_ax = plt.axes([0.65, 0.15, 0.25, 0.05])
        self.mc_button = Button(mc_button_ax, 'Calculate Monte Carlo')
        self.mc_button.on_clicked(self.calculate_monte_carlo)
        
        # Buttons
        reset_ax = plt.axes([0.05, 0.05, 0.15, 0.05])
        self.reset_button = Button(reset_ax, 'Reset Polygon')
        self.reset_button.on_clicked(self.reset_polygon)
        
        add_point_ax = plt.axes([0.25, 0.05, 0.15, 0.05])
        self.add_point_button = Button(add_point_ax, 'Add Vertex')
        self.add_point_button.on_clicked(self.add_vertex)
        
        remove_point_ax = plt.axes([0.45, 0.05, 0.15, 0.05])
        self.remove_point_button = Button(remove_point_ax, 'Remove Vertex')
        self.remove_point_button.on_clicked(self.remove_vertex)
        
        # Monte Carlo results storage
        self.mc_results = {'interior': None, 'boundary': None}
        
        # Initialize interaction
        self.selected_vertex = None
        self.dragging_direction_start = False
        self.dragging_direction_end = False
        
        # Connect events
        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        
        # Initial results calculation
        self.update_results()
        
        # Add a title for the application
        self.fig.suptitle('Interactive Polygon Projection Calculator', fontsize=16)
        plt.subplots_adjust(top=0.9, bottom=0.1, left=0.05, right=0.95)
    
    def update_vertex_points(self):
        """Update the visual representation of vertices."""
        for p in self.vertex_points:
            p.remove()
        
        self.vertex_points = []
        for x, y in self.vertices:
            point = Line2D([x], [y], marker='o', color='blue', markersize=8, zorder=3)
            self.ax_canvas.add_line(point)
            self.vertex_points.append(point)
    
    def update_polygon(self):
        """Update the polygon patch with current vertices."""
        self.polygon_patch.set_xy(self.vertices)
        self.fig.canvas.draw_idle()
    
    def update_direction(self):
        """Update the direction vector and all related visual elements."""
        # Update direction vector
        self.direction = self.direction_end - self.direction_start
        
        # Update direction line
        self.direction_line.set_data(
            [self.direction_start[0], self.direction_end[0]], 
            [self.direction_start[1], self.direction_end[1]]
        )
        
        # Update direction endpoints
        self.direction_start_point.set_data([self.direction_start[0]], [self.direction_start[1]])
        self.direction_end_point.set_data([self.direction_end[0]], [self.direction_end[1]])
        
        # Update the projection axis to align with the segment
        norm = np.linalg.norm(self.direction)
        if norm > 0:
            # Calculate a perpendicular vector to the direction
            # This is used to extend the projection axis along the direction
            normalized = self.direction / norm
            
            # Calculate extension distance (6 units total)
            extension = 3.0
            
            # Get the midpoint of the direction segment
            midpoint = (self.direction_start + self.direction_end) / 2
            
            # Set the projection axis to pass through the midpoint of the segment
            # and extend in both directions
            self.projection_axis.set_data(
                [midpoint[0] - extension * normalized[0], 
                 midpoint[0] + extension * normalized[0]], 
                [midpoint[1] - extension * normalized[1], 
                 midpoint[1] + extension * normalized[1]]
            )
        else:
            # If direction is zero, draw a point
            self.projection_axis.set_data([0, 0], [0, 0])
        
        self.fig.canvas.draw_idle()
    
    def update_region(self, label):
        """Update the selected region."""
        self.region = PolygonRegion.INTERIOR if label == 'Interior' else PolygonRegion.BOUNDARY
        self.update_results()
    
    def calculate_monte_carlo(self, event):
        """Calculate Monte Carlo approximations for both regions."""
        try:
            if len(self.vertices) >= 3 and np.linalg.norm(self.direction) > 0:
                # Get number of samples from logarithmic slider (10^slider_val)
                log_val = self.samples_slider.val
                self.n_samples = int(10**log_val)
                
                # Calculate for interior
                interior_proj = PolygonProjection(
                    self.vertices, self.direction, 
                    region=PolygonRegion.INTERIOR,
                    seed=42
                )
                interior_mc = interior_proj.monte_carlo_expected_distance(self.n_samples)
                self.mc_results['interior'] = interior_mc
                
                # Calculate for boundary
                boundary_proj = PolygonProjection(
                    self.vertices, self.direction, 
                    region=PolygonRegion.BOUNDARY,
                    seed=42
                )
                boundary_mc = boundary_proj.monte_carlo_expected_distance(self.n_samples)
                self.mc_results['boundary'] = boundary_mc
                
                # Update display
                self.update_results()
                
        except Exception as e:
            print(f"Monte Carlo calculation error: {e}")
    
    def update_results(self):
        """Calculate and display exact projection results."""
        try:
            # If we have at least 3 vertices and a valid direction
            if len(self.vertices) >= 3 and np.linalg.norm(self.direction) > 0:
                # Create projection objects for both regions
                interior_proj = PolygonProjection(
                    self.vertices, self.direction, 
                    region=PolygonRegion.INTERIOR
                )
                boundary_proj = PolygonProjection(
                    self.vertices, self.direction, 
                    region=PolygonRegion.BOUNDARY
                )
                
                # Calculate exact expected distances
                interior_exact = interior_proj.exact_expected_distance()
                boundary_exact = boundary_proj.exact_expected_distance()
                
                # Get currently selected value
                current_exact = interior_exact if self.region == PolygonRegion.INTERIOR else boundary_exact
                
                # Format results for display
                results_str = (
                    f"Polygon Info:\n"
                    f"- {len(self.vertices)} vertices\n"
                    f"- Direction: [{self.direction[0]:.3f}, {self.direction[1]:.3f}]\n\n"
                    f"Exact Expected Distances (E[|X-Y|]):\n\n"
                    f"Interior: {interior_exact:.6f}\n\n"
                    f"Boundary: {boundary_exact:.6f}\n\n"
                    f"Selected: {self.region.value.capitalize()}\n"
                    f"Value: {current_exact:.6f}"
                )
                
                # Add Monte Carlo results if available
                if self.mc_results['interior'] is not None and self.mc_results['boundary'] is not None:
                    interior_mc = self.mc_results['interior']
                    boundary_mc = self.mc_results['boundary']
                    
                    # Get current Monte Carlo result
                    current_mc = interior_mc if self.region == PolygonRegion.INTERIOR else boundary_mc
                    
                    # Calculate errors
                    interior_error = abs(interior_exact - interior_mc) / interior_exact
                    boundary_error = abs(boundary_exact - boundary_mc) / boundary_exact
                    
                    results_str += (
                        f"\n\n====== Monte Carlo Results ======\n"
                        f"Samples: {self.n_samples:,}\n\n"
                        f"Interior: {interior_mc:.6f} (error: {interior_error:.2%})\n\n"
                        f"Boundary: {boundary_mc:.6f} (error: {boundary_error:.2%})"
                    )
                
            else:
                results_str = (
                    "Need at least 3 vertices and\n"
                    "a non-zero direction vector."
                )
                
        except Exception as e:
            results_str = f"Error calculating results:\n{str(e)}"
        
        # Update results text
        self.results_text.set_text(results_str)
        self.fig.canvas.draw_idle()
    
    def reset_polygon(self, event):
        """Reset to a regular hexagon."""
        self.vertices = generate_regular_polygon(n_sides=6, radius=2.0)
        self.update_vertex_points()
        self.update_polygon()
        # Reset Monte Carlo results
        self.mc_results = {'interior': None, 'boundary': None}
        self.update_results()
    
    def add_vertex(self, event):
        """Add a new vertex to the polygon."""
        if len(self.vertices) < 20:  # Limit to 20 vertices
            # Add new vertex at midpoint of the last and first vertex
            if len(self.vertices) >= 2:
                last = self.vertices[-1]
                first = self.vertices[0]
                new_vertex = (last + first) / 2
                
                # Insert at the end (before closing the loop)
                self.vertices = np.vstack([self.vertices, new_vertex])
                
                # Update visualization
                self.update_vertex_points()
                self.update_polygon()
                # Reset Monte Carlo results since polygon changed
                self.mc_results = {'interior': None, 'boundary': None}
                self.update_results()
    
    def remove_vertex(self, event):
        """Remove the last vertex from the polygon."""
        if len(self.vertices) > 3:  # Maintain at least a triangle
            self.vertices = self.vertices[:-1]
            self.update_vertex_points()
            self.update_polygon()
            # Reset Monte Carlo results since polygon changed
            self.mc_results = {'interior': None, 'boundary': None}
            self.update_results()
    
    def on_press(self, event):
        """Handle mouse button press event."""
        if event.inaxes != self.ax_canvas or event.button != 1:
            return
        
        # Check if a vertex was clicked
        for i, point in enumerate(self.vertex_points):
            contains, _ = point.contains(event)
            if contains:
                self.selected_vertex = i
                return
        
        # Check if direction start point was clicked
        contains, _ = self.direction_start_point.contains(event)
        if contains:
            self.dragging_direction_start = True
            return
            
        # Check if direction end point was clicked
        contains, _ = self.direction_end_point.contains(event)
        if contains:
            self.dragging_direction_end = True
            return
    
    def on_release(self, event):
        """Handle mouse button release event."""
        if event.button != 1:
            return
        
        if self.selected_vertex is not None:
            # If vertex was being dragged, reset Monte Carlo results
            self.mc_results = {'interior': None, 'boundary': None}
            self.selected_vertex = None
            # Update results
            self.update_results()
        elif (self.dragging_direction_start or self.dragging_direction_end):
            # If direction was being changed, reset Monte Carlo results
            self.mc_results = {'interior': None, 'boundary': None}
            self.dragging_direction_start = False
            self.dragging_direction_end = False
            # Update results
            self.update_results()
    
    def on_motion(self, event):
        """Handle mouse motion event."""
        if event.inaxes != self.ax_canvas:
            return
        
        update_needed = False
        
        if self.selected_vertex is not None:
            # Update vertex position
            self.vertices[self.selected_vertex] = [event.xdata, event.ydata]
            self.vertex_points[self.selected_vertex].set_data(
                [event.xdata], [event.ydata]
            )
            self.update_polygon()
            update_needed = True
        
        elif self.dragging_direction_start:
            # Update direction start point
            self.direction_start = np.array([event.xdata, event.ydata])
            self.update_direction()
            update_needed = True
            
        elif self.dragging_direction_end:
            # Update direction end point
            self.direction_end = np.array([event.xdata, event.ydata])
            self.update_direction()
            update_needed = True
            
        # Update results while dragging for real-time feedback
        if update_needed:
            self.update_results()
    
    def show(self):
        """Display the interactive figure."""
        plt.show()


if __name__ == "__main__":
    app = InteractivePolygonProjection()
    app.show()