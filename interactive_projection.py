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
from matplotlib.patches import Polygon as MplPolygon, Rectangle
from matplotlib.lines import Line2D
from matplotlib.widgets import Button, RadioButtons, Slider
import matplotlib.gridspec as gridspec
from polygon_projection import PolygonProjection, PolygonRegion, generate_regular_polygon
import threading
import time


class SimpleProgressBar:
    """A simple progress bar widget for matplotlib."""
    
    def __init__(self, ax, valmin=0, valmax=100):
        """
        Initialize a simple progress bar.
        
        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw the progress bar on.
        valmin : float, optional
            Minimum value of the progress bar. Default is 0.
        valmax : float, optional
            Maximum value of the progress bar. Default is 100.
        """
        self.ax = ax
        self.valmin = valmin
        self.valmax = valmax
        self.val = valmin
        
        # Create the progress bar rectangle
        self.rect = Rectangle((0, 0), 0, 1, color='blue', alpha=0.6)
        self.ax.add_patch(self.rect)
        
        # Create the background rectangle
        self.background = Rectangle((0, 0), 1, 1, color='gray', alpha=0.3)
        self.ax.add_patch(self.background)
        
        # Remove tick marks and spines
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for spine in self.ax.spines.values():
            spine.set_visible(False)
        
        # Set the limits
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        
        # Add a percentage label
        self.text = self.ax.text(0.5, 0.5, "0%", 
                                ha='center', va='center',
                                fontweight='bold')
    
    def set_val(self, val):
        """
        Set the value of the progress bar.
        
        Parameters
        ----------
        val : float
            Value between valmin and valmax to set the progress bar to.
        """
        self.val = max(self.valmin, min(self.valmax, val))
        
        # Convert to a fraction between 0 and 1
        fraction = (self.val - self.valmin) / (self.valmax - self.valmin)
        
        # Update the width of the rectangle
        self.rect.set_width(fraction)
        
        # Update the percentage text
        percentage = int(fraction * 100)
        self.text.set_text(f"{percentage}%")
        
        # Redraw the canvas
        self.ax.figure.canvas.draw_idle()


def is_convex(vertices):
    """
    Check if a polygon is convex by verifying that all interior angles are less than 180 degrees.
    Uses the cross product to check if all vertices are making "right turns" or all making "left turns".
    
    Parameters
    ----------
    vertices : array_like, shape (n,2)
        CCW-ordered polygon vertices.
        
    Returns
    -------
    bool
        True if the polygon is convex, False otherwise.
    """
    # Need at least 3 vertices to form a polygon
    if len(vertices) < 3:
        return True  # Degenerate case
        
    # Compute the cross product for each set of three consecutive vertices
    n = len(vertices)
    # Initialize sign of first cross product
    dx1 = vertices[1][0] - vertices[0][0]
    dy1 = vertices[1][1] - vertices[0][1]
    dx2 = vertices[2][0] - vertices[1][0]
    dy2 = vertices[2][1] - vertices[1][1]
    cross_z = dx1 * dy2 - dy1 * dx2
    sign = 1 if cross_z > 0 else -1 if cross_z < 0 else 0
    
    # Check all other cross products
    for i in range(1, n):
        dx1 = vertices[(i+1) % n][0] - vertices[i][0]
        dy1 = vertices[(i+1) % n][1] - vertices[i][1]
        dx2 = vertices[(i+2) % n][0] - vertices[(i+1) % n][0]
        dy2 = vertices[(i+2) % n][1] - vertices[(i+1) % n][1]
        cross_z = dx1 * dy2 - dy1 * dx2
        
        # If sign changes, the polygon is not convex
        curr_sign = 1 if cross_z > 0 else -1 if cross_z < 0 else 0
        if curr_sign != 0 and sign != 0 and curr_sign != sign:
            return False
        # Update sign if it was previously 0
        if sign == 0 and curr_sign != 0:
            sign = curr_sign
    
    return True


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
        
        # We no longer need a region selector as we'll show both interior and boundary values
        
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
        
        # Progress bar for Monte Carlo calculations
        progress_ax = plt.axes([0.25, 0.15, 0.35, 0.03])
        self.progress_bar = SimpleProgressBar(progress_ax, 0, 100)
        self.progress_bar.set_val(0)
        
        # Monte Carlo calculation button
        mc_button_ax = plt.axes([0.65, 0.15, 0.15, 0.05])
        self.mc_button = Button(mc_button_ax, 'Calculate Monte Carlo')
        self.mc_button.on_clicked(self.calculate_monte_carlo)
        
        # Stop button for Monte Carlo calculations
        stop_button_ax = plt.axes([0.82, 0.15, 0.08, 0.05])
        self.stop_button = Button(stop_button_ax, 'Stop')
        self.stop_button.on_clicked(self.stop_monte_carlo)
        
        # Calculation status
        self.calculating = False
        self.stop_requested = False
        
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
        
        # Monte Carlo results storage - includes both projection and average distances
        self.mc_results = {
            'interior_proj': None, 
            'boundary_proj': None,
            'interior_avg': None,
            'boundary_avg': None
        }
        
        # Exact average distance results
        self.avg_distances = {
            'interior': None,
            'boundary': None
        }
        
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
    
    def calculate_average_distances(self):
        """Calculate exact average distances for both interior and boundary."""
        if len(self.vertices) >= 3:
            # Interior
            interior_proj = PolygonProjection(
                self.vertices, self.direction, 
                region=PolygonRegion.INTERIOR
            )
            # Boundary
            boundary_proj = PolygonProjection(
                self.vertices, self.direction, 
                region=PolygonRegion.BOUNDARY
            )
            
            # Store the results
            self.avg_distances['interior'] = interior_proj.average_distance_exact()
            self.avg_distances['boundary'] = boundary_proj.average_distance_exact()
    
    def monte_carlo_worker(self):
        """Background worker thread to run Monte Carlo calculations with progress updates."""
        try:
            # Reset progress
            self.progress_bar.set_val(0)
            self.fig.canvas.draw_idle()
            
            # Get number of samples from logarithmic slider (10^slider_val)
            log_val = self.samples_slider.val
            self.n_samples = int(10**log_val)
            
            # Pre-allocate result containers
            batch_size = min(10000, max(100, self.n_samples // 20))  # 5% increments or at least 100
            num_batches = self.n_samples // batch_size
            
            # Initialize PolygonProjection objects
            interior_proj = PolygonProjection(
                self.vertices, self.direction, 
                region=PolygonRegion.INTERIOR,
                seed=42
            )
            boundary_proj = PolygonProjection(
                self.vertices, self.direction, 
                region=PolygonRegion.BOUNDARY,
                seed=43
            )
            
            # For averaging results across batches
            interior_proj_total = 0
            interior_avg_total = 0
            boundary_proj_total = 0
            boundary_avg_total = 0
            samples_processed = 0
            
            # Process in batches to update progress
            for i in range(num_batches):
                if self.stop_requested:
                    # Early termination if stop was requested
                    self.stop_requested = False
                    self.calculating = False
                    self.progress_bar.set_val(0)
                    return
                
                # Process a batch
                samples_in_batch = min(batch_size, self.n_samples - samples_processed)
                
                # Interior projections
                distances = interior_proj.monte_carlo_expected_distance(samples_in_batch)
                interior_proj_total += distances * samples_in_batch
                
                # Interior average distances
                avg_distances = interior_proj.average_distance_monte_carlo(samples_in_batch)
                interior_avg_total += avg_distances * samples_in_batch
                
                # Boundary projections
                distances = boundary_proj.monte_carlo_expected_distance(samples_in_batch)
                boundary_proj_total += distances * samples_in_batch
                
                # Boundary average distances
                avg_distances = boundary_proj.average_distance_monte_carlo(samples_in_batch)
                boundary_avg_total += avg_distances * samples_in_batch
                
                # Update progress and processed count
                samples_processed += samples_in_batch
                progress = int(100 * samples_processed / self.n_samples)
                
                # Update the progress bar safely from the worker thread
                def update_progress():
                    self.progress_bar.set_val(progress)
                    # Also update results as we go for real-time feedback
                    if i % 2 == 0:  # Update every other batch to avoid excessive redrawing
                        self.mc_results['interior_proj'] = interior_proj_total / samples_processed
                        self.mc_results['interior_avg'] = interior_avg_total / samples_processed
                        self.mc_results['boundary_proj'] = boundary_proj_total / samples_processed
                        self.mc_results['boundary_avg'] = boundary_avg_total / samples_processed
                        self.update_results()
                
                # Schedule the UI update on the main thread
                plt.gcf().canvas.flush_events()
                plt.gcf().canvas.draw_idle()
                update_progress()
                plt.gcf().canvas.flush_events()
                
                # Slight delay to allow UI updates
                time.sleep(0.01)
            
            # Store final results properly averaged
            if samples_processed > 0:
                self.mc_results['interior_proj'] = interior_proj_total / samples_processed
                self.mc_results['interior_avg'] = interior_avg_total / samples_processed
                self.mc_results['boundary_proj'] = boundary_proj_total / samples_processed
                self.mc_results['boundary_avg'] = boundary_avg_total / samples_processed
            
            # Final update
            self.progress_bar.set_val(100)
            self.update_results()
            
        except Exception as e:
            print(f"Monte Carlo calculation error: {e}")
        finally:
            # Always reset the calculation state
            self.calculating = False
            self.stop_requested = False
            # Enable all interactive elements
            self.enable_interactive_elements()
    
    def calculate_monte_carlo(self, event):
        """Start Monte Carlo calculations in a background thread with progress tracking."""
        if self.calculating:
            return  # Already calculating
        
        if len(self.vertices) >= 3 and np.linalg.norm(self.direction) > 0:
            # Set calculation state
            self.calculating = True
            self.stop_requested = False
            
            # Disable all interactive elements during calculation
            self.disable_interactive_elements()
            
            # Reset progress
            self.progress_bar.set_val(0)
            
            # Start calculation in a background thread
            thread = threading.Thread(target=self.monte_carlo_worker)
            thread.daemon = True  # Thread will be terminated when main program exits
            thread.start()
    
    def stop_monte_carlo(self, event):
        """Stop the ongoing Monte Carlo calculation."""
        if self.calculating:
            self.stop_requested = True
            # Note: The worker thread will handle cleanup and UI updates
    
    def disable_interactive_elements(self):
        """Disable all interactive elements during calculation."""
        self.mc_button.active = False
        self.reset_button.active = False
        self.add_point_button.active = False
        self.remove_point_button.active = False
        self.samples_slider.active = False
        # Color changes to indicate disabled state
        self.mc_button.color = 'lightgray'
        self.reset_button.color = 'lightgray'
        self.add_point_button.color = 'lightgray'
        self.remove_point_button.color = 'lightgray'
        # Refresh the figure
        self.fig.canvas.draw_idle()
    
    def enable_interactive_elements(self):
        """Re-enable all interactive elements after calculation completes."""
        self.mc_button.active = True
        self.reset_button.active = True
        self.add_point_button.active = True
        self.remove_point_button.active = True
        self.samples_slider.active = True
        # Reset colors
        self.mc_button.color = '0.85'
        self.reset_button.color = '0.85'
        self.add_point_button.color = '0.85'
        self.remove_point_button.color = '0.85'
        # Refresh the figure
        self.fig.canvas.draw_idle()
    
    def update_results(self):
        """Calculate and display projection and average distance results."""
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
                
                # Calculate exact expected distances for projection
                interior_exact_proj = interior_proj.exact_expected_distance()
                boundary_exact_proj = boundary_proj.exact_expected_distance()
                
                # Calculate or update average distances
                if self.avg_distances['interior'] is None or self.avg_distances['boundary'] is None:
                    self.calculate_average_distances()
                
                interior_exact_avg = self.avg_distances['interior']
                boundary_exact_avg = self.avg_distances['boundary']
                
                # Calculate polygon properties
                area = interior_proj.get_area()
                perimeter = interior_proj.get_perimeter()
                
                # Format results for display
                results_str = (
                    f"Polygon Info:\n"
                    f"- {len(self.vertices)} vertices\n"
                    f"- Area: {area:.4f}\n"
                    f"- Perimeter: {perimeter:.4f}\n"
                    f"- Direction: [{self.direction[0]:.3f}, {self.direction[1]:.3f}]\n\n"
                    f"Exact Results\n"
                    f"-------------\n"
                    f"Projected Distances:\n"
                    f"- Interior: {interior_exact_proj:.6f}\n"
                    f"- Boundary: {boundary_exact_proj:.6f}\n\n"
                    f"Average Distances:\n"
                    f"- Interior: {interior_exact_avg:.6f}\n"
                    f"- Boundary: {boundary_exact_avg:.6f}"
                )
                
                # Add Monte Carlo results if available
                has_mc_results = (
                    self.mc_results['interior_proj'] is not None and 
                    self.mc_results['boundary_proj'] is not None and
                    self.mc_results['interior_avg'] is not None and
                    self.mc_results['boundary_avg'] is not None
                )
                
                if has_mc_results:
                    # Get Monte Carlo values
                    interior_mc_proj = self.mc_results['interior_proj']
                    boundary_mc_proj = self.mc_results['boundary_proj']
                    interior_mc_avg = self.mc_results['interior_avg']
                    boundary_mc_avg = self.mc_results['boundary_avg']
                    
                    # Calculate errors for projected distances
                    interior_proj_error = abs(interior_exact_proj - interior_mc_proj) / interior_exact_proj
                    boundary_proj_error = abs(boundary_exact_proj - boundary_mc_proj) / boundary_exact_proj
                    
                    # Calculate errors for average distances
                    interior_avg_error = abs(interior_exact_avg - interior_mc_avg) / interior_exact_avg
                    boundary_avg_error = abs(boundary_exact_avg - boundary_mc_avg) / boundary_exact_avg
                    
                    results_str += (
                        f"\n\nMonte Carlo Results\n"
                        f"------------------\n"
                        f"Samples: {self.n_samples:,}\n\n"
                        f"Projected Distances:\n"
                        f"- Interior: {interior_mc_proj:.6f} (error: {interior_proj_error:.2%})\n"
                        f"- Boundary: {boundary_mc_proj:.6f} (error: {boundary_proj_error:.2%})\n\n"
                        f"Average Distances:\n"
                        f"- Interior: {interior_mc_avg:.6f} (error: {interior_avg_error:.2%})\n"
                        f"- Boundary: {boundary_mc_avg:.6f} (error: {boundary_avg_error:.2%})"
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
        self.mc_results = {
            'interior_proj': None, 
            'boundary_proj': None,
            'interior_avg': None,
            'boundary_avg': None
        }
        # Reset average distances
        self.avg_distances = {
            'interior': None,
            'boundary': None
        }
        self.update_results()
    
    def add_vertex(self, event):
        """
        Add a new vertex to the polygon while preserving convexity.
        The method attempts to find the edge with the longest length 
        and inserts a new vertex at its midpoint.
        """
        if len(self.vertices) < 20:  # Limit to 20 vertices
            if len(self.vertices) >= 3:
                # Find the longest edge to add a new vertex to
                n = len(self.vertices)
                max_length = 0
                max_idx = 0
                
                # Iterate through all edges to find the longest one
                for i in range(n):
                    p1 = self.vertices[i]
                    p2 = self.vertices[(i+1) % n]
                    length = np.linalg.norm(p2 - p1)
                    
                    if length > max_length:
                        max_length = length
                        max_idx = i
                
                # Create a new vertex at the midpoint of the longest edge
                p1 = self.vertices[max_idx]
                p2 = self.vertices[(max_idx+1) % n]
                new_vertex = (p1 + p2) / 2
                
                # Insert the new vertex between the endpoints of the longest edge
                new_vertices = np.insert(
                    self.vertices, 
                    (max_idx+1) % n, 
                    new_vertex, 
                    axis=0
                )
                
                # Ensure the result is convex (should always be the case with this method)
                if is_convex(new_vertices):
                    self.vertices = new_vertices
                    
                    # Update visualization
                    self.update_vertex_points()
                    self.update_polygon()
                    # Reset Monte Carlo results since polygon changed
                    self.mc_results = {
                        'interior_proj': None, 
                        'boundary_proj': None,
                        'interior_avg': None,
                        'boundary_avg': None
                    }
                    # Reset average distances
                    self.avg_distances = {
                        'interior': None,
                        'boundary': None
                    }
                    self.update_results()
                else:
                    print("Warning: Adding vertex failed to maintain convexity.")
            elif len(self.vertices) == 2:
                # For a line segment, just add a third point to make a triangle
                p1 = self.vertices[0]
                p2 = self.vertices[1]
                midpoint = (p1 + p2) / 2
                # Create a new point perpendicular to the line
                direction = p2 - p1
                perpendicular = np.array([-direction[1], direction[0]])  # 90-degree rotation
                perpendicular = perpendicular / np.linalg.norm(perpendicular) * 0.5  # normalize and scale
                
                new_vertex = midpoint + perpendicular
                self.vertices = np.vstack([self.vertices, new_vertex])
                
                # Update visualization
                self.update_vertex_points()
                self.update_polygon()
                # Reset results
                self.mc_results = {
                    'interior_proj': None, 
                    'boundary_proj': None,
                    'interior_avg': None,
                    'boundary_avg': None
                }
                self.avg_distances = {
                    'interior': None,
                    'boundary': None
                }
                self.update_results()
    
    def remove_vertex(self, event):
        """Remove the last vertex from the polygon."""
        if len(self.vertices) > 3:  # Maintain at least a triangle
            self.vertices = self.vertices[:-1]
            self.update_vertex_points()
            self.update_polygon()
            # Reset Monte Carlo results since polygon changed
            self.mc_results = {
                'interior_proj': None, 
                'boundary_proj': None,
                'interior_avg': None,
                'boundary_avg': None
            }
            # Reset average distances
            self.avg_distances = {
                'interior': None,
                'boundary': None
            }
            self.update_results()
    
    def on_press(self, event):
        """Handle mouse button press event."""
        # Block all interactions during Monte Carlo calculations
        if self.calculating:
            return
            
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
        # Block all interactions during Monte Carlo calculations
        if self.calculating:
            return
            
        if event.button != 1:
            return
        
        if self.selected_vertex is not None:
            # If vertex was being dragged, reset Monte Carlo results and average distances
            self.mc_results = {
                'interior_proj': None, 
                'boundary_proj': None,
                'interior_avg': None,
                'boundary_avg': None
            }
            # Reset average distances because polygon changed
            self.avg_distances = {
                'interior': None,
                'boundary': None
            }
            self.selected_vertex = None
            # Update results
            self.update_results()
        elif (self.dragging_direction_start or self.dragging_direction_end):
            # If direction was being changed, reset only Monte Carlo results
            # (average distances don't need to be recalculated since they're independent of direction)
            self.mc_results = {
                'interior_proj': None, 
                'boundary_proj': None,
                'interior_avg': None,
                'boundary_avg': None
            }
            self.dragging_direction_start = False
            self.dragging_direction_end = False
            # Update results
            self.update_results()
    
    def on_motion(self, event):
        """Handle mouse motion event."""
        # Block all interactions during Monte Carlo calculations
        if self.calculating:
            return
            
        if event.inaxes != self.ax_canvas:
            return
        
        update_needed = False
        reset_avg_distances = False
        
        if self.selected_vertex is not None:
            # Store the original vertex position in case we need to revert
            original_position = self.vertices[self.selected_vertex].copy()
            
            # Try updating the vertex position
            temp_vertices = self.vertices.copy()
            temp_vertices[self.selected_vertex] = [event.xdata, event.ydata]
            
            # Check if the new polygon would be convex
            if is_convex(temp_vertices):
                # Update vertex position since the result would be convex
                self.vertices[self.selected_vertex] = [event.xdata, event.ydata]
                self.vertex_points[self.selected_vertex].set_data(
                    [event.xdata], [event.ydata]
                )
                self.update_polygon()
                update_needed = True
                # When polygon vertices change, we need to reset average distances
                reset_avg_distances = True
            else:
                # Revert to original position
                self.vertices[self.selected_vertex] = original_position
        
        elif self.dragging_direction_start:
            # Update direction start point
            self.direction_start = np.array([event.xdata, event.ydata])
            self.update_direction()
            update_needed = True
            # Direction change doesn't require recalculating average distances
            
        elif self.dragging_direction_end:
            # Update direction end point
            self.direction_end = np.array([event.xdata, event.ydata])
            self.update_direction()
            update_needed = True
            # Direction change doesn't require recalculating average distances
            
        # Update results while dragging for real-time feedback
        if update_needed:
            # If polygon changes, reset average distances to force recalculation
            if reset_avg_distances:
                self.avg_distances = {
                    'interior': None,
                    'boundary': None
                }
            self.update_results()
    
    def show(self):
        """Display the interactive figure."""
        plt.show()


if __name__ == "__main__":
    app = InteractivePolygonProjection()
    app.show()