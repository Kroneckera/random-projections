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
import sys
sys.path.append('/Users/azimin/Programming/random_proj/polygon-projection/python')
from polygon_projection import Polygon, Region, ProjectionCalculator, is_convex
import threading
import time
import queue


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


# Using the is_convex function directly from the API


class InteractivePolygonProjection:
    """
    Interactive visualization tool for polygon projections.
    """
    
    def __init__(self):
        # Initialize figure and layout
        self.fig = plt.figure(figsize=(12, 8))
        self.fig.canvas.manager.set_window_title('Interactive Polygon Projection')
        
        # Create a GridSpec layout with original proportions
        gs = gridspec.GridSpec(3, 3)
        
        # Main canvas with increased size and centered position
        # Original size was 0.5 × 0.5, now 20% larger: 0.6 × 0.6
        canvas_width = 0.6
        canvas_height = 0.6
        
        # Calculate the position to center it in the left half
        left_half_width = 0.5  # Left half width (out of 1.0 total figure width)
        left_half_center = 0.25  # Center of left half (0.0 to 0.5)
        
        # Calculate position: centered horizontally in left half, vertically in upper area
        canvas_left = left_half_center - canvas_width/2
        canvas_bottom = 0.3  # Position vertically to fit in the figure
        
        # Create figure-level axes with the new position and size
        self.ax_canvas = self.fig.add_axes([canvas_left, canvas_bottom, canvas_width, canvas_height])  # [left, bottom, width, height]
        self.ax_canvas.set_xlim(-3, 3)
        self.ax_canvas.set_ylim(-3, 3)
        self.ax_canvas.set_aspect('equal')
        self.ax_canvas.grid(True)
        self.ax_canvas.set_title('Polygon Editor - Drag vertices or direction')
        
        # Results panel - adjust position to account for larger canvas
        # Move it to right and adjust size accordingly
        results_left = canvas_left + canvas_width + 0.02  # Position right after canvas with small gap
        results_width = 0.98 - results_left  # Use remaining space to right edge
        self.ax_results = self.fig.add_axes([results_left, canvas_bottom, results_width, canvas_height])  # [left, bottom, width, height]
        self.ax_results.axis('off')
        self.results_text = self.ax_results.text(
            0.05, 0.95, 'Results will appear here',
            verticalalignment='top', 
            transform=self.ax_results.transAxes
        )
        
        # Controls panel - using add_axes for consistent positioning
        self.ax_controls = self.fig.add_axes([0.05, 0.05, 0.9, 0.25])  # [left, bottom, width, height]
        self.ax_controls.axis('off')
        
        # Polygon initialization - standard regular hexagon
        polygon = Polygon.regular(sides=6, radius=2.0)
        self.vertices = polygon.vertices
        self.polygon_patch = MplPolygon(
            self.vertices, alpha=0.2, fc='b', ec='blue', zorder=1
        )
        self.ax_canvas.add_patch(self.polygon_patch)
        
        # Direction initialization - original values
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
        
        # Define controls layout with clear separation
        # Simple split at the middle (50/50)
        panel_width = 0.40    # Width for each panel (polygon and Monte Carlo)
        
        # Starting positions
        left_start = 0.05     # Left panel start
        right_start = 0.55    # Right panel start
        button_height = 0.04  # Standard button height
        control_spacing = 0.07  # Vertical spacing between controls
        
        # Polygon controls - left half, positioned further to the center
        reset_button_width = 0.10
        vertex_button_width = 0.10
        button_spacing = 0.01
        
        # Clear row positions used for all controls
        top_row_y = 0.22      # Top row of controls (slider)
        middle_row_y = 0.16   # Middle row (progress bar)
        bottom_row_y = 0.08   # Bottom row (buttons)
        
        # Position polygon buttons precisely in the center of the left half
        left_half_width = 0.45  # Width of the left half (from 0.05 to 0.5)
        total_polygon_buttons_width = reset_button_width + vertex_button_width + vertex_button_width + 3*button_spacing
        
        # Calculate position to center buttons beneath the polygon canvas
        # Use the canvas center position as reference
        canvas_center = canvas_left + canvas_width/2
        polygon_left_margin = canvas_center - total_polygon_buttons_width/2  # Center directly below canvas
        
        # All polygon buttons in a single row at the middle row height
        reset_ax = plt.axes([polygon_left_margin, middle_row_y, reset_button_width, button_height])
        self.reset_button = Button(reset_ax, 'Reset')
        self.reset_button.on_clicked(self.reset_polygon)
        
        # Add vertex button
        add_point_ax = plt.axes([polygon_left_margin + reset_button_width + button_spacing, 
                               middle_row_y, vertex_button_width, button_height])
        self.add_point_button = Button(add_point_ax, 'Add Vertex')
        self.add_point_button.on_clicked(self.add_vertex)
        
        # Remove vertex button
        remove_point_ax = plt.axes([polygon_left_margin + reset_button_width + vertex_button_width + 2*button_spacing, 
                                  middle_row_y, vertex_button_width, button_height])
        self.remove_point_button = Button(remove_point_ax, 'Remove Vertex')
        self.remove_point_button.on_clicked(self.remove_vertex)
        
        # Monte Carlo controls - right half
        mc_controls_width = 0.35
        
        # MC samples slider (logarithmic scale) - aligned with top row
        samples_ax = plt.axes([right_start, top_row_y, mc_controls_width, 0.03])
        self.samples_slider = Slider(
            samples_ax, 'MC Samples', 2, 6, 
            valinit=5
        )
        
        # Format function for showing actual sample count
        def format_func(x, pos):
            return f"{int(10**x):,}"
        
        self.samples_slider.valtext.set_text(format_func(5, None))
        
        # Update function for logarithmic scale
        def logarithmic_update(val):
            samples = int(10**val)
            self.samples_slider.valtext.set_text(format_func(val, None))
        
        self.samples_slider.on_changed(logarithmic_update)
        self.n_samples = 100000
        
        # Progress bar - aligned with middle row
        progress_ax = plt.axes([right_start, middle_row_y, mc_controls_width, 0.03])
        self.progress_bar = SimpleProgressBar(progress_ax, 0, 100)
        self.progress_bar.set_val(0)
        
        # Calculation buttons - aligned with bottom row
        calc_button_width = 0.25
        stop_button_width = 0.09
        
        mc_button_ax = plt.axes([right_start, bottom_row_y, calc_button_width, button_height])
        self.mc_button = Button(mc_button_ax, 'Calculate MC')
        self.mc_button.on_clicked(self.calculate_monte_carlo)
        
        stop_button_ax = plt.axes([right_start + calc_button_width + button_spacing, bottom_row_y, stop_button_width, button_height])
        self.stop_button = Button(stop_button_ax, 'Stop')
        self.stop_button.on_clicked(self.stop_monte_carlo)
        
        # Calculation status and communication
        self.calculating = False
        self.stop_requested = False
        self.results_queue = queue.Queue()
        
        # Set up timer for UI updates from worker thread
        self.timer = self.fig.canvas.new_timer(interval=100)  # 100ms interval
        self.timer.add_callback(self.process_results_queue)
        self.timer.start()
        
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
            # Create polygon
            polygon = Polygon(self.vertices)
            
            # Interior
            interior_calc = ProjectionCalculator(
                polygon=polygon,
                direction=self.direction,
                region=Region.INTERIOR
            )
            # Boundary
            boundary_calc = ProjectionCalculator(
                polygon=polygon,
                direction=self.direction,
                region=Region.BOUNDARY
            )
            
            # Store the results
            self.avg_distances['interior'] = interior_calc.average_distance(exact=True)
            self.avg_distances['boundary'] = boundary_calc.average_distance(exact=True)
    
    def process_results_queue(self):
        """Process any results in the queue and update the UI safely (called by timer on main thread)"""
        try:
            # Process all available messages in the queue
            while not self.results_queue.empty():
                msg = self.results_queue.get_nowait()
                
                if msg['type'] == 'progress':
                    # Update progress bar
                    self.progress_bar.set_val(msg['value'])
                
                elif msg['type'] == 'results':
                    # Update the results
                    self.mc_results['interior_proj'] = msg['interior_proj']
                    self.mc_results['interior_avg'] = msg['interior_avg']
                    self.mc_results['boundary_proj'] = msg['boundary_proj']
                    self.mc_results['boundary_avg'] = msg['boundary_avg']
                    self.update_results()
                
                elif msg['type'] == 'done':
                    # Calculation is complete
                    self.calculating = False
                    self.enable_interactive_elements()
                    
                elif msg['type'] == 'error':
                    # Handle error
                    print(f"Monte Carlo calculation error: {msg['error']}")
                    self.calculating = False
                    self.progress_bar.set_val(0)
                    self.enable_interactive_elements()
                
                self.results_queue.task_done()
                
        except Exception as e:
            print(f"Error processing results queue: {e}")
    
    def monte_carlo_worker(self):
        """Background worker thread to run Monte Carlo calculations using generators for better progress reporting."""
        try:
            # Get number of samples from logarithmic slider (10^slider_val)
            log_val = self.samples_slider.val
            self.n_samples = int(10**log_val)
            
            # Create polygon
            polygon = Polygon(self.vertices.copy())
            
            # Initialize ProjectionCalculator objects
            interior_calc = ProjectionCalculator(
                polygon=polygon,
                direction=self.direction.copy(),
                region=Region.INTERIOR,
                random_seed=42
            )
            boundary_calc = ProjectionCalculator(
                polygon=polygon,
                direction=self.direction.copy(),
                region=Region.BOUNDARY,
                random_seed=43
            )
            
            # Calculate batch size based on total samples
            batch_size = min(1000, max(100, self.n_samples // 100))  # Dynamic batch size
            yield_interval = 2  # Yield after processing every 2 batches
            
            # Create generators for interior and boundary calculations
            interior_generator = interior_calc.monte_carlo_analysis(
                samples=self.n_samples, 
                batch_size=batch_size,
                yield_interval=yield_interval
            )
            boundary_generator = boundary_calc.monte_carlo_analysis(
                samples=self.n_samples,
                batch_size=batch_size,
                yield_interval=yield_interval
            )
            
            # Process both generators in parallel - interleaving their results
            interior_done = False
            boundary_done = False
            interior_results = None
            boundary_results = None
            
            # Loop until both generators are done
            while not (interior_done and boundary_done):
                if self.stop_requested:
                    # Early termination if stop was requested
                    break
                
                # Get next results from interior generator
                if not interior_done:
                    try:
                        interior_results = next(interior_generator)
                    except StopIteration:
                        interior_done = True
                
                # Get next results from boundary generator
                if not boundary_done:
                    try:
                        boundary_results = next(boundary_generator)
                    except StopIteration:
                        boundary_done = True
                
                # Update UI if we have results from both
                if interior_results and boundary_results:
                    # Calculate overall progress as average of both generators
                    total_processed = (interior_results['samples_processed'] + 
                                      boundary_results['samples_processed'])
                    total_samples = (interior_results['total_samples'] + 
                                    boundary_results['total_samples'])
                    progress = int(100 * total_processed / total_samples)
                    
                    # Queue the new progress for UI update
                    self.results_queue.put({
                        'type': 'progress',
                        'value': progress
                    })
                    
                    # Queue the current calculation results for UI update
                    self.results_queue.put({
                        'type': 'results',
                        'interior_proj': interior_results['projected_distance'],
                        'interior_avg': interior_results['average_distance'],
                        'boundary_proj': boundary_results['projected_distance'],
                        'boundary_avg': boundary_results['average_distance']
                    })
                
                # Small sleep to prevent CPU hogging
                time.sleep(0.01)
            
            # Final progress update if not stopped
            if not self.stop_requested:
                # Ensure progress bar shows 100%
                self.results_queue.put({
                    'type': 'progress',
                    'value': 100
                })
            
        except Exception as e:
            # Queue error message
            self.results_queue.put({'type': 'error', 'error': str(e)})
        finally:
            # Always queue 'done' message to clean up
            self.results_queue.put({'type': 'done'})
            self.stop_requested = False
    
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
            # Set stop flag for worker thread
            self.stop_requested = True
            # Change the button appearance to indicate stop request is in process
            self.stop_button.label.set_text("Stopping...")
            self.fig.canvas.draw_idle()
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
        # Reset stop button text
        self.stop_button.label.set_text("Stop")
        # Refresh the figure
        self.fig.canvas.draw_idle()
    
    def update_results(self):
        """Calculate and display projection and average distance results."""
        try:
            # If we have at least 3 vertices and a valid direction
            if len(self.vertices) >= 3 and np.linalg.norm(self.direction) > 0:
                # Create polygon
                polygon = Polygon(self.vertices)
                
                # Create projection calculators for both regions
                interior_calc = ProjectionCalculator(
                    polygon=polygon,
                    direction=self.direction,
                    region=Region.INTERIOR
                )
                boundary_calc = ProjectionCalculator(
                    polygon=polygon,
                    direction=self.direction,
                    region=Region.BOUNDARY
                )
                
                # Calculate exact expected distances for projection
                interior_exact_proj = interior_calc.projected_distance(exact=True)
                boundary_exact_proj = boundary_calc.projected_distance(exact=True)
                
                # Calculate or update average distances
                if self.avg_distances['interior'] is None or self.avg_distances['boundary'] is None:
                    self.calculate_average_distances()
                
                interior_exact_avg = self.avg_distances['interior']
                boundary_exact_avg = self.avg_distances['boundary']
                
                # Calculate polygon properties
                area = polygon.area
                perimeter = polygon.perimeter
                
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
        polygon = Polygon.regular(sides=6, radius=2.0)
        self.vertices = polygon.vertices
        
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
                # Using the API function to check convexity
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
            # Using the API function to check convexity
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