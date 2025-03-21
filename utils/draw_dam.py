import argparse
from PIL import Image, ImageDraw, ImageTk
import tkinter as tk
from tkinter import Canvas, Button, Frame, Scale, HORIZONTAL

class ImageEditor:
    def __init__(self, image_path, threshold):
        self.threshold = threshold
        self.image_path = image_path
        self.image = Image.open(image_path)
        self.original_image = self.image.copy()  # Store original image for reset
        self.draw = ImageDraw.Draw(self.image)
        self.root = tk.Tk()
        self.root.title("Dam Editor")
        
        # Scale factor for display (doesn't affect actual image)
        self.scale_factor = 1.0
        
        # Create frame for canvas
        self.frame = Frame(self.root)
        self.frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        # Create canvas with adjusted size for scaling
        canvas_width = min(self.image.width * self.scale_factor, 1200)
        canvas_height = min(self.image.height * self.scale_factor, 800)
        self.canvas = Canvas(self.frame, width=canvas_width, height=canvas_height)
        self.canvas.pack(fill=tk.BOTH, expand=True)
        
        # Load and display the image
        self.update_canvas()
        
        # Bind events
        self.canvas.bind("<Button-1>", self.on_click)
        # Add mousewheel zoom
        self.canvas.bind("<MouseWheel>", self.on_mousewheel)  # Windows
        self.canvas.bind("<Button-4>", self.on_mousewheel)    # Linux scroll up
        self.canvas.bind("<Button-5>", self.on_mousewheel)    # Linux scroll down
        
        # Create control panel
        self.control_panel = Frame(self.root)
        self.control_panel.pack(side=tk.BOTTOM, fill=tk.X)
        
        # Add scale slider
        self.scale_label = tk.Label(self.control_panel, text="Zoom:")
        self.scale_label.pack(side=tk.LEFT, padx=5, pady=5)
        self.scale_slider = Scale(self.control_panel, from_=0.5, to=5.0, resolution=0.1,
                                 orient=HORIZONTAL, command=self.update_scale, length=200)
        self.scale_slider.set(self.scale_factor)
        self.scale_slider.pack(side=tk.LEFT, padx=5, pady=5)
        
        # Add buttons for different tools
        self.fill_button = Button(self.control_panel, text="Fill Mode", command=self.toggle_fill_mode)
        self.fill_button.pack(side=tk.LEFT, padx=5, pady=5)
        
        self.reset_button = Button(self.control_panel, text="Reset Image", command=self.reset_image)
        self.reset_button.pack(side=tk.LEFT, padx=5, pady=5)
        
        self.save_button = Button(self.control_panel, text="Save Image", command=self.save_current_image)
        self.save_button.pack(side=tk.RIGHT, padx=5, pady=5)
        
        self.points = []
        self.fill_mode = False
        self.output_path = None

    def on_click(self, event):
        # Convert display coordinates back to image coordinates
        img_x = int(event.x / self.scale_factor)
        img_y = int(event.y / self.scale_factor)
        
        if self.fill_mode:
            self.flood_fill(img_x, img_y, (255, 0, 0, 255))  # Red color for fill
        else:
            self.points.append((img_x, img_y))
            # Draw a point marker on the canvas (not on the image)
            self.canvas.create_oval(event.x-3, event.y-3, event.x+3, event.y+3, fill='red')
            
            if len(self.points) == 2:
                self.draw_line()
                self.points = []

    def toggle_fill_mode(self):
        self.fill_mode = not self.fill_mode
        if self.fill_mode:
            self.fill_button.config(text="Line Mode")
            self.canvas.config(cursor="crosshair")
        else:
            self.fill_button.config(text="Fill Mode")
            self.canvas.config(cursor="arrow")

    def on_mousewheel(self, event):
        """Handle mousewheel events for zooming"""
        delta = 0
        # Different systems report scroll events differently
        if event.num == 4 or event.delta > 0:
            delta = 0.1
        elif event.num == 5 or event.delta < 0:
            delta = -0.1
            
        new_scale = max(0.5, min(5.0, self.scale_factor + delta))
        self.scale_slider.set(new_scale)
        self.update_scale(new_scale)

    def update_scale(self, value):
        """Update the scale factor and redraw"""
        self.scale_factor = float(value)
        self.update_canvas()

    def draw_line(self):
        self.draw.line(self.points, fill="green", width=2)
        self.update_canvas()

    def flood_fill(self, x, y, fill_color):
        # Convert x,y to PIL image coordinates in case of canvas scrolling
        width, height = self.image.size
        if x < 0 or y < 0 or x >= width or y >= height:
            return  # Out of bounds

        # Get the color at the clicked position
        try:
            old_color = self.image.getpixel((x, y))
        except IndexError:
            return
        
        if old_color == fill_color:
            return

        # Use PIL's flood fill instead of recursive implementation to avoid stack overflow
        ImageDraw.floodfill(self.image, (x, y), fill_color, thresh=self.threshold)
        self.update_canvas()

    def update_canvas(self):
        # Create a resized version for display
        display_width = int(self.image.width * self.scale_factor)
        display_height = int(self.image.height * self.scale_factor)
        
        # Resize for display only (temporary)
        display_image = self.image.resize((display_width, display_height), Image.LANCZOS)
        
        # Convert PIL Image to PhotoImage
        self.photo = ImageTk.PhotoImage(display_image)
        
        # Adjust canvas size if needed
        canvas_width = min(display_width, 1200)
        canvas_height = min(display_height, 800)
        self.canvas.config(width=canvas_width, height=canvas_height)
        
        # Update the image on canvas
        self.canvas.delete("all")
        self.canvas.create_image(0, 0, anchor=tk.NW, image=self.photo)
        self.canvas.config(scrollregion=(0, 0, display_width, display_height))

    def save_current_image(self):
        if self.output_path:
            self.image.save(self.output_path)
            print(f"Image saved to {self.output_path}")
            self.root.quit()

    def save_image(self, output_path):
        self.output_path = output_path

    def reset_image(self):
        """Reset the image to its original state."""
        self.image = self.original_image.copy()
        self.draw = ImageDraw.Draw(self.image)
        self.update_canvas()
        self.points = []  # Clear any selected points

    def run(self):
        self.root.mainloop()

def main():
    parser = argparse.ArgumentParser(description="Image Editor")
    parser.add_argument("input_image", help="Path to the input image")
    parser.add_argument("output_image", help="Path to save the edited image")
    parser.add_argument("threshold", help="Flood fill threshold")
    args = parser.parse_args()

    editor = ImageEditor(args.input_image, int(args.threshold))
    editor.save_image(args.output_image)
    editor.run()

if __name__ == "__main__":
    main()