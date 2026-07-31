import tkinter as tk
from tkinter import ttk, filedialog
from tkinterdnd2 import DND_FILES, TkinterDnD
from PIL import Image, ImageTk
import os
import json
from pathlib import Path

class FrameLabeler:
    def __init__(self, frames_dir=None):
        # Create GUI with drag-and-drop support
        self.root = TkinterDnD.Tk()
        self.root.title("NORT Frame Labeler")
        self.root.geometry("1200x900")
        
        self.frames_dir = None
        self.output_file = None
        self.image_files = []
        self.current_idx = 0
        self.labels = {}
        
        if frames_dir:
            self.load_directory(frames_dir)
        else:
            self.show_folder_selection()
    
    def show_folder_selection(self):
        """Show initial screen for folder selection"""
        self.selection_frame = tk.Frame(self.root, bg='#f0f0f0')
        self.selection_frame.pack(fill=tk.BOTH, expand=True)
        
        # Title
        tk.Label(self.selection_frame, text="NORT Frame Labeler", 
                font=('Arial', 24, 'bold'), bg='#f0f0f0').pack(pady=30)
        
        # Drag and drop area
        drop_frame = tk.Frame(self.selection_frame, bg='white', 
                             relief=tk.SOLID, borderwidth=2)
        drop_frame.pack(pady=20, padx=50, fill=tk.BOTH, expand=True)
        
        drop_label = tk.Label(drop_frame, 
                             text="Drag and Drop Frames Folder Here\n\n📁\n\nor", 
                             font=('Arial', 16), bg='white', fg='#666')
        drop_label.pack(pady=50)
        
        # Browse button
        ttk.Button(drop_frame, text="Browse for Folder", 
                  command=self.browse_folder).pack(pady=20)
        
        # Enable drag and drop
        drop_frame.drop_target_register(DND_FILES)
        drop_frame.dnd_bind('<<Drop>>', self.on_drop)
        drop_label.drop_target_register(DND_FILES)
        drop_label.dnd_bind('<<Drop>>', self.on_drop)
    
    def on_drop(self, event):
        """Handle drag and drop event"""
        path = event.data
        # Clean up the path (remove curly braces if present)
        path = path.strip('{}')
        
        if os.path.isdir(path):
            self.load_directory(path)
        else:
            tk.messagebox.showerror("Error", "Please drop a folder, not a file")
    
    def browse_folder(self):
        """Open folder browser dialog"""
        folder = filedialog.askdirectory(title="Select Frames Folder")
        if folder:
            self.load_directory(folder)
    
    def load_directory(self, frames_dir):
        """Load the frames directory and set up the labeler"""
        self.frames_dir = Path(frames_dir)
        
        # Create output filename based on folder name
        folder_name = self.frames_dir.name
        self.output_file = self.frames_dir / f"{folder_name}_labels.json"
        
        self.image_files = sorted([f for f in os.listdir(frames_dir) 
                                  if f.endswith(('.jpg', '.png', '.jpeg'))])
        
        if not self.image_files:
            tk.messagebox.showerror("Error", "No image files found in the selected folder")
            return
        
        self.current_idx = 0
        self.labels = {}
        
        # Load existing labels if available
        if os.path.exists(self.output_file):
            with open(self.output_file, 'r') as f:
                self.labels = json.load(f)
            print(f"Loaded {len(self.labels)} existing labels")
        
        # Remove selection screen if it exists
        if hasattr(self, 'selection_frame'):
            self.selection_frame.destroy()
        
        # Set up the labeling interface
        self.setup_labeling_interface()
        
    def setup_labeling_interface(self):
        """Create the labeling interface"""
        # Bind keyboard shortcuts
        self.root.bind('e', lambda event: self.label_frame('exploratory'))
        self.root.bind('n', lambda event: self.label_frame('non-exploratory'))
        self.root.bind('u', lambda event: self.label_frame('uncertain'))
        self.root.bind('<Left>', lambda event: self.prev_frame())
        self.root.bind('<Right>', lambda event: self.next_frame())
        self.root.bind('s', lambda event: self.save_labels())
        
        # Folder info label
        folder_info = tk.Label(self.root, 
                              text=f"Working on: {self.frames_dir}", 
                              font=('Arial', 10), fg='#666')
        folder_info.pack(pady=5)
        
        # Progress label
        self.progress_label = tk.Label(self.root, text="", font=('Arial', 12))
        self.progress_label.pack(pady=10)
        
        # Image display
        self.image_label = tk.Label(self.root)
        self.image_label.pack(pady=10)
        
        # Current file label
        self.file_label = tk.Label(self.root, text="", font=('Arial', 10))
        self.file_label.pack(pady=5)
        
        # Button frame
        button_frame = tk.Frame(self.root)
        button_frame.pack(pady=10)
        
        # Buttons
        ttk.Button(button_frame, text="← Previous (Left Arrow)", 
                   command=self.prev_frame).grid(row=0, column=0, padx=5)
        ttk.Button(button_frame, text="Exploratory (E)", 
                   command=lambda: self.label_frame('exploratory')).grid(row=0, column=1, padx=5)
        ttk.Button(button_frame, text="Non-Exploratory (N)", 
                   command=lambda: self.label_frame('non-exploratory')).grid(row=0, column=2, padx=5)
        ttk.Button(button_frame, text="Uncertain (U)", 
                   command=lambda: self.label_frame('uncertain')).grid(row=0, column=3, padx=5)
        ttk.Button(button_frame, text="Next → (Right Arrow)", 
                   command=self.next_frame).grid(row=0, column=4, padx=5)
        
        # Save and change folder buttons
        bottom_frame = tk.Frame(self.root)
        bottom_frame.pack(pady=10)
        
        ttk.Button(bottom_frame, text="Save Labels (S)", 
                   command=self.save_labels).grid(row=0, column=0, padx=5)
        ttk.Button(bottom_frame, text="Change Folder", 
                   command=self.change_folder).grid(row=0, column=1, padx=5)
        
        # Instructions
        instructions = """
        Instructions:
        - Press 'E' or click button to label as EXPLORATORY (mouse sniffing/investigating object)
        - Press 'N' or click button to label as NON-EXPLORATORY (mouse not engaging with object)
        - Press 'U' if UNCERTAIN
        - Use Arrow Keys or buttons to navigate
        - Press 'S' to save progress
        """
        tk.Label(self.root, text=instructions, justify=tk.LEFT, font=('Arial', 9)).pack(pady=10)
        
        # Status label
        self.status_label = tk.Label(self.root, text="", font=('Arial', 10, 'bold'), fg='green')
        self.status_label.pack(pady=5)
        
        # Load first frame
        self.load_frame()
    
    def change_folder(self):
        """Save current work and switch to a different folder"""
        if self.labels:
            self.save_labels()
        
        # Clear current interface
        for widget in self.root.winfo_children():
            widget.destroy()
        
        # Show folder selection again
        self.frames_dir = None
        self.show_folder_selection()
        
    def load_frame(self):
        if self.current_idx >= len(self.image_files):
            self.status_label.config(text="All frames labeled!", fg='blue')
            return
        
        # Update progress
        labeled_count = len(self.labels)
        self.progress_label.config(
            text=f"Frame {self.current_idx + 1} / {len(self.image_files)} | Labeled: {labeled_count}"
        )
        
        # Load and display image
        img_path = os.path.join(self.frames_dir, self.image_files[self.current_idx])
        img = Image.open(img_path)
        
        # Resize if too large
        max_size = (1000, 700)
        img.thumbnail(max_size, Image.Resampling.LANCZOS)
        
        photo = ImageTk.PhotoImage(img)
        self.image_label.config(image=photo)
        self.image_label.image = photo  # Keep reference
        
        # Show filename and current label
        current_file = self.image_files[self.current_idx]
        current_label = self.labels.get(current_file, "UNLABELED")
        self.file_label.config(text=f"{current_file}\nCurrent label: {current_label}")
        
    def label_frame(self, label):
        current_file = self.image_files[self.current_idx]
        self.labels[current_file] = label
        self.status_label.config(text=f"Labeled as: {label}", fg='green')
        
        # Auto-advance to next frame
        self.next_frame()
        
    def next_frame(self):
        if self.current_idx < len(self.image_files) - 1:
            self.current_idx += 1
            self.load_frame()
            self.status_label.config(text="")
        else:
            self.status_label.config(text="End of frames. Press 'S' to save.", fg='orange')
    
    def prev_frame(self):
        if self.current_idx > 0:
            self.current_idx -= 1
            self.load_frame()
            self.status_label.config(text="")
    
    def save_labels(self):
        with open(self.output_file, 'w') as f:
            json.dump(self.labels, f, indent=2)
        self.status_label.config(text=f"Saved {len(self.labels)} labels!", fg='blue')
        print(f"Labels saved to {self.output_file}")
    
    def run(self):
        self.root.mainloop()

if __name__ == "__main__":
    labeler = FrameLabeler()
    labeler.run()