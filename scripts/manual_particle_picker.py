import mrcfile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def pick_particles_with_input():
    # --- 1. User Inputs ---
    file_path = input("Enter the path to your MRC file: ")
    try:
        box_size = int(input("Enter the box size for picking (e.g., 64, 128, 256): "))
    except ValueError:
        print("Invalid input. Using default box size of 64.")
        box_size = 64

    # --- 2. Open and Prepare Data ---
    try:
        with mrcfile.open(file_path, mode='r') as mrc:
            data = mrc.data
            # If it's a 3D volume, we'll look at the middle slice
            if data.ndim == 3:
                img = data[data.shape[0] // 2]
            else:
                img = data
    except Exception as e:
        print(f"Error: Could not open file. {e}")
        return

    # --- 3. Interactive Plotting ---
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(img, cmap='gray', origin='lower')
    ax.set_title(f"Box Size: {box_size}px | Click to pick | Close to finish")

    picked_coords = []
    bottom_left_coords = []

    def on_click(event):
        if event.inaxes != ax:
            return
        
        x, y = event.xdata, event.ydata
        picked_coords.append((x, y))
        
        # Draw box centered on the click
        rect = Rectangle((x - box_size/2, y - box_size/2), 
                         box_size, box_size, 
                         linewidth=1, edgecolor='yellow', facecolor='none')
        ax.add_patch(rect)
        plt.draw()
        #print(f"Selected particle at X={x:.1f}, Y={y:.1f}")

        # Calculate bottom-left corner and convert to integer
        # (x_click - half_box, y_click - half_box)
        x_bl = int(round(event.xdata - box_size / 2))
        y_bl = int(round(event.ydata - box_size / 2))

        bottom_left_coords.append((x_bl, y_bl))

        # Visual confirmation using the calculated integer coordinates
        rect = Rectangle((x_bl, y_bl), box_size, box_size,
                         linewidth=1, edgecolor='cyan', facecolor='none')
        ax.add_patch(rect)
        plt.draw()
        #print(f"Bottom-Left Corner Picked: X={x_bl}, Y={y_bl}")


    fig.canvas.mpl_connect('button_press_event', on_click)
    plt.show()

    # --- 4. Final Output ---
    print(f"\nFinished! Total particles picked: {len(picked_coords)}")
    for i, (px, py) in enumerate(picked_coords):
        print(f"Particle {i+1}: X={px:.2f}, Y={py:.2f}")

    # Save boxes to file
    output_file = "boxfile.box"
    with open(output_file, "w") as f:
        for x, y in bottom_left_coords:
            f.write(f"{x}\t{y}\t{box_size}\t{box_size}\n")

if __name__ == "__main__":
    pick_particles_with_input()
