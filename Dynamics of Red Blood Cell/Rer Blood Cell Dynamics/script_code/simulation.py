import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def load_frames_from_file(file_path):
    frames = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
    current_frame = []
    for line in lines:
        line = line.strip()
        if line.startswith("# Frame"):
            if current_frame:
                frames.append(np.array(current_frame))
                current_frame = []
        elif line:
            current_frame.append(list(map(float, line.split())))
    if current_frame:
        frames.append(np.array(current_frame))
    return frames

def visualize_animation(frames, latitude_divisions, longitude_divisions, interval=100):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    def update(frame_index):
        ax.clear()
        frame = frames[frame_index]
        x, y, z = frame[:, 0], frame[:, 1], frame[:, 2]
        x = x.reshape((latitude_divisions + 1, longitude_divisions + 1))
        y = y.reshape((latitude_divisions + 1, longitude_divisions + 1))
        z = z.reshape((latitude_divisions + 1, longitude_divisions + 1))
        ax.plot_surface(x, y, z, color='red', edgecolor='none')
        ax.view_init(elev=20, azim= 2/3 * frame_index)
        ax.set_title(f"Frame {frame_index + 1}")
        ax.set_zlim(-1,1)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

    from matplotlib.animation import FuncAnimation
    anim = FuncAnimation(fig, update, frames=len(frames), interval=interval)
    anim.save('red_blood_cell_dynamics.gif', writer='pillow', fps=30)
    plt.show()

# Load and visualize
latitude_divisions = 200
longitude_divisions = 200
frames = load_frames_from_file("red_blood_cell_dynamics_data.txt")
visualize_animation(frames, latitude_divisions, longitude_divisions)
